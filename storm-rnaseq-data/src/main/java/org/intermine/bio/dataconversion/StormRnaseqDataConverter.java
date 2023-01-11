package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2021-2022 STORM Therapeutics Limited
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Read results from RNA-seq data analysis
 *
 * @author Adrian Bazaga, Hendrik Weisser
 */
public class StormRnaseqDataConverter extends BioDirectoryConverter
{
    private static final String DATASET_TITLE = "STORM RNA-Seq Data";
    private static final String DATA_SOURCE_NAME = "STORM Therapeutics";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    // mapping: Ensembl ID -> gene (InterMine Item)
    private Map<String, Item> geneItems = new HashMap<String, Item>();
    // mapping: NCBI (primary) ID -> Ensembl (secondary) ID
    private Map<String, String> geneIDs = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormRnaseqDataConverter.class);

    public StormRnaseqDataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {
        Map<String, File> directories = readDirectoriesInDir(dataDir);

        // Get all JSON config files in the directory and process one by one
        File[] jsonFiles = dataDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".json");
            }
        });

        for (File jsonFile : jsonFiles) {
            // read experiment metadata:
            StormOmicsMetadata meta = new StormOmicsMetadata(this);
            meta.processJSONFile(jsonFile);

            LOG.info("Processing RNA-seq data for experiment: " + meta.experimentShortName);
            Item experiment = createItem("StormRNASeqExperiment");
            experiment.setReference("metadata", meta.experiment);

            // find matching results:
            File experimentDir = new File(dataDir.getAbsolutePath(), meta.experimentShortName);
            Map<String, File> filesInDir = readFilesInDir(experimentDir);

            // Process differential expression results
            LOG.debug("StormRnaseqDataConverter [process] - processing comparisons: " + meta.experimentShortName);
            for (StormOmicsMetadata.ConditionsPair comparison : meta.comparisons) {
                String fileName = comparison.treatment + "_vs_" + comparison.control + "_DESeq2.tsv";
                File deSeq2File = filesInDir.get(fileName);
                if (deSeq2File == null) {
                    LOG.info("Failed to find DESeq2 file: " + fileName);
                    continue;
                }
                processRNASeqComparison(deSeq2File, experiment, meta, comparison.treatment, comparison.control);
            }

            // Process feature counts
            LOG.debug("StormRnaseqDataConverter [process] - processing feature counts: " + meta.experimentShortName);
            String fileName = "salmon.merged.gene_counts.tsv";
            File geneCountsFile = filesInDir.get(fileName);
            if (geneCountsFile != null) {
                processRNASeqGeneCounts(geneCountsFile, experiment, meta);
            } else {
                LOG.info("Failed to find feature counts file: " + fileName);
                // continue; // abort or store experiment without counts?
            }

            // Store the experiment
            LOG.info("Storing experiment: " + meta.experimentShortName);
            try {
                store(experiment);
            } catch (Exception e) {
                throw new RuntimeException("Error storing StormRNASeqExperiment ", e);
            }
        }
    }

    private Map<String, Integer> getColumnIndexes(String[] header) {
        // example formats:
        // ensembl	entrez	symbol	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        // Ensembl	Entrez	Gene	Description	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        // Gene	Ensembl	Entrez	Description	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        ArrayList<String> headerLower = new ArrayList<String>(header.length);
        for (String entry : header) {
            headerLower.add(entry.toLowerCase());
        }
        Map<String, Integer> indexes = new HashMap<String, Integer>();
        // special case - "symbol" or "Gene" for gene symbol:
        int symbolIndex = -1;
        if (header.length == 9) {
            symbolIndex = headerLower.indexOf("symbol");
        }
        else if (header.length == 10) {
            symbolIndex = headerLower.indexOf("gene");
        }
        else {
            throw new RuntimeException("Unexpected number of columns in DESeq2 file");
        }
        if (symbolIndex == -1) {
            throw new RuntimeException("Column 'Gene'/'symbol' not found in DESeq2 file");
        }
        indexes.put("symbol", symbolIndex);
        // check remaining columns (ignore optional "Description"):
        String[] columns = {"ensembl", "entrez", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"};
        for (String column : columns) {
            int index = headerLower.indexOf(column.toLowerCase());
            if (index == -1) {
                throw new RuntimeException("Column '" + column + "' not found in DESeq2 file");
            }
            indexes.put(column, index);
        }
        return indexes;
    }


    private void processRNASeqComparison(File DESeq2File, Item experiment, StormOmicsMetadata meta, String treatmentName, String controlName) throws ObjectStoreException, IOException {
        String fileName = DESeq2File.getName();
        String fileAbsPath = DESeq2File.getAbsolutePath();

        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] header = lineIter.next();

        // what information is in which column?
        Map<String, Integer> columnIndexes = getColumnIndexes(header);
        int symbolIndex = columnIndexes.get("symbol");
        int ensemblIndex = columnIndexes.get("ensembl");
        int entrezIndex = columnIndexes.get("entrez");
        String[] valueColumns = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"};
        ArrayList<Integer> valueIndexes = new ArrayList<Integer>(valueColumns.length);
        for (String column : valueColumns) {
            valueIndexes.add(columnIndexes.get(column));
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                LOG.error("Unexpected number of items per line: " + line.length + " vs. " + header.length);
                continue;
            }

            Item gene = getGene(line[ensemblIndex], line[entrezIndex], line[symbolIndex]);
            if (gene == null) // could not find matching gene
                continue;

            Item integratedItem = createItem("StormRNASeqComparison");
            // integratedItem.setAttribute("geneEnsemblId", ensemblId);
            integratedItem.setReference("gene", gene);
            integratedItem.setReference("experiment", experiment);
            integratedItem.setReference("control", meta.conditions.get(controlName));
            integratedItem.setReference("treatment", meta.conditions.get(treatmentName));

            for (int i = 0; i < valueColumns.length; i++) {
                String key = valueColumns[i];
                int index = valueIndexes.get(i);
                String value = line[index];
                if (!value.equals("NA")) {
                    integratedItem.setAttribute(key, value);
                }
            }

            store(integratedItem);
        }
    }


    private void processRNASeqGeneCounts(File geneCountsFile, Item experiment, StormOmicsMetadata meta) throws ObjectStoreException, IOException {
        String fileAbsPath = geneCountsFile.getAbsolutePath();
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] header = lineIter.next();
        // two possible formats:
        int geneColumns; // number of columns with gene information
        if ((header[0].equals("Geneid") || header[0].equals("gene_id")) && header[1].equals("gene_name")) {
            geneColumns = 2;
        }
        else if (header[0].equals("gene_id")) {
            geneColumns = 1;
        }
        else {
            throw new RuntimeException("Unexpected header in gene counts file");
        }

        // condition names may appear in altered form:
        Map<String, String> conditionNames = new HashMap<String, String>();
        for (String conditionName : meta.conditions.keySet()) {
            conditionNames.put(conditionName.replaceAll("-", "."), conditionName);
        }

        // create a template Item for each data column:
        Map<Integer, Item> columnItems = new HashMap<Integer, Item>();
        Map<String, String> bioReplicates = null; // reverse mapping of 'meta.bioReplicates' (created if needed)
        for (int i = geneColumns; i < header.length; i++) {
            Item columnItem = createItem("StormRNASeqFeatureCount");
            columnItem.setReference("experiment", experiment);
            // column names can be based on samples (legacy) or conditions (current):
            String name = header[i];
            if (name.matches(".*_R\\d+")) {
                // current format - biological replicates: "[condition]_R1", "[condition]_R2" etc.
                int index = name.lastIndexOf('_');
                String conditionName = name.substring(0, index);
                String replicate = name.substring(index + 2);
                // replace "altered" condition name with original one, if applicable:
                conditionName = conditionNames.getOrDefault(conditionName, conditionName);
                Item condition = meta.conditions.get(conditionName);
                if (condition == null) {
                    LOG.error("Could not find condition corresonding to column in gene counts file: " +
                              conditionName);
                    continue;
                }
                columnItem.setReference("condition", condition);
                columnItem.setAttribute("replicate", replicate);
                ArrayList<String> sampleNames = meta.bioReplicates.get(name);
                if (sampleNames != null) {
                    for (String sampleName : sampleNames) {
                        columnItem.addToCollection("samples", meta.samples.get(sampleName));
                    }
                }
                else {
                    LOG.error("Could not find corresponding samples for column in gene counts file: " + name);
                }
                columnItems.put(i, columnItem);
            }
            else { // legacy format: sample names
                Item sample = meta.samples.get(name);
                if ((sample == null) && name.startsWith("X")) { // try look-up again without "X" prefix
                    name = name.substring(1);
                    sample = meta.samples.get(name);
                }
                if (sample != null) {
                    if (bioReplicates == null) { // create mapping: sample -> biological replicate
                        bioReplicates = new HashMap<String, String>();
                        for (Map.Entry<String, ArrayList<String>> entry : meta.bioReplicates.entrySet()) {
                            String bioRep = entry.getKey();
                            int index = bioRep.lastIndexOf('_');
                            String rep = bioRep.substring(index + 2);
                            for (String sampleName : entry.getValue()) {
                                bioReplicates.put(sampleName, rep);
                            }
                        }
                    }
                    columnItem.addToCollection("samples", sample);
                    columnItem.setReference("condition", sample.getReference("condition").getRefId());
                    columnItem.setAttribute("replicate", bioReplicates.get(name));
                    columnItems.put(i, columnItem);
                }
                else {
                    LOG.error("Could not find sample corresponding to column in gene counts file: " + name);
                }
            }
        }
        if (columnItems.isEmpty()) {
            LOG.error("Could not map any conditions/samples in gene counts file - aborting");
            return; // TODO: throw an exception instead?
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                LOG.error("Unexpected number of items per line: " + line.length + " vs. " + header.length);
                continue;
            }
            // find gene:
            String symbol = null;
            if (geneColumns == 2)
                symbol = line[1];
            Item gene = getGene(line[0], null, symbol); // no Entrez ID in this file
            if (gene == null) // could not find matching gene
                continue;

            for (Map.Entry<Integer, Item> entry : columnItems.entrySet()) {
                String value = line[entry.getKey()];
                if (value.equals("0") || value.equals("0.0"))
                    continue; // don't store zero counts
                // create new Items to get new identifiers (avoid "duplicate identifier" errors);
                // alternative solution using 'newId()' to assign new id. gives strange SQL error...
                Item dummy = createItem("StormRNASeqFeatureCount");
                Item integratedItem = entry.getValue();
                integratedItem.setReference("gene", gene);
                // if (!geneEnsemblId.isEmpty()) {
                //     integratedItem.setAttribute("geneEnsemblId", geneEnsemblId);
                // }
                integratedItem.setAttribute("count", value);
                store(integratedItem);
                integratedItem.setIdentifier(dummy.getIdentifier());
                // integratedItem.setIdentifier(newId("StormRNASeqFeatureCount")); // update id. for next round
            }
        }
    }


    private Map<String, File> readDirectoriesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            if(file.isDirectory()) {
                files.put(file.getName(), file);
            }
        }
        return files;
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }


    /**
     * Look up a gene by its Ensembl ID (using Entrez ID and/or gene symbol to resolve ambiguities)
     *
     * Return the Item for the gene, or 'null' if not found.
     * Create a new entry (possibly 'null') in 'geneItems' the first time an Ensembl ID is looked up.
     */
    private Item getGene(String ensemblId, String entrezId, String symbol) throws ObjectStoreException {
        ensemblId = ensemblId.split("\\.")[0]; // remove version number (if any)
        // have we looked up this ID before?
        if (geneItems.containsKey(ensemblId)) {
            return geneItems.get(ensemblId);
        }
        // look up primary ID, create new map entry:
        Set<String> resolved = rslv.resolveId(TAXON_ID, "gene", ensemblId);
        boolean multipleMatches = resolved.size() > 1;
        if (multipleMatches) { // use additional information to choose the gene
            if ((entrezId != null) && !entrezId.isEmpty() && !entrezId.equals("NA")) { // try to use Entrez ID
                if (resolved.contains(entrezId)) {
                    resolved = Collections.singleton(entrezId);
                }
            }
            else if ((symbol != null) && !symbol.isEmpty() && !symbol.equals("NA")) { // try to use gene symbol
                symbol = symbol.split("\\.")[0]; // remove version number (if any)
                Set<String> resolvedBySymbol = rslv.resolveId(TAXON_ID, "gene", symbol);
                resolved.retainAll(resolvedBySymbol);
            }
        }
        Item gene = null;
        if (resolved.size() == 1) { // success
            String primaryId = resolved.iterator().next();
            // problem: multiple Ensembl genes can match to the same Entrez/NCBI ID,
            // causing errors ("duplicate objects") during integration;
            // make sure not to store gene Items with the same ID twice:
            String existing = geneIDs.putIfAbsent(primaryId, ensemblId);
            if (existing == null) { // no gene resolved to this primary ID yet
                gene = createItem("Gene");
                gene.setAttribute("primaryIdentifier", primaryId);
                store(gene);
            }
            else { // conflict
                LOG.error("Multiple matches to gene primary ID " + primaryId +
                          ": " + existing + " (kept), " + ensemblId + " (skipped)");
            }
        }
        else if (multipleMatches) {
            LOG.error("Failed to resolve multiple gene matches for Ensembl ID " + ensemblId);
        }
        else {
            LOG.error("No gene found for Ensembl ID " + ensemblId);
        }
        geneItems.put(ensemblId, gene);
        return gene;
    }
}
