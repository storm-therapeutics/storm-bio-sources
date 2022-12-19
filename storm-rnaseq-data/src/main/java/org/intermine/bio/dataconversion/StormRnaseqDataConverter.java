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
                processRNASeqComparison(deSeq2File, meta, comparison.treatment, comparison.control);
            }

            // Process feature counts
            LOG.debug("StormRnaseqDataConverter [process] - processing feature counts: " + meta.experimentShortName);
            String fileName = "salmon.merged.gene_counts.tsv";
            File geneCountsFile = filesInDir.get(fileName);
            if (geneCountsFile != null) {
                processRNASeqGeneCounts(geneCountsFile, meta);
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
        // example formats (experiments: "Lexogen_Dec19", "Lexogen_May20"):
        // ensembl	entrez	symbol	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        // Ensembl	Entrez	Gene	Description	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        Map<String, Integer> indexes = new HashMap<String, Integer>();
        // two possible formats (legacy and current) for initial columns:
        List<String> columns;
        if (header.length == 9) { // legacy format
            columns = List.of("ensembl", "entrez", "symbol");
        }
        else if (header.length == 10) { // current format
            columns = List.of("Ensembl", "Entrez", "Gene", "Description");
        }
        else {
            throw new RuntimeException("Unexpected number of columns in DESeq2 file");
        }
        // check that column names match expectations:
        for (int i = 0; i < columns.size(); i++) {
            if (!header[i].equals(columns.get(i))) {
                throw new RuntimeException("Unexpected column name in DESeq2 file: '" +
                                           header[i] + "' vs. '" + columns.get(i) + "'");
            }
        }
        // check remaining columns (same for both formats):
        String[] moreColumns = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"};
        for (int i = header.length - moreColumns.length; i < header.length; i++) {
            if (!header[i].equals(moreColumns[i])) {
                throw new RuntimeException("Unexpected column name in DESeq2 file: '" +
                                           header[i] + "' vs. '" + moreColumns[i] + "'");
            }
            indexes.put(header[i], i);
        }

        return indexes;
    }


    private void processRNASeqComparison(File DESeq2File, StormOmicsMetadata meta, String treatmentName, String controlName) throws ObjectStoreException, IOException {
        String experimentShortName = meta.experimentShortName;
        String fileName = DESeq2File.getName();
        String fileAbsPath = DESeq2File.getAbsolutePath();

        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] header = lineIter.next();
        Map<String, Integer> columnIndexes = getColumnIndexes(header);

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            Item gene = getGene(line[0], line[1], line[2]);
            if (gene == null) // could not find matching gene
                continue;

            Item integratedItem = createItem("StormRNASeqComparison");
            // integratedItem.setAttribute("geneEnsemblId", ensemblId);
            integratedItem.setReference("gene", gene);
            integratedItem.setReference("experiment", meta.experiment);
            integratedItem.setReference("control", meta.conditions.get(controlName));
            integratedItem.setReference("treatment", meta.conditions.get(treatmentName));

            for (Map.Entry<String, Integer> entry : columnIndexes.entrySet()) {
                String key = entry.getKey();
                int index = entry.getValue();
                String value = line[index];
                if (!value.equals("NA")) {
                    integratedItem.setAttribute("key", value);
                }
            }

            store(integratedItem);
        }
    }


    private void processRNASeqGeneCounts(File geneCountsFile, StormOmicsMetadata meta) throws ObjectStoreException, IOException {
        String experimentShortName = meta.experimentShortName;
        String fileAbsPath = geneCountsFile.getAbsolutePath();
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] firstLine = lineIter.next();

        ArrayList<String> runs = new ArrayList<String>();
        for (int i = 2; i < firstLine.length; i++) {
            String run = firstLine[i];
            runs.add(run);
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            Item gene = getGene(line[0], null, line[1]); // no Entrez ID in this file
            if (gene == null) // could not find matching gene
                continue;

            try {
                for (int i = 2; i < line.length; i++) {
                    String count = line[i];
                    String runForThisItem = runs.get(i-2);
                    Item integratedItem = createItem("StormRNASeqFeatureCounts");
                    integratedItem.setReference("experiment", meta.experiment);
                    integratedItem.setReference("gene", gene);
                    // if (!geneEnsemblId.isEmpty()) {
                    //     integratedItem.setAttribute("geneEnsemblId", geneEnsemblId);
                    // }
                    integratedItem.setAttribute("run", runForThisItem);
                    // TODO: exclude zero counts?
                    integratedItem.setAttribute("count", count);
                    store(integratedItem);
                }
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqGeneCount with gene: " + gene + " - " + e.getMessage());
                continue;
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


    private Item getGene(String ensemblId, String entrezId, String symbol) throws ObjectStoreException {
        ensemblId = ensemblId.split("\\.")[0]; // remove version number (if any)
        // have we looked up this ID before?
        if (geneItems.containsKey(ensemblId)) {
            return geneItems.get(ensemblId);
        }
        // look up ID, create new map entry:
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
            gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", resolved.iterator().next());
            store(gene);
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
