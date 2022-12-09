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
    private static final String DATA_SOURCE_NAME = "STORM RNA-Seq Data";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();

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
            // find matching results:
            File experimentDir = new File(dataDir.getAbsolutePath(), meta.experimentShortName);
            Map<String, File> filesInDir = readFilesInDir(experimentDir);

            // Process differential expression results
            for (StormOmicsMetadata.ConditionsPair comparison : meta.comparisons) {
                String fileName = comparison.treatment + "_vs_" + comparison.control + "_DESeq2.tsv";

                if (filesInDir.get(fileName) != null) {
                    File deSeq2File = filesInDir.get(fileName);
                    processRNASeqExperimentComparison(deSeq2File, meta, comparison.treatment, comparison.control);
                } else {
                    LOG.info("Failed to find DESeq2 file: " + fileName);
                    continue;
                    //throw new RuntimeException("Failed to find DESeq2 file: " + fileName);
                }
            }

            // Process feature counts
            LOG.info("StormRnaseqDataConverter [processConfigFile] - Processing 5: " + meta.experimentShortName);
            String fileName = "salmon.merged.gene_counts.tsv";
            File geneCountsFile = filesInDir.get(fileName);
            if (geneCountsFile != null) {
                processRNASeqExperimentGeneCount(geneCountsFile, meta);
            } else {
                LOG.info("Failed to find counts file: " + fileName);
                continue;
                //throw new RuntimeException("Failed to find DESeq2 file: " + geneCountsFile);
            }

            // Store the experiment
            LOG.info("Storing experiment: " + meta.experimentShortName);
            try {
                store(meta.experiment);
            } catch (Exception e) {
                throw new RuntimeException("Error storing StormRNASeqExperiment ", e);
            }
        }
    }

    private Map<String, Integer> getColumnIndexes(String[] header, String fileType) {
        Map<String, Integer> indexes = new HashMap<String, Integer>();

        if(fileType == "DESEQ2") {
            //20: Ensembl	Entrez	Gene	Description	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
            //19: ensembl	entrez	symbol	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
            for(int i = 0; i < header.length; i++) {
                String column = header[i];
                if(column.toLowerCase().contains("ensembl")) {
                    indexes.put("gene", i);
                    continue;
                } else if(column.toLowerCase().contains("basemean")) {
                    indexes.put("baseMean", i);
                    continue;
                } else if(column.toLowerCase().contains("log2foldchange")) {
                    indexes.put("log2FoldChange", i);
                    continue;
                } else if(column.toLowerCase().contains("lfcse")) {
                    indexes.put("lfcSE", i);
                    continue;
                } else if(column.toLowerCase().contains("stat")) {
                    indexes.put("stat", i);
                    continue;
                } else if(column.toLowerCase().contains("pvalue")) {
                    indexes.put("pvalue", i);
                    continue;
                } else if(column.toLowerCase().contains("padj")) {
                    indexes.put("padj", i);
                    continue;
                }
            }
        } else if(fileType == "SampleInfo") {
            //20: Run,Sample,Condition,Celline,IFN_gamma,Compound,Concentration µM,Replicate,Timepoint,Lexogen_ID,Filename
            //19: Run,Condition,Celline,IFN gamma,Compound,Concentration,Replicate,Timepoint,Sample number,Sample,File name
            for(int i = 0; i < header.length; i++) {
                String column = header[i];
                if(column.toLowerCase().contains("run")) {
                    indexes.put("run", i);
                    continue;
                } else if(column.toLowerCase().equals("sample")) {
                    indexes.put("sample", i);
                    continue;
                } else if(column.toLowerCase().contains("condition") || column.toLowerCase().contains("treatment")) {
                    indexes.put("treatment", i);
                    continue;
                } else if(column.toLowerCase().contains("cellline") || column.toLowerCase().contains("celline")) {
                    indexes.put("cellLine", i);
                    continue;
                } else if(column.toLowerCase().contains("ifn_gamma") || column.toLowerCase().contains("ifn gamma")) {
                    indexes.put("IFN_gamma", i);
                    continue;
                } else if(column.toLowerCase().contains("compound")) {
                    indexes.put("compound", i);
                    continue;
                } else if(column.toLowerCase().contains("concentration")) {
                    indexes.put("concentration", i);
                    continue;
                } else if(column.toLowerCase().contains("replicate")) {
                    indexes.put("replicate", i);
                    continue;
                } else if(column.toLowerCase().contains("timepoint")) {
                    indexes.put("timepoint", i);
                    continue;
                }
            }
        } else {
            throw new BuildException("Unknown fileType for getColumnIndexes: " + fileType);
        }

        return indexes;
    }




    private void processRNASeqExperimentComparison(File DESeq2File, StormOmicsMetadata meta, String treatmentName, String controlName) throws ObjectStoreException, IOException {
        String experimentShortName = meta.experimentShortName;
        String fileName = DESeq2File.getName();

        if (fileName.endsWith("_DESeq2.tsv")) {
            String fileAbsPath = DESeq2File.getAbsolutePath();

            Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
            String[] firstLine = (String[]) lineIter.next();
            Map<String, Integer> columnIndexes = getColumnIndexes(firstLine, "DESEQ2");

            while (lineIter.hasNext()) {
                String[] line = (String[]) lineIter.next();

                String gene = line[columnIndexes.get("gene").intValue()].split("\\.")[0];
                String geneEnsemblId = line[0].split("\\.")[0];
                String baseMean = line[columnIndexes.get("baseMean").intValue()];
                String log2FoldChange = line[columnIndexes.get("log2FoldChange").intValue()];
                String lfcSE = line[columnIndexes.get("lfcSE").intValue()];
                String stat = line[columnIndexes.get("stat").intValue()];
                String pvalue = line[columnIndexes.get("pvalue").intValue()];
                String padj = line[columnIndexes.get("padj").intValue()];

                Item integratedItem = createItem("RNASeqExperimentComparison");

                if (!gene.isEmpty()) {
                    if(unresolvableGenes.get(gene) != null) {
                        continue;
                    }
                    String geneId = getGeneId(gene);
                    if(geneId == null) {
                        continue;
                    }
                    integratedItem.setReference("gene", geneId);
                } else {
                    throw new BuildException("[processRNASeqExperimentComparison] gene was empty: " + fileName);
                }

                if (!geneEnsemblId.isEmpty()) {
                    integratedItem.setAttribute("geneEnsemblId", geneEnsemblId);
                }

                if (meta.conditions.containsKey(controlName)) {
                    integratedItem.setReference("control", meta.conditions.get(controlName));
                } else {
                    continue;
                }
                if (meta.conditions.containsKey(treatmentName)) {
                    integratedItem.setReference("treatment", meta.conditions.get(treatmentName));
                } else {
                    continue;
                }

                if (!StringUtils.isEmpty(baseMean) && isDouble(baseMean)) {
                    integratedItem.setAttribute("baseMean", baseMean);
                }
                if (!StringUtils.isEmpty(log2FoldChange) && isDouble(log2FoldChange)) {
                    integratedItem.setAttribute("log2FoldChange", log2FoldChange);
                }
                if (!StringUtils.isEmpty(lfcSE) && isDouble(lfcSE)) {
                    integratedItem.setAttribute("lfcSE", lfcSE);
                }
                if (!StringUtils.isEmpty(stat) && isDouble(stat)) {
                    integratedItem.setAttribute("stat", stat);
                }
                if (!StringUtils.isEmpty(pvalue) && isDouble(pvalue)) {
                    integratedItem.setAttribute("pvalue", pvalue);
                }
                if (!StringUtils.isEmpty(padj) && isDouble(padj)) {
                    integratedItem.setAttribute("padj", padj);
                }

                integratedItem.setReference("experiment", meta.experiment);
                store(integratedItem);
            }
        }
    }


    private void processRNASeqExperimentGeneCount(File geneCountsFile, StormOmicsMetadata meta) throws ObjectStoreException, IOException {
        String experimentShortName = meta.experimentShortName;
        String fileAbsPath = geneCountsFile.getAbsolutePath();
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] firstLine = (String[]) lineIter.next();

        ArrayList<String> runs = new ArrayList<String>();
        for (int i = 2; i < firstLine.length; i++) {
            String run = firstLine[i];
            runs.add(run);
        }

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            String gene = line[1];
            String geneEnsemblId = line[0].split("\\.")[0];
            try {
                for (int i = 2; i < line.length; i++) {
                    String count = line[i];
                    String runForThisItem = runs.get(i-2);
                    Item integratedItem = createItem("RNASeqExperimentFeatureCounts");
                    if (!gene.isEmpty()) {
                        if (unresolvableGenes.get(gene) != null) {
                            continue;
                        }
                        String geneId = getGeneId(gene);
                        if (geneId == null) {
                            continue;
                        }
                        integratedItem.setReference("gene", geneId);
                    } else {
                        continue;
                    }

                    if (!geneEnsemblId.isEmpty()) {
                        integratedItem.setAttribute("geneEnsemblId", geneEnsemblId);
                    }
                    if (!StringUtils.isEmpty(runForThisItem)) {
                        integratedItem.setAttribute("run", runForThisItem);
                    } else {
                        continue;
                    }
                    if (!StringUtils.isEmpty(count) && isDouble(count)) {
                        integratedItem.setAttribute("count", count);
                    }

                    integratedItem.setReference("experiment", meta.experiment);
                    store(integratedItem);
                }
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentGeneCount with gene: " + gene + " - " + e.getMessage());
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

    private String getGeneId(String identifier) throws ObjectStoreException {
        String geneId = null;
        try {
            String resolvedIdentifier = resolveGene(identifier);
            if(resolvedIdentifier != null) {
                geneId = genes.get(resolvedIdentifier);
                if (geneId == null) {
                    Item gene = createItem("Gene");
                    gene.setAttribute("primaryIdentifier", resolvedIdentifier);
                    store(gene);
                    geneId = gene.getIdentifier();
                    genes.put(resolvedIdentifier, geneId);
                }
                return geneId;
            } else {
                return resolvedIdentifier;
            }
        } catch (Exception e) {
            LOG.info("getGeneId: failed to resolve gene: " + identifier);
            return null;
        }
    }

    private String resolveGene(String identifier) {
        String id = null;

        if(resolvedGenes.get(identifier) != null) {
            id = resolvedGenes.get(identifier);
        } else {
            if (rslv != null && rslv.hasTaxon(TAXON_ID)) {
                int resCount = rslv.countResolutions(TAXON_ID, identifier);
                if (resCount != 1) {
                    unresolvableGenes.put(identifier, identifier);
                    return null;
                }
                id = rslv.resolveId(TAXON_ID, identifier).iterator().next();
                resolvedGenes.put(identifier, id);
            }
        }
        return id;
    }

    private boolean isDouble(String str) {
        try {
            double x = Double.parseDouble(str);
            //if (x == (int) x)
            //    return false;
            return true;
        }
        catch(NumberFormatException e) {
            return false;
        }
    }
}
