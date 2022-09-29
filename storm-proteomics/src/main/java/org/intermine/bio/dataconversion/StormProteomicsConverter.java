package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2022 STORM Therapeutics Limited
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileReader;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;
import org.intermine.util.FormattedTextParser;

/**
 *
 * @author Hendrik Weisser
 */
public class StormProteomicsConverter extends BioDirectoryConverter
{
    private static final String DATASET_TITLE = "STORM Proteomics Data";
    private static final String DATA_SOURCE_NAME = "STORM Therapeutics";

    // expected input file names:
    private static final String MZTAB_NAME = "out.mzTab";
    private static final String MSSTATS_NAME = "msstats_comparisons.csv"; // actually tab-separated

    private static final Logger LOG = Logger.getLogger(StormProteomicsConverter.class);

    Item experiment;

    /**
     * Constructor
     */
    public StormProteomicsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        // loop through subdirectories (corresponding to experiments):
        for (File subDir : dataDir.listFiles()) {
            if (subDir.isDirectory()) {
                File mztabFile = new File(subDir, MZTAB_NAME);
                if (mztabFile.exists()) {
                    LOG.info("Processing mzTab file in " + subDir.getName());
                    experiment = createItem("StormProteomicsExperiment");
                    experiment.setAttribute("shortName", subDir.getName());
                    store(experiment);
                    processMzTab(mztabFile);
                    File msstatsFile = new File(subDir, MSSTATS_NAME);
                    if (msstatsFile.exists()) {
                        LOG.info("Processing MSstats file in " + subDir.getName());
                        // processMSstats(msstatsFile);
                    }
                }
            }
        }
    }


    private Item makeProtein(String accession) throws Exception {
        // parse accession:
        String[] parts = accession.split("\\|");
        if ((parts.length != 3) || (!parts[0].equals("sp") && !parts[0].equals("tr"))) {
            throw new RuntimeException("Unexpected accession format in mzTab file: " + accession);
        }
        Item protein = createItem("Protein");
        protein.setAttribute("primaryIdentifier", parts[2]);
        protein.setAttribute("primaryAccession", parts[1]);
        store(protein);
        return protein;
    }


    private void processMzTab(File file) throws Exception {
        FileReader reader = new FileReader(file);
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        int nCol = 0;
        ArrayList<Item> samples = new ArrayList<Item>();
        Map<String, Integer> proteinHeader = new HashMap<String, Integer>();
        // parse metadata section and protein section header:
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if ((line.length == 0) || (line[0].equals("COM"))) // skip empty lines and comments
                continue;
            if (line[0].equals("MTD")) { // metadata section
                if (line[1].equals("ms_run[" + (samples.size() + 1) + "]-location")) {
                    // remove file extension to get sample name:
                    String sampleName = line[2].substring(0, line[2].lastIndexOf('.'));
                    if (sampleName.startsWith("file://")) {
                        sampleName = sampleName.substring(7);
                    }
                    Item sample = createItem("StormProteomicsSample");
                    sample.setAttribute("name", sampleName);
                    sample.setReference("experiment", experiment);
                    store(sample);
                    samples.add(sample);
                }
            }
            else if (line[0].equals("PRH")) { // protein section header
                for (int i = 1; i < line.length; i++) {
                    proteinHeader.put(line[i], i);
                }
                nCol = line.length;
                break;
            }
        }
        if (samples.isEmpty()) {
            throw new RuntimeException("No sample information ('MTD ms_run[...]-location') found in mzTab file");
        }
        if (proteinHeader.isEmpty() || (nCol == 0)) {
            throw new RuntimeException("No protein section found in mzTab file");
        }
        // look up column indexes in protein header:
        int resultTypeIndex = proteinHeader.getOrDefault("opt_global_result_type", -1);
        int[] studyVariableIndexes = new int[samples.size()];
        for (int i = 0; i < samples.size(); i++) {
            String colName = "protein_abundance_study_variable[" + (i + 1) + "]";
            studyVariableIndexes[i] = proteinHeader.get(colName);
        }
        // parse protein section:
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if ((line.length == 0) || (line[0].equals("COM"))) // skip empty lines and comments
                continue;
            if (line[0].equals("PRT")) { // protein section
                if (line.length != nCol) {
                    throw new RuntimeException("Unexpected number of columns in mzTab protein section: " + line.length + " (PRT), " + nCol + " (PRH)");
                }
                if ((resultTypeIndex >= 0) && (line[resultTypeIndex].equals("protein_details"))) {
                    continue; // no quantification - skip
                }
                Item protein = makeProtein(line[1]);
                for (int i = 0; i < samples.size(); i++) {
                    String value = line[studyVariableIndexes[i]];
                    // treat "0" as missing value (see https://github.com/OpenMS/OpenMS/issues/6363):
                    if (!value.equals("null") && !value.equals("0.0")) {
                        Item abundance = createItem("ProteinLFQAbundance");
                        abundance.setAttribute("abundance", value);
                        abundance.setReference("sample", samples.get(i));
                        abundance.setReference("protein", protein);
                        store(abundance);
                    }
                }
            }
            else { // protein section finished
                break;
            }
        }
    }


    private void processMSstats(File file) throws Exception {
        FileReader reader = new FileReader(file);
        Iterator<String[]> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
    }
}
