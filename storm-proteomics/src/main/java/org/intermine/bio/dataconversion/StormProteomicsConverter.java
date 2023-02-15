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

import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.AbstractMap.SimpleEntry; // pair class

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
    StormOmicsMetadata metadata;
    Map<String, Item> samples = new HashMap<String, Item>();
    Map<String, Item> conditions = new HashMap<String, Item>();
    Map<String, Item> proteinGroups = new HashMap<String, Item>();

    /**
     * Constructor
     */
    public StormProteomicsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * Process directories containing proteomics results
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        // loop through subdirectories (corresponding to experiments):
        for (File subDir : dataDir.listFiles()) {
            if (subDir.isDirectory()) {
                File mztabFile = new File(subDir, MZTAB_NAME);
                if (!mztabFile.exists()) {
                    LOG.info("No results (mzTab) file found in directory '" + subDir.getName() + "' - skipping");
                    continue;
                }
                // read full metadata from JSON file, or sample/condition information from TSV:
                File jsonFile = new File(subDir, subDir.getName() + ".json");
                if (!jsonFile.exists()) {
                    LOG.info("No metadata (JSON) file found in directory '" + subDir.getName() + "' - skipping");
                    continue;
                }
                LOG.info("Processing proteomics results in directory '" + subDir.getName() + "'");
                samples.clear();
                conditions.clear();
                proteinGroups.clear();
                LOG.info("Processing metadata (JSON) file");
                metadata = new StormOmicsMetadata(this);
                metadata.processJSONFile(jsonFile);
                experiment = createItem("StormProteomicsExperiment");
                experiment.setReference("metadata", metadata.experiment);
                // raw protein abundances:
                LOG.info("Processing mzTab file");
                processMzTab(mztabFile);
                // differential protein expression (optional):
                File msstatsFile = new File(subDir, MSSTATS_NAME);
                if (msstatsFile.exists()) {
                    LOG.info("Processing MSstats file");
                    processMSstats(msstatsFile);
                }
                store(experiment);
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


    private Item makeProteinGroup(String accessions, String sep, String score) throws Exception {
        String[] parts = accessions.replace(" ", "").split(sep);
        // generate unique key for the group:
        Arrays.sort(parts);
        String sortedAccessions = String.join(";", parts);
        if (proteinGroups.containsKey(sortedAccessions)) {
            throw new RuntimeException("Duplicate protein group: " + sortedAccessions);
        }
        Item group = createItem("ProteomicsProteinGroup");
        group.setAttribute("accessions", sortedAccessions);
        group.setReference("experiment", experiment);
        if (score != null)
            group.setAttribute("identificationQValue", score);
        for (String accession : parts) {
            Item protein = makeProtein(accession);
            group.addToCollection("proteins", protein);
        }
        store(group);
        proteinGroups.put(sortedAccessions, group);
        return group;
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
                    Item sample = metadata.samples.get(sampleName);
                    if (sample == null) {
                        throw new RuntimeException("Sample name not found in metadata: " + sampleName);
                    }
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
        int scoreIndex = proteinHeader.get("best_search_engine_score[1]");
        int membersIndex = proteinHeader.get("ambiguity_members");
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
                Item group = makeProteinGroup(line[membersIndex], ",", line[scoreIndex]);
                for (int i = 0; i < samples.size(); i++) {
                    String value = line[studyVariableIndexes[i]];
                    // treat "0" as missing value (see https://github.com/OpenMS/OpenMS/issues/6363):
                    if (!value.equals("null") && !value.equals("0.0")) {
                        Item abundance = createItem("ProteomicsLFQAbundance");
                        abundance.setAttribute("abundance", value);
                        abundance.setReference("sample", samples.get(i));
                        abundance.setReference("proteinGroup", group);
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
        // this header is followed by additional condition-specific columns:
        String[] expectedHeader = {"Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue",
                                   "adj.pvalue", "issue", "MissingPercentage", "ImputationPercentage"};
        String[] header = lineIter.next();
        if (header.length < expectedHeader.length) {
            throw new RuntimeException("Unexpected header in MSstats file");
        }
        for (int i = 0; i < expectedHeader.length; i++) {
            if (!header[i].equals(expectedHeader[i])) {
                throw new RuntimeException("Unexpected column name in MSstats file: '" +
                                           header[i] + "' vs. '" + expectedHeader[i] + "'");
            }
        }
        // trailing column names are condition names:
        for (int i = expectedHeader.length; i < header.length; i++) {
            // spaces in condition names are replaced when metadata is read from YAML:
            String condition = header[i].replace(' ', '_');
            if (!metadata.conditions.containsKey(condition)) {
                throw new RuntimeException("Condition name not found in metadata: " + condition);
            }
        }
        // "Label" column contains conditions (treatment/control) separated by "-";
        // since "-" may also occur in condition names, splitting doesn't work;
        // -> generate all possible pairs and perform look-up:
        Map<String, SimpleEntry<Item, Item>> conditionPairs = new HashMap<>();
        for (Map.Entry<String, Item> treatment : metadata.conditions.entrySet()) {
            for (Map.Entry<String, Item> control : metadata.conditions.entrySet()) {
                if (!treatment.getKey().equals(control.getKey())) {
                    SimpleEntry<Item, Item> pair = new SimpleEntry(treatment.getValue(), control.getValue());
                    conditionPairs.put(treatment.getKey() + "-" + control.getKey(), pair);
                }
            }
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }
            // @TODO: check if "issue" column is non-NA instead?
            if (line[6].equals("NA"))
                continue; // no p-value due to missing data
            // sort list of accessions to look up protein group:
            String accessions = line[0];
            String[] parts = accessions.replace(" ", "").split(";");
            Arrays.sort(parts);
            String groupKey = String.join(";", parts);
            Item group = proteinGroups.get(groupKey);
            if (group == null) {
                throw new RuntimeException("Protein group not found: " + groupKey);
            }
            Item result = createItem("ProteomicsMSstatsResult");
            result.setReference("experiment", experiment);
            result.setReference("proteinGroup", group);
            String conditionLabel = line[1].replace(' ', '_');
            SimpleEntry<Item, Item> conditionPair = conditionPairs.get(conditionLabel);
            if (conditionPair == null) {
                throw new RuntimeException("Failed to resolve conditions in 'Label' column: " + conditionLabel);
            }
            result.setReference("conditionTreatment", conditionPair.getKey());
            result.setReference("conditionControl", conditionPair.getValue());
            result.setAttribute("log2FoldChange", line[2]);
            result.setAttribute("rawPValue", line[6]);
            result.setAttribute("adjPValue", line[7]); // Benjamini-Hochberg (FDR) adjusted
            // fractions are relative to total number of quant. features (over all runs):
            result.setAttribute("missingFraction", line[9]);
            result.setAttribute("imputedFraction", line[10]);
            // @TODO: include any of these values?
            // result.setAttribute("standardError", line[3]); // SE of the log2FC
            // result.setAttribute("tTestValue", line[4]); // Student T-test statistic
            // result.setAttribute("degreesFreedom", line[5]); // DF of the T-test
            store(result);
        }
    }
}
