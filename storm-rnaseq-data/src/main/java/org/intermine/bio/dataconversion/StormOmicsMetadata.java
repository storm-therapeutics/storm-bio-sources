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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

// for troubleshooting (see below):
// import java.lang.ClassLoader;
// import java.lang.System;
// import java.net.URL;

import org.apache.log4j.Logger;

import org.json.*;

import org.intermine.dataconversion.DataConverter;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Read omics experiment metadata from JSON and store in InterMine
 *
 * The metadata file is validated against the schema and additional constraints before export to JSON, so we can assume it is "correct" and do not need to check.
 *
 * @author Adrian Bazaga, Hendrik Weisser
 */
public class StormOmicsMetadata
{
    public class ConditionsPair
    {
        public String control;
        public String treatment;

        public ConditionsPair(String treatment, String control) {
            this.treatment = treatment;
            this.control = control;
        }
    }

    // TODO: make variables private and add getters
    public String experimentShortName;
    public Item experiment;
    public String species;
    public Map<String, Item> materials = new HashMap<>();
    public Map<String, Item> treatments = new HashMap<>();
    public Map<String, Item> conditions = new HashMap<>();
    public Map<String, Item> samples = new HashMap<>();
    public Map<String, ArrayList<String>> bioReplicates = new HashMap<>(); // replicate name -> list of sample names
    public ArrayList<ConditionsPair> comparisons = new ArrayList<>();

    private DataConverter converter; // needed to create and store Items
    private static final Logger LOG = Logger.getLogger(StormOmicsMetadata.class);

    public StormOmicsMetadata(DataConverter converter) {
        this.converter = converter;

        // investigate "java.lang.NoSuchMethodError: 'java.util.Set org.json.JSONObject.keySet()'"
        // (caused by different versions of org.json package):
    //     ClassLoader classloader = org.json.JSONObject.class.getClassLoader();
    //     URL res = classloader.getResource("org/json/JSONObject.class");
    //     String path = res.getPath();
    //     System.out.println("Core JSONObject came from: " + path);
    }


    /// Convert a string (YAML entry) into an attribute name for an InterMine Item
    public static String convertToAttributeName(String name) {
        if (name.isEmpty())
            return name;
        StringBuilder builder = new StringBuilder();
        builder.append(name.substring(0, 1).toLowerCase());
        boolean afterSpace = false;
        for (int i = 1; i < name.length(); ++i) {
            char currentChar = name.charAt(i);
            if (currentChar == ' ') { // skip space, capitalize next letter
                afterSpace = true;
            }
            else if (afterSpace) {
                builder.append(Character.toUpperCase(currentChar));
                afterSpace = false;
            }
            else {
                builder.append(currentChar);
            }
        }
        return builder.toString();
    }


    private void extractAttributesFromJSON(JSONObject json, Iterable<String> expectedEntries, Item item) {
        extractAttributesFromJSON(json, expectedEntries, item, null);
    }


    private void extractAttributesFromJSON(JSONObject json, Iterable<String> expectedEntries,
                                           Item item, String prefix) {
        for (String entry : expectedEntries) {
            if (json.has(entry)) {
                String value = json.getString(entry);
                String name;
                if (prefix == null) {
                    name = convertToAttributeName(entry);
                }
                else {
                    name = convertToAttributeName(prefix + " " + entry);
                }
                item.setAttribute(name, value);
            }
        }
    }


    public void processJSONFile(File jsonFile) throws ObjectStoreException, IOException {
        LOG.info("Reading metadata for omics experiment from file: " + jsonFile.getName());
        try (Reader reader = new FileReader(jsonFile)) {
            JSONTokener tokener = new JSONTokener(reader);
            JSONObject json = new JSONObject(tokener);

            // Get the experiment metadata
            JSONObject experimentJson = json.getJSONObject("experiment");
            experimentShortName = experimentJson.getString("short name");
            LOG.debug("StormOmicsMetadata [processJSONFile] - processing experiment: " + experimentShortName);
            experiment = converter.createItem("StormOmicsExperiment");
            experiment.setAttribute("shortName", experimentShortName);
            species = experimentJson.optString("species", "human").toLowerCase();
            experiment.setAttribute("species", species);
            List<String> expectedEntries = List.of("name", "project", "contact person", "date",
                                                   "provider", "sequencing", "Dotmatics reference");
            extractAttributesFromJSON(experimentJson, expectedEntries, experiment);
            converter.store(experiment);

            // Process materials
            LOG.debug("StormOmicsMetadata [processJSONFile] - processing materials: " + experimentShortName);
            JSONObject materialsJson = json.getJSONObject("materials");
            processMetadataMaterials(materialsJson);

            // Process treatments
            LOG.debug("StormOmicsMetadata [processJSONFile] - processing treatments: " + experimentShortName);
            JSONObject treatmentsJson = json.getJSONObject("treatments");
            processMetadataTreatments(treatmentsJson);

            // Process conditions
            LOG.debug("StormOmicsMetadata [processJSONFile] - processing conditions: " + experimentShortName);
            JSONObject conditionsJson = json.getJSONObject("conditions");
            processMetadataConditions(conditionsJson);

            // Process each comparison individually
            // TODO: comparisons don't get stored in the database - do we need to process them at all?
            LOG.debug("StormOmicsMetadata [processJSONFile] - processing comparisons: " + experimentShortName);
            JSONArray experimentComparisons = json.getJSONArray("comparisons");
            for (int i = 0; i < experimentComparisons.length(); i++) {
                JSONObject comparison = experimentComparisons.getJSONObject(i);
                JSONObject treatment = comparison.getJSONObject("treatment");
                JSONObject control = comparison.getJSONObject("control");
                comparisons.add(new ConditionsPair(treatment.getString("name"),
                                                   control.getString("name")));
            }
        }
        catch (JSONException err) {
            throw new RuntimeException("Failed to read JSON from file. Error was: ", err);
        }
    }


    private void processMetadataMaterials(JSONObject materialsJson) {
        // Iterate over each material
        for (String materialName : materialsJson.keySet()) {
            try {
                JSONObject material = materialsJson.getJSONObject(materialName);
                // There should only be one key under this - what type of material:
                String materialType = material.keys().next();

                // Save the item
                Item materialItem = converter.createItem("StormOmicsMaterial");
                materialItem.setAttribute("name", materialName);
                materialItem.setAttribute("materialType", materialType);
                materialItem.setReference("experiment", experiment);

                JSONObject materialDetails = material.getJSONObject(materialType);
                extractAttributesFromJSON(materialDetails, List.of("tissue"), materialItem);
                if (materialType.equals("cell line")) {
                    extractAttributesFromJSON(materialDetails, List.of("name"), materialItem, "cell line");
                }
                else if (materialType.equals("tumour")) {
                    extractAttributesFromJSON(materialDetails, List.of("primary disease", "disease subtype"),
                                              materialItem, "tumour");
                }

                converter.store(materialItem);
                materials.put(materialName, materialItem);
            }
            catch (Exception e) {
                LOG.error("Exception in processMetadataMaterials with key: " +
                          materialName + " - " + e.getMessage());
                continue;
            }
        }
    }


    private void processMetadataTreatments(JSONObject treatmentsJson) {
        // Iterate over each treatment
        for (String treatmentName : treatmentsJson.keySet()) {
            try {
                JSONObject treatment = treatmentsJson.getJSONObject(treatmentName);
                // There should only be one key under this - what type of treatment:
                String treatmentType = treatment.keys().next();

                // Save the item
                Item treatmentItem = converter.createItem("StormOmicsTreatment");
                treatmentItem.setAttribute("name", treatmentName);
                treatmentItem.setAttribute("treatmentType", treatmentType);
                treatmentItem.setReference("experiment", experiment);

                ArrayList<String> expectedEntries = new ArrayList<>(List.of("name", "time point"));
                if (treatmentType.equals("inhibitor") || treatmentType.equals("activator")) {
                    expectedEntries.add("target gene");
                    expectedEntries.add("Dotmatics reference");
                }
                // nothing to add for types "knock-down", "overexpression", "untargeted" - but see below

                JSONObject treatmentDetails = treatment.getJSONObject(treatmentType);
                extractAttributesFromJSON(treatmentDetails, expectedEntries, treatmentItem);
                // some entries/attributes need special handling because the names don't match exactly:
                if (treatmentDetails.has("dose")) { // type "inhibitor" or "activator"
                    treatmentItem.setAttribute("dose_concentration", treatmentDetails.getString("dose"));
                }
                else if (treatmentDetails.has("concentration")) { // other types
                    treatmentItem.setAttribute("dose_concentration", treatmentDetails.getString("concentration"));
                }
                // type "knock-down" or "overexpression":
                extractAttributesFromJSON(treatmentDetails, List.of("type"), treatmentItem, "perturbation");

                converter.store(treatmentItem);
                treatments.put(treatmentName, treatmentItem);
            }
            catch (Exception e) {
                LOG.error("Exception in processMetadataTreatments with key: " +
                          treatmentName + " - " + e.getMessage());
                continue;
            }
        }
    }


    private void processMetadataConditions(JSONObject conditionsJson) {
        // Iterate over each condition
        for (String conditionName : conditionsJson.keySet()) {
            try {
                JSONObject condition = conditionsJson.getJSONObject(conditionName);

                // Save the item
                Item conditionItem = converter.createItem("StormOmicsCondition");
                conditionItem.setAttribute("name", conditionName);

                // exactly one material name per condition:
                String materialName = condition.getString("material");
                conditionItem.setReference("material", materials.get(materialName));

                // zero or more treatment names:
                JSONArray treatmentsJson = condition.optJSONArray("treatments");
                if (treatmentsJson != null) {
                    for (int i = 0; i < treatmentsJson.length(); i++) {
                        String treatmentName = treatmentsJson.getString(i);
                        conditionItem.addToCollection("treatments", treatments.get(treatmentName));
                    }
                }

                // Samples
                JSONObject samplesJson = condition.getJSONObject("samples");
                int repCounter = 1; // counter for biological replicates
                for (String sampleName : samplesJson.keySet()) {
                    JSONObject sample = samplesJson.getJSONObject(sampleName);
                    String label = sample.optString("label");
                    // names for biological replicates as used for Nextflow RNA-seq analysis pipeline
                    // (see "validate_yaml.py" script):
                    String bioReplicate = conditionName + "_R" + repCounter;
                    repCounter++;

                    ArrayList<String> replicates = new ArrayList<String>(); // technical replicates
                    if (sample.has("file")) { // legacy format
                        replicates.add(sample.getString("file"));
                    }
                    else {
                        JSONArray replicatesJson = sample.getJSONArray("replicates");
                        for (int i = 0; i < replicatesJson.length(); i++) {
                            replicates.add(replicatesJson.getString(i));
                        }
                    }
                    // flatten the "sample: replicates" structure (JSON) and just store replicates as samples:
                    for (int i = 0; i < replicates.size(); i++) {
                        String replicateFile = replicates.get(i);
                        String replicateName = replicateFile.split("\\.")[0]; // remove file extension(s)
                        Item sampleItem = converter.createItem("StormOmicsSample");
                        sampleItem.setReference("experiment", experiment);
                        sampleItem.setReference("condition", conditionItem);
                        sampleItem.setAttribute("name", replicateName);
                        sampleItem.setAttribute("file", replicateFile);
                        sampleItem.setAttribute("bioReplicate", sampleName);
                        if (!label.isEmpty()) {
                            sampleItem.setAttribute("label", label);
                        }
                        converter.store(sampleItem);
                        samples.put(replicateName, sampleItem);
                        replicates.set(i, replicateName); // for use in 'bioReplicates'
                    }
                    bioReplicates.put(bioReplicate, replicates);
                }

                converter.store(conditionItem);
                conditions.put(conditionName, conditionItem);
            }
            catch (Exception e) {
                LOG.error("Exception in processMetadataConditions with key: " +
                          conditionName + " - " + e.getMessage());
                continue;
            }
        }
    }
}
