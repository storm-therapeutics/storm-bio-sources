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

import org.json.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;

import org.intermine.dataconversion.DataConverter;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Read experiment metadata from JSON and store in InterMine
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
    public Map<String, Item> materials = new HashMap<>();
    public Map<String, Item> treatments = new HashMap<>();
    public Map<String, Item> conditions = new HashMap<>();
    public ArrayList<ConditionsPair> comparisons = new ArrayList<>();

    private DataConverter converter; // needed to create and store Items
    private static final Logger LOG = Logger.getLogger(StormOmicsMetadata.class);

    public StormOmicsMetadata(DataConverter converter) {
        this.converter = converter;
    }


    public void processJSONFile(File jsonFile) throws ObjectStoreException, IOException {
        Reader reader = new FileReader(jsonFile);
        // TODO: read file content directly into JSON object (without BufferedReader/single-line String)
        try (BufferedReader br = new BufferedReader(reader)) {
            String line;
            while ((line = br.readLine()) != null) {
                try {
                    JSONObject jsonObject = new JSONObject(line);

                    // Get the experiment key
                    JSONObject experimentJson = jsonObject.getJSONObject("experiment");

                    // Get the experiment metadata
                    experimentShortName = experimentJson.getString("short name");
                    String experimentName = experimentJson.getString("name");
                    String experimentProject = experimentJson.getString("project");
                    String experimentContactPerson = experimentJson.getString("contact person");
                    String experimentDate = experimentJson.getString("date");
                    String experimentSequencing = experimentJson.getString("sequencing");
                    String experimentProvider = experimentJson.getString("provider");
                    String experimentDotmaticsReference = experimentJson.getString("Dotmatics reference");

                    LOG.info("StormRnaseqDataConverter [processJSONFile] - Processing 1: " + experimentShortName);

                    // Save the item
                    experiment = converter.createItem("RNASeqExperimentMetadata");

                    if (!experimentName.isEmpty()) {
                        experiment.setAttribute("name", experimentName);
                    }
                    if (!experimentShortName.isEmpty()) {
                        experiment.setAttribute("shortName", experimentShortName);
                    }
                    if (!experimentProject.isEmpty()) {
                        experiment.setAttribute("project", experimentProject);
                    }
                    if (!experimentContactPerson.isEmpty()) {
                        experiment.setAttribute("contactPerson", experimentContactPerson);
                    }
                    if (!experimentDate.isEmpty()) {
                        experiment.setAttribute("date", experimentDate);
                    }
                    if (!experimentSequencing.isEmpty()) {
                        experiment.setAttribute("sequencing", experimentSequencing);
                    }
                    if (!experimentProvider.isEmpty()) {
                        experiment.setAttribute("provider", experimentProvider);
                    }
                    if (!experimentDotmaticsReference.isEmpty()) {
                        experiment.setAttribute("dotmaticsReference", experimentDotmaticsReference);
                    }
                    //store(experiment);

                    // Process materials
                    LOG.info("StormRnaseqDataConverter [processJSONFile] - Processing 2: " + experimentShortName);
                    JSONObject materialsJson = jsonObject.getJSONObject("materials");
                    processRNASeqExperimentMaterials(materialsJson);

                    // Process treatments
                    LOG.info("StormRnaseqDataConverter [processJSONFile] - Processing 3: " + experimentShortName);
                    JSONObject treatmentsJson = jsonObject.getJSONObject("treatments");
                    processRNASeqExperimentTreatments(treatmentsJson);

                    // Process conditions
                    JSONObject conditionsJson = jsonObject.getJSONObject("conditions");
                    processRNASeqExperimentConditions(conditionsJson);

                    // Process each comparison individually
                    LOG.info("StormRnaseqDataConverter [processJSONFile] - Processing 4: " + experimentShortName);
                    JSONArray experimentComparisons = (JSONArray)jsonObject.get("comparisons");
                    for (int i = 0; i < experimentComparisons.length(); i++) {
                        JSONObject comparison = experimentComparisons.getJSONObject(i);
                        JSONObject treatmentObject = comparison.getJSONObject("treatment");
                        JSONObject controlObject = comparison.getJSONObject("control");
                        comparisons.add(new ConditionsPair(treatmentObject.getString("name"),
                                                           controlObject.getString("name")));
                    }
                }
                catch (JSONException err) {
                    throw new RuntimeException("Failed to read the following JSON: " + line, err);
                }
            }
        }
    }


    private void processRNASeqExperimentMaterials(JSONObject materialsJson) {
        // Iterate over each material
        Iterator<String> keys = materialsJson.keys();

        while (keys.hasNext()) {
            String key = keys.next();
            try {
                // Now we have the material
                String materialName = key;

                // There should only be one key under this
                JSONObject materialTypeKeysJSON = materialsJson.getJSONObject(key);

                if (materialTypeKeysJSON.has("cell line")) {
                    String materialType = "cell line";
                    // Get the material object
                    JSONObject materialObject = materialTypeKeysJSON.getJSONObject(materialType);

                    String cellLineName = "";
                    if (materialObject.has("name")) {
                        cellLineName = materialObject.getString("name");
                    }
                    String cellLineTissue = "";
                    if (materialObject.has("tissue")) {
                        cellLineTissue = materialObject.getString("tissue");
                    }
                    String cellLineSpecies = "";
                    if (materialObject.has("species")) {
                        cellLineSpecies = materialObject.getString("species");
                    }

                    // Save the item
                    Item materialMetadataItem = converter.createItem("RNASeqExperimentMaterial");

                    if (!materialType.isEmpty()) {
                        materialMetadataItem.setAttribute("materialType", materialType);
                    }
                    if (!cellLineName.isEmpty()) {
                        materialMetadataItem.setAttribute("name", cellLineName);
                    }
                    if (!cellLineTissue.isEmpty()) {
                        materialMetadataItem.setAttribute("tissue", cellLineTissue);
                    }
                    if (!cellLineSpecies.isEmpty()) {
                        materialMetadataItem.setAttribute("species", cellLineSpecies);
                    }

                    materialMetadataItem.setReference("experiment", experiment);
                    converter.store(materialMetadataItem);

                    if (!materials.containsKey(materialName)) {
                        materials.put(materialName, materialMetadataItem);
                    }
                } else if (materialTypeKeysJSON.has("tumour")) {
                    String materialType = "tumour";
                    // Get the material object
                    JSONObject materialObject = materialTypeKeysJSON.getJSONObject(materialType);

                    String tumourPrimaryDisease = "";
                    if (materialObject.has("primary disease")) {
                        tumourPrimaryDisease = materialObject.getString("primary disease");
                    }
                    String tumourDiseaseSubtype = "";
                    if (materialObject.has("disease subtype")) {
                        tumourDiseaseSubtype = materialObject.getString("disease subtype");
                    }
                    String tumourTissue = "";
                    if (materialObject.has("tissue")) {
                        tumourTissue = materialObject.getString("tissue");
                    }
                    String tumourSpecies = "";
                    if (materialObject.has("species")) {
                        tumourSpecies = materialObject.getString("species");
                    }

                    // Save the item
                    Item materialMetadataItem = converter.createItem("RNASeqExperimentMaterial");

                    if (!materialType.isEmpty()) {
                        materialMetadataItem.setAttribute("materialType", materialType);
                    }
                    if (!tumourPrimaryDisease.isEmpty()) {
                        materialMetadataItem.setAttribute("primaryDisease", tumourPrimaryDisease);
                    }
                    if (!tumourDiseaseSubtype.isEmpty()) {
                        materialMetadataItem.setAttribute("diseaseSubtype", tumourDiseaseSubtype);
                    }
                    if (!tumourTissue.isEmpty()) {
                        materialMetadataItem.setAttribute("tissue", tumourTissue);
                    }
                    if (!tumourSpecies.isEmpty()) {
                        materialMetadataItem.setAttribute("species", tumourSpecies);
                    }

                    materialMetadataItem.setReference("experiment", experiment);
                    converter.store(materialMetadataItem);

                    if (!materials.containsKey(materialName)) {
                        materials.put(materialName, materialMetadataItem);
                    }
                } else if (materialTypeKeysJSON.has("tissue")) {
                    String materialType = "tissue";
                    // Get the material object
                    JSONObject materialObject = materialTypeKeysJSON.getJSONObject(materialType);

                    String tissueTissue = "";
                    if (materialObject.has("tissue")) {
                        tissueTissue = materialObject.getString("tissue");
                    }
                    String tissueSpecies = "";
                    if (materialObject.has("species")) {
                        tissueSpecies = materialObject.getString("species");
                    }

                    // Save the item
                    Item materialMetadataItem = converter.createItem("RNASeqExperimentMaterial");

                    if (!materialType.isEmpty()) {
                        materialMetadataItem.setAttribute("materialType", materialType);
                    }
                    if (!tissueTissue.isEmpty()) {
                        materialMetadataItem.setAttribute("tissue", tissueTissue);
                    }
                    if (!tissueSpecies.isEmpty()) {
                        materialMetadataItem.setAttribute("species", tissueSpecies);
                    }

                    materialMetadataItem.setReference("experiment", experiment);
                    converter.store(materialMetadataItem);

                    if (!materials.containsKey(materialName)) {
                        materials.put(materialName, materialMetadataItem);
                    }
                }

            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentMaterials with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }


    private void processRNASeqExperimentTreatments(JSONObject treatmentsJson) {
        // Iterate over each material
        Iterator<String> keys = treatmentsJson.keys();

        while(keys.hasNext()) {
            String key = keys.next();
            try {
                // Now we have the material
                String treatmentName = key;

                // There should only be one key under this
                JSONObject treatmentTypeKeysJSON = treatmentsJson.getJSONObject(key);

                if(treatmentTypeKeysJSON.has("inhibitor")) {
                    String treatmentType = "inhibitor";
                    // Get the treatment object
                    JSONObject treatmentObject = treatmentTypeKeysJSON.getJSONObject(treatmentType);

                    String inhibitorName = "";
                    if (treatmentObject.has("name")) {
                        inhibitorName = treatmentObject.getString("name");
                    }
                    String targetGene = "";
                    if (treatmentObject.has("target gene")) {
                        targetGene = treatmentObject.getString("target gene");
                    }
                    String dotmaticsReference = "";
                    if (treatmentObject.has("Dotmatics reference")) {
                        dotmaticsReference = treatmentObject.getString("Dotmatics reference");
                    }
                    String dose = "";
                    if (treatmentObject.has("dose")) {
                        dose = treatmentObject.getString("dose");
                    }
                    String timePoint = "";
                    if (treatmentObject.has("time point")) {
                        timePoint = treatmentObject.getString("time point");
                    }

                    // Save the item
                    Item treatmentMetadataItem = converter.createItem("RNASeqExperimentTreatment");
                    treatmentMetadataItem.setAttribute("treatmentType", treatmentType);

                    if (!key.isEmpty()) {
                        treatmentMetadataItem.setAttribute("name", key);
                    }
                    if (!targetGene.isEmpty()) {
                        treatmentMetadataItem.setAttribute("targetGene", targetGene);
                    }
                    if (!dotmaticsReference.isEmpty()) {
                        treatmentMetadataItem.setAttribute("dotmaticsReference", dotmaticsReference);
                    }
                    if (!StringUtils.isEmpty(dose)) {
                        treatmentMetadataItem.setAttribute("dose_concentration", dose);
                    }
                    if (!timePoint.isEmpty()) {
                        treatmentMetadataItem.setAttribute("timePoint", timePoint);
                    }

                    treatmentMetadataItem.setReference("experiment", experiment);
                    converter.store(treatmentMetadataItem);

                    if (!treatments.containsKey(treatmentName)) {
                        treatments.put(treatmentName, treatmentMetadataItem);
                    }

                } else if(treatmentTypeKeysJSON.has("knock-down")) {
                    String treatmentType = "knock-down";
                    // Get the treatment object
                    JSONObject treatmentObject = treatmentTypeKeysJSON.getJSONObject(treatmentType);

                    String inhibitorName = "";
                    if (treatmentObject.has("name")) {
                        inhibitorName = treatmentObject.getString("name");
                    }
                    String targetGene = "";
                    if (treatmentObject.has("target gene")) {
                        targetGene = treatmentObject.getString("target gene");
                    }
                    String concentration = "";
                    if (treatmentObject.has("concentration")) {
                        concentration = treatmentObject.getString("concentration");
                    }
                    String type = "";
                    if (treatmentObject.has("type")) {
                        type = treatmentObject.getString("type");
                    }
                    String timePoint = "";
                    if (treatmentObject.has("time point")) {
                        timePoint = treatmentObject.getString("time point");
                    }

                    // Save the item
                    Item treatmentMetadataItem = converter.createItem("RNASeqExperimentTreatment");
                    treatmentMetadataItem.setAttribute("treatmentType", treatmentType);

                    if (!key.isEmpty()) {
                        treatmentMetadataItem.setAttribute("name", key);
                    }
                    if (!targetGene.isEmpty()) {
                        treatmentMetadataItem.setAttribute("targetGene", targetGene);
                    }
                    if (!StringUtils.isEmpty(concentration)) {
                        treatmentMetadataItem.setAttribute("dose_concentration", concentration);
                    }
                    if (!type.isEmpty()) {
                        treatmentMetadataItem.setAttribute("type", type);
                    }
                    if (!timePoint.isEmpty()) {
                        treatmentMetadataItem.setAttribute("timePoint", timePoint);
                    }

                    treatmentMetadataItem.setReference("experiment", experiment);
                    converter.store(treatmentMetadataItem);

                    if (!treatments.containsKey(treatmentName)) {
                        treatments.put(treatmentName, treatmentMetadataItem);
                    }
                }
                else if (treatmentTypeKeysJSON.has("untargeted")) {
                    String treatmentType = "untargeted";
                    // Get the treatment object
                    JSONObject treatmentObject = treatmentTypeKeysJSON.getJSONObject(treatmentType);

                    String inhibitorName = "";
                    if (treatmentObject.has("name")) {
                        inhibitorName = treatmentObject.getString("name");
                    }
                    String targetGene = "";
                    if (treatmentObject.has("target gene")) {
                        targetGene = treatmentObject.getString("target gene");
                    }
                    String concentration = "";
                    if (treatmentObject.has("concentration")) {
                        concentration = treatmentObject.getString("concentration");
                    }
                    String type = "";
                    if (treatmentObject.has("type")) {
                        type = treatmentObject.getString("type");
                    }
                    String timePoint = "";
                    if (treatmentObject.has("time point")) {
                        timePoint = treatmentObject.getString("time point");
                    }

                    // Save the item
                    Item treatmentMetadataItem = converter.createItem("RNASeqExperimentTreatment");
                    treatmentMetadataItem.setAttribute("treatmentType", treatmentType);

                    if (!key.isEmpty()) {
                        treatmentMetadataItem.setAttribute("name", key);
                    }
                    if (!StringUtils.isEmpty(concentration)) {
                        treatmentMetadataItem.setAttribute("dose_concentration", concentration);
                    }
                    if (!timePoint.isEmpty()) {
                        treatmentMetadataItem.setAttribute("timePoint", timePoint);
                    }

                    treatmentMetadataItem.setReference("experiment", experiment);
                    converter.store(treatmentMetadataItem);

                    if (!treatments.containsKey(treatmentName)) {
                        treatments.put(treatmentName, treatmentMetadataItem);
                    }
                }
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentTreatments with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }


    private void processRNASeqExperimentConditions(JSONObject conditionsJson) {
        // Iterate over each condition
        Iterator<String> keys = conditionsJson.keys();

        while(keys.hasNext()) {
            String key = keys.next();
            try {
                // Now we have the condition
                String conditionName = key;

                // There should only be one key under this
                JSONObject conditionsKeysJSON = conditionsJson.getJSONObject(conditionName);

                String materialName = conditionsKeysJSON.getString("material");

                // Samples
                ArrayList<String> samplesArray = new ArrayList<String>();
                JSONObject samplesJson = conditionsKeysJSON.getJSONObject("samples");
                Iterator<String> samplesKeys = samplesJson.keys();
                while(samplesKeys.hasNext()) {
                    String sampleKey = samplesKeys.next();
                    samplesArray.add(sampleKey);
                }

                String samples = String.join(", ", samplesArray);

                // Treatments
                ArrayList<String> treatmentsArray = new ArrayList<String>();
                JSONArray treatmentsJson = conditionsKeysJSON.getJSONArray("treatments");
                for(int i = 0; i < treatmentsJson.length(); i++) {
                    treatmentsArray.add(treatmentsJson.getString(i));
                }

                //String treatments = String.join(", ", treatmentsArray);

                // Save the item
                Item conditionMetadataItem = converter.createItem("RNASeqExperimentCondition");

                if (!treatmentsArray.isEmpty()) {
                    conditionMetadataItem.setAttribute("name", conditionName);
                } else {
                    continue;
                }

                if (!treatments.isEmpty()) {
                    //conditionMetadataItem.setAttribute("treatments", treatments);
                    List<String> treatmentsIds = new ArrayList<>();

                    for(int i = 0; i < treatmentsArray.size(); i++) {
                        String treatmentNameToAdd = treatmentsArray.get(i);
                        if (treatments.containsKey(treatmentNameToAdd)) {
                            String treatmentIdToAdd = treatments.get(treatmentNameToAdd).getIdentifier();
                            treatmentsIds.add(treatmentIdToAdd);
                        }
                    }

                    conditionMetadataItem.setCollection("treatments", treatmentsIds);
                }

                if (!samples.isEmpty()) {
                    conditionMetadataItem.setAttribute("samples", samples);
                }

                if (!materialName.isEmpty()) {
                    if (materials.containsKey(materialName)) {
                        conditionMetadataItem.setReference("material", materials.get(materialName));
                    }
                }

                conditionMetadataItem.setReference("experiment", experiment);

                converter.store(conditionMetadataItem);

                if (!conditions.containsKey(conditionName)) {
                    conditions.put(conditionName, conditionMetadataItem);
                }

            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentConditions with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }
}
