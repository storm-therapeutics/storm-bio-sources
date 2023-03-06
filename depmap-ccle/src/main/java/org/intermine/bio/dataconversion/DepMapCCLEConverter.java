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
import java.io.Reader;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Read data from DepMap/CCLE
 * @author Hendrik Weisser
 */
public class DepMapCCLEConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap/CCLE data";
    private static final String DATA_SOURCE_NAME = "DepMap";
    private static final String HUMAN_TAXON_ID = "9606";

    private static final String CELL_LINES_FILE = "Model.csv";
    private static final String CRISPR_FILE = "CRISPRGeneEffect.csv";
    private static final String RNAI_FILE = "D2_combined_gene_dep_scores.csv";
    private static final String EXPRESSION_FILE = "OmicsExpressionProteinCodingGenesTPMLogp1.csv";
    private static final String COPY_NUMBER_FILE = "OmicsCNGene.csv";
    private static final String DAMAGING_MUT_FILE = "OmicsSomaticMutationsMatrixDamaging.csv";
    private static final String HOTSPOT_MUT_FILE = "OmicsSomaticMutationsMatrixHotspot.csv";
    private static final String CRISPR_PRED_FILE = "Chronos_Combined_predictability_results.csv";
    private static final String RNAI_PRED_FILE = "RNAi_merged_predictability_results.csv";


    private static final Logger LOG = Logger.getLogger(DepMapCCLEConverter.class);

    protected GeneLookup geneLookup;
    // mapping: cell line -> gene -> data ('DepMapCCLEData' item)
    protected Map<String, HashMap<String, Item>> cellLineData = new HashMap<String, HashMap<String, Item>>();
    protected Map<String, Item> cellLines = new HashMap<String, Item>();
    // mapping: cell line CCLE name (used in RNAi data) -> DepMap ID (used everywhere else)
    protected Map<String, String> cellLineCCLENames = new HashMap<String, String>();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepMapCCLEConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * Process directory containing DepMap/CCLE files
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        geneLookup = new GeneLookup(this);
        readCellLineInfo(new File(dataDir, CELL_LINES_FILE));
        readDataMatrix(new File(dataDir, EXPRESSION_FILE), "geneExpression", null);
        readDataMatrix(new File(dataDir, COPY_NUMBER_FILE), "copyNumber", null);
        readDataMatrix(new File(dataDir, CRISPR_FILE), "crisprGeneEffect", null);
        // TODO: convert mutation values (expected: 0.0, 1.0, 2.0) to integer
        readDataMatrix(new File(dataDir, DAMAGING_MUT_FILE), "hasDamagingMutation", "0.0");
        readDataMatrix(new File(dataDir, HOTSPOT_MUT_FILE), "hasHotspotMutation", "0.0");
        readRNAiData(new File(dataDir, RNAI_FILE));
        readPredictabilityData(new File(dataDir, CRISPR_PRED_FILE), "CRISPR");
        readPredictabilityData(new File(dataDir, RNAI_PRED_FILE), "RNAi");

        // store data:
        LOG.info("Storing DepMap/CCLE data...");
        for (Item cellLine : cellLines.values()) {
            store(cellLine);
        }
        for (Map.Entry<String, HashMap<String, Item>> entry1 : cellLineData.entrySet()) {
            Item cellLine = cellLines.get(entry1.getKey());
            for (Map.Entry<String, Item> entry2 : entry1.getValue().entrySet()) {
                // this creates and stores the gene Item (if it doesn't exist yet):
                Item gene = geneLookup.getGene(entry2.getKey());
                Item data = entry2.getValue();
                data.setReference("cellLine", cellLine);
                data.setReference("gene", gene);
                store(data);
            }
            entry1.getValue().clear(); // free up space
        }
        // do this last to include genes that were just looked up above:
        geneLookup.storeAllGeneItems();
    }


    private void readCellLineInfo(File inputPath) throws Exception {
        if (!inputPath.exists()) {
            throw new RuntimeException("DepMap/CCLE cell line info file not found: " + inputPath);
        }
        LOG.info("Processing DepMap/CCLE cell line info file: " + inputPath);
        FileReader reader = new FileReader(inputPath);
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // file header:
        // ModelID,PatientID,CellLineName,StrippedCellLineName,Age,SourceType,SangerModelID,RRID,DepmapModelType,GrowthPattern,MolecularSubtype,PrimaryOrMetastasis,SampleCollectionSite,Sex,SourceDetail,CatalogNumber,CCLEName,COSMICID,PublicComments,WTSIMasterCellID,OncotreeCode,OncotreeSubtype,OncotreePrimaryDisease,OncotreeLineage
        // we're only interested in some of the columns:
        Map<Integer, String> columnIndexes = new HashMap<Integer, String>(15);
        columnIndexes.put(0, "ModelID");
        columnIndexes.put(2, "CellLineName");
        columnIndexes.put(3, "StrippedCellLineName");
        columnIndexes.put(4, "Age");
        columnIndexes.put(7, "RRID");
        columnIndexes.put(8, "DepmapModelType");
        columnIndexes.put(9, "GrowthPattern");
        columnIndexes.put(10, "MolecularSubtype");
        columnIndexes.put(11, "PrimaryOrMetastasis");
        columnIndexes.put(12, "SampleCollectionSite");
        columnIndexes.put(13, "Sex");
        columnIndexes.put(16, "CCLEName");
        columnIndexes.put(21, "OncotreeSubtype");
        columnIndexes.put(22, "OncotreePrimaryDisease");
        columnIndexes.put(23, "OncotreeLineage");
        String[] header = lineIter.next();
        if (header.length < 24) {
            throw new RuntimeException("Unexpected end of header in DepMap/CCLE cell line info file");
        }
        for (Map.Entry<Integer, String> entry : columnIndexes.entrySet()) {
            int index = entry.getKey();
            if (!header[index].equals(entry.getValue())) {
                throw new RuntimeException("Unexpected entry at pos. " + index +
                                           " in header of DepMap/CCLE cell line info file: " + header[index]);
            }
        }
        // replace column names with attribute names:
        columnIndexes.replace(0, "depmapID");
        columnIndexes.replace(2, "name");
        columnIndexes.replace(3, "strippedName");
        columnIndexes.put(7, "rrid");
        columnIndexes.put(16, "ccleName");
        // lowercase first character:
        columnIndexes.replaceAll((key, value) ->
                                 value.substring(0, 1).toLowerCase() + value.substring(1));

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }
            Item cellLine = createItem("CellLine");
            for (Map.Entry<Integer, String> entry : columnIndexes.entrySet()) {
                int index = entry.getKey();
                String value = line[index];
                if (!value.isEmpty()) {
                    String name = entry.getValue();
                    cellLine.setAttribute(name, value);
                }
            }
            cellLines.put(line[0], cellLine);
            if (!line[16].isEmpty()) { // "CCLEName" entry
                cellLineCCLENames.put(line[16], line[0]);
            }
        }

        /*
        String[] expectedColumns = {"ModelID", "CellLineName", "StrippedCellLineName", "Age", "RRID",
                                    "DepmapModelType", "GrowthPattern", "MolecularSubtype", "PrimaryOrMetastasis",
                                    "SampleCollectionSite", "Sex", "CCLEName", "OncotreeSubtype",
                                    "OncotreePrimaryDisease", "OncotreeLineage"};
        Map<String, Integer> columnIndexes = new HashMap<String, Integer>(expectedColumns.length);
        for (String colName : expectedColumns) {
            columnIndexes.put(colName, null);
        }
        String[] header = lineIter.next();
        for (int i = 0; i < header.length; i++) {
            columnIndexes.replace(header[i], i);
        }
        // any missing columns?
        Iterator<Map.Entry<String, Integer>> entryIter = columnIndexes.entrySet().iterator();
        while (entryIter.hasNext()) {
            Map.Entry<String, Integer> entry = entryIter.next();
            if (entry.getValue() == null) {
                LOG.warn("Column name not found in DepMap/CCLE cell line info file: " + entry.getKey());
                entryIter.remove();
            }
        }
        */
    }


    /**
     * Process DepMap/CCLE file containing numeric data
     *
     * Cell lines are in rows, genes are in columns.
     */
    private void readDataMatrix(File inputPath, String outputAttribute, String exclude) throws Exception {
        if (!inputPath.exists()) {
            throw new RuntimeException("DepMap/CCLE input file not found: " + inputPath);
        }
        LOG.info("Processing DepMap/CCLE input file: " + inputPath);
        FileReader reader = new FileReader(inputPath);
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // parse header with gene symbols and IDs:
        String[] header = lineIter.next();
        // NCBI IDs for genes corresponding to columns in the data matrix; null if gene not found:
        // (first column is DepMap IDs of cell lines)
        String[] geneIDs = new String[header.length - 1];
        for (int i = 1; i < header.length; i++) {
            // expected format: "symbol (ID)", or just "ID"
            String[] parts = header[i].split(" ", 2); // "2" means "max. one split"
            String symbol = null;
            String id = null;
            if (parts.length == 2) {
                symbol = parts[0];
                id = parts[1].substring(1, parts[1].length() - 1); // remove brackets
                if (id.equals("nan"))
                    id = null;
            }
            else { // only one part
                if (parts[0].startsWith("ENSG"))
                    id = parts[0];
                else
                    symbol = parts[0];
            }
            if ((id != null) && id.matches("\\d+")) { // all numeric -> NCBI ID
                geneIDs[i - 1] = id;
                continue;
            }
            Item geneItem = geneLookup.getGene(null, id, symbol);
            if (geneItem != null) {
                geneIDs[i - 1] = geneItem.getAttribute("primaryIdentifier").getValue();
            }
            else {
                LOG.warn("Could not resolve gene: " + header[i] + " - skipping");
            }
        }
        // check for duplicate gene IDs:
        Set<String> geneIDSet = new HashSet<String>();
        for (int i = 0; i < geneIDs.length; i++) {
            if ((geneIDs[i] != null) && !geneIDSet.add(geneIDs[i])) {
                LOG.warn("Duplicate gene ID found: " + geneIDs[i] + " - removing");
                geneIDs[i] = null;
            }
        }

        // parse rows with cell line data:
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }
            String cellLine = line[0];
            if (!cellLines.containsKey(cellLine)) {
                LOG.warn("Cell line not found: " + cellLine + " - skipping");
                continue;
            }
            HashMap<String, Item> geneData = cellLineData.get(cellLine);
            if (geneData == null) {
                geneData = new HashMap<String, Item>();
                cellLineData.put(cellLine, geneData);
            }

            for (int i = 1; i < line.length; i++) {
                String geneID = geneIDs[i - 1];
                // skip missing genes and empty/"uninformative" values:
                if ((geneID != null) && !line[i].isEmpty() && ((exclude == null) || !line[i].equals(exclude))) {
                    Item dataItem = geneData.get(geneID);
                    if (dataItem == null) {
                        dataItem = createItem("DepMapCCLEData");
                        geneData.put(geneID, dataItem);
                    }
                    dataItem.setAttribute(outputAttribute, line[i]);
                    // cell line and gene references are set later (befor storing)
                }
            }
        }
    }


    /**
     * Process RNAi data file
     *
     * Cell lines are in columns, genes are in rows.
     */
    private void readRNAiData(File inputPath) throws Exception {
        if (!inputPath.exists()) {
            LOG.warn("DepMap/CCLE RNAi file not found: " + inputPath);
            return;
        }
        LOG.info("Processing DepMap/CCLE RNAi input file: " + inputPath);
        FileReader reader = new FileReader(inputPath);
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        String[] header = lineIter.next();
        // gene data maps for cell lines in columns:
        List<HashMap<String, Item>> rnaiGeneData = new ArrayList<HashMap<String, Item>>(header.length - 1);
        for (int i = 1; i < header.length; i++) {
            String ccleName = header[i];
            // special cases where names don't quite match the sample info file:
            // (note: "KP1NL_PANCREAS" is not the same as "KP1N_PANCREAS")
            if (ccleName.equals("AZ521_STOMACH"))
                ccleName = "AZ521_SMALL_INTESTINE";
            else if (ccleName.equals("GISTT1_GASTROINTESTINAL_TRACT"))
                ccleName = "GISTT1_STOMACH";
            else if (ccleName.equals("MB157_BREAST"))
                ccleName = "MDAMB157_BREAST";
            else if (ccleName.equals("SW527_BREAST"))
                ccleName = "SW527_LARGE_INTESTINE";

            String cellLine = cellLineCCLENames.get(ccleName); // null if not found
            if (cellLine != null) {
                rnaiGeneData.add(cellLineData.get(cellLine));
            }
            else {
                LOG.warn("Cell line not found: " + ccleName + " - skipping");
                rnaiGeneData.add(null);
            }
        }

        Set<String> geneIDSet = new HashSet<String>();
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }

            String geneID = line[0].split(" ")[1];
            geneID = geneID.substring(1, geneID.length() - 1); // strip brackets
            if (geneID.contains("&")) {
                LOG.warn("Fusion gene not supported: " + line[0] + " - skipping");
                continue;
            }
            if (!geneIDSet.add(geneID)) {
                LOG.warn("Duplicate gene ID in: " + line[0] + " - skipping");
                continue;
            }

            for (int i = 1; i < line.length; i++) {
                if (!line[i].isEmpty() && !line[i].equals("NA")) {
                    HashMap<String, Item> geneData = rnaiGeneData.get(i - 1);
                    if (geneData == null)
                        continue; // not interested in cell lines that have only RNAi data
                    Item dataItem = geneData.get(geneID);
                    if (dataItem != null) { // not interested in genes that have only RNAi data
                        dataItem.setAttribute("rnaiGeneEffect", line[i]);
                    }
                }
            }
        }
    }


    private void readPredictabilityData(File inputPath, String type) throws Exception {
        if (!inputPath.exists()) {
            LOG.warn("DepMap/CCLE " + type + " predictability file not found: " + inputPath);
            return;
        }
        LOG.info("Processing DepMap/CCLE " + type + " predictability input file: " + inputPath);
        String attribute = "depmap" + type.substring(0, 1).toUpperCase() +
            type.substring(1).toLowerCase() + "Predictability";

        FileReader reader = new FileReader(inputPath);
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        String[] header = lineIter.next();
        if ((header.length < 4) || !header[0].equals("gene") || !header[1].equals("model") ||
            !header[2].equals("pearson") || !header[3].equals("best")) {
                throw new RuntimeException("Unexpected header format in DepMap/CCLE predictability file");
        }
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }
            if (!line[3].equals("True")) // is this the best model for the gene?
                continue;
            // gene column format: "symbol (ID)"
            String geneID = line[0].split(" ")[1];
            geneID = geneID.substring(1, geneID.length() - 1); // remove brackets
            Item gene = geneLookup.getGene(geneID);
            if (gene != null)
                gene.setAttribute(attribute, line[2]);
        }
    }
}
