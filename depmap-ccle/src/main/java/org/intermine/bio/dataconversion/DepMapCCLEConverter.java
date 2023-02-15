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

    private static final String SAMPLE_INFO_FILE = "sample_info.csv";
    private static final String EXPRESSION_FILE = "CCLE_expression.csv";
    private static final String COPY_NUMBER_FILE = "CCLE_gene_cn.csv";
    private static final String CRISPR_FILE = "CRISPR_gene_effect.csv";
    private static final String RNAI_FILE = "D2_combined_gene_dep_scores.csv";
    private static final String MUTATIONS_FILE = "CCLE_mutations.csv";

    private static final Logger LOG = Logger.getLogger(DepMapCCLEConverter.class);

    private class MutationCounts {
        public int deleterious;
        public int tcga;
        public int cosmic;
        public int damaging;
        public int nonconserving;
        public int conserving;
        public int silent;
    };

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
        readSampleInfo(new File(dataDir, SAMPLE_INFO_FILE));
        readDataMatrix(new File(dataDir, EXPRESSION_FILE), "", "geneExpression");
        readDataMatrix(new File(dataDir, COPY_NUMBER_FILE), "", "copyNumber");
        readDataMatrix(new File(dataDir, CRISPR_FILE), "DepMap_ID", "crisprGeneEffect");
        readRNAiData(new File(dataDir, RNAI_FILE));
        readMutationData(new File(dataDir, MUTATIONS_FILE));

        // store data (cell lines are already done):
        Set<String> geneIDSet = new HashSet<String>();
        for (Map.Entry<String, HashMap<String, Item>> entry1 : cellLineData.entrySet()) {
            Item cellLine = cellLines.get(entry1.getKey());
            for (Map.Entry<String, Item> entry2 : entry1.getValue().entrySet()) {
                String geneID = entry2.getKey();
                Item gene = geneLookup.getGene(geneID);
                // store each gene only once:
                if (!geneIDSet.add(geneID)) {
                    store(gene);
                }
                Item data = entry2.getValue();
                data.setReference("cellLine", cellLine);
                data.setReference("gene", gene);
                store(data);
            }
            entry1.getValue().clear(); // free up space
        }
    }


    private void readSampleInfo(File inputPath) throws Exception {
        if (!inputPath.exists()) {
            throw new RuntimeException("DepMap/CCLE sample info file not found: " + inputPath);
        }
        LOG.info("Processing DepMap/CCLE sample information file: " + inputPath);
        FileReader reader = new FileReader(inputPath);
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // we're only interested in some of the columns:
        String[] expectedColumns = {"DepMap_ID", "cell_line_name", "stripped_cell_line_name", "CCLE_Name", "alias",
                                    "sex", "RRID", "sample_collection_site", "primary_or_metastasis",
                                    "primary_disease", "Subtype", "age", "lineage", "lineage_subtype",
                                    "lineage_sub_subtype", "lineage_molecular_subtype", "default_growth_pattern",
                                    "parent_depmap_id", "Cellosaurus_NCIt_disease"};
        Map<String, Integer> columnIndexes = new HashMap<String, Integer>(expectedColumns.length);
        for (String colName : expectedColumns) {
            columnIndexes.put(colName, null);
        }
        String[] header = lineIter.next();
        for (int i = 0; i < header.length; i++) {
            columnIndexes.replace(header[i], i);
        }
        for (Map.Entry<String, Integer> entry : columnIndexes.entrySet()) {
            if (entry.getValue() == null) {
                LOG.warn("Column name not found in DepMap/CCLE sample info file: " + entry.getKey());
            }
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }
            Item cellLine = createItem("cellLine");
            for (Map.Entry<String, Integer> entry : columnIndexes.entrySet()) {
                int index = entry.getValue();
                String value = line[index];
                if (!value.isEmpty()) {
                    if (value.startsWith("\"")) { // strip quotes, if necessary
                        value = value.substring(1, value.length() - 1);
                    }
                    String name = entry.getKey();
                    // use better attribute names:
                    if (name.equals("cell_line_name"))
                        name = "name";
                    else if (name.equals("stripped_cell_line_name"))
                        name = "strippedName";
                    else if (name.equals("CCLE_Name"))
                        name = "ccleName";
                    else if (name.equals("RRID"))
                        name = "rrid";
                    else
                        name = StormOmicsMetadata.convertToAttributeName(name.replace('_', ' '));
                    cellLine.setAttribute(name, value);
                }
            }
            cellLines.put(line[0], cellLine);
            int index = columnIndexes.get("CCLE_Name");
            if (!line[index].isEmpty()) {
                cellLineCCLENames.put(line[index], line[0]);
            }
        }

        for (Item cellLine : cellLines.values()) {
            if (cellLine.hasAttribute("parentDepmapId")) { // create proper reference to parent cell line
                String parentId = cellLine.getAttribute("parentDepmapId").getValue();
                Item parentLine = cellLines.get(parentId);
                if (parentLine != null) {
                    cellLine.setReference("parentCellLine", parentLine);
                }
                else {
                    LOG.warn("Parent cell line not found: " + parentId);
                }
                cellLine.removeAttribute("parentDepmapId");
            }
            store(cellLine);
        }
    }


    /**
     * Process DepMap/CCLE file containing numeric data
     *
     * Cell lines are in rows, genes are in columns.
     */
    private void readDataMatrix(File inputPath, String headerStart, String outputAttribute) throws Exception {
        if (!inputPath.exists()) {
            throw new RuntimeException("DepMap/CCLE input file not found: " + inputPath);
        }
        LOG.info("Processing DepMap/CCLE input file: " + inputPath);
        FileReader reader = new FileReader(inputPath);
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // parse header with gene symbols and IDs:
        String[] header = lineIter.next();
        if ((headerStart != null) && !header[0].equals(headerStart)) {
            throw new RuntimeException("Unexpected item at start of input file");
        }
        // NCBI IDs for genes corresponding to columns in the data matrix; null if gene not found:
        String[] geneIDs = new String[header.length - 1];
        for (int i = 1; i < header.length; i++) {
            // expected format: "symbol (ID)", or just "ID"
            String[] parts = header[i].split(" ", 2); // "2" means "max. one split"
            String symbol = null;
            String id = null;
            if (parts.length == 2) {
                symbol = parts[0];
                id = parts[1].substring(1, parts[1].length() - 1); // remove brackets
            }
            else {
                id = parts[0];
            }
            if (id.matches("\\d+")) { // all numeric -> NCBI ID
                geneIDs[i - 1] = id;
                continue;
            }
            Item geneItem;
            if (id.startsWith("ENSG")) {
                geneItem = geneLookup.getGene(null, id, symbol);
                geneIDs[i - 1] = geneItem.getAttribute("primaryIdentifier").getValue();
            }
            else {
                LOG.warn("Unexpected gene ID format: " + id + " - skipping");
            }
        }
        // check for duplicate gene IDs:
        Set<String> geneIDSet = new HashSet<String>();
        for (int i = 0; i < geneIDs.length; i++) {
            if (!geneIDSet.add(geneIDs[i])) {
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
                if ((geneID != null) && !line[i].isEmpty()) {
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
            String ccleName = header[i].substring(1, header[i].length() - 1); // strip quotes
            String cellLine = cellLineCCLENames.get(ccleName); // null if not found
            if (cellLine != null) {
                rnaiGeneData.add(cellLineData.get(cellLine));
            }
            else {
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
            geneID = geneID.substring(1, geneID.length() - 2); // strip brackets and closing quote
            if (!geneIDSet.add(geneID)) {
                LOG.warn("Duplicate gene ID found: " + geneID + " - skipping");
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


    private void readMutationData(File inputPath) throws Exception {
        if (!inputPath.exists()) {
            LOG.warn("DepMap/CCLE mutations file not found: " + inputPath);
            return;
        }
        LOG.info("Processing DepMap/CCLE mutations input file: " + inputPath);
        FileReader reader = new FileReader(inputPath);
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        List<String> headerList = Arrays.asList(lineIter.next());
        int geneIndex = headerList.indexOf("Entrez_Gene_Id");
        int cellLineIndex = headerList.indexOf("DepMap_ID");
        int deleteriousIndex = headerList.indexOf("isDeleterious");
        int tcgaIndex = headerList.indexOf("isTCGAhotspot");
        int cosmicIndex = headerList.indexOf("isCOSMIChotspot");
        int variantIndex = headerList.indexOf("Variant_annotation");

        if ((geneIndex < 0) || (cellLineIndex < 0) || (deleteriousIndex < 0) || (tcgaIndex < 0) ||
            (cosmicIndex < 0) || (variantIndex < 0)) {
            throw new RuntimeException("Unexpected header in DepMap/CCLE mutations file");
        }

        // mapping: cell line -> gene -> mutation counts
        Map<String, HashMap<String, MutationCounts>> cellLineCounts = new HashMap<String, HashMap<String, MutationCounts>>();
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != headerList.size()) {
                throw new RuntimeException("Unexpected number of items per line");
            }
            String geneID = line[geneIndex];
            if (geneID.equals("0"))
                continue;
            String cellLine = line[cellLineIndex];
            HashMap<String, MutationCounts> geneCounts = cellLineCounts.get(cellLine);
            if (geneCounts == null) {
                geneCounts = new HashMap<String, MutationCounts>();
                cellLineCounts.put(cellLine, geneCounts);
            }
            MutationCounts counts = geneCounts.get(geneID);
            if (counts == null) {
                counts = new MutationCounts();
                geneCounts.put(geneID, counts);
            }
            if (line[deleteriousIndex].equals("True"))
                counts.deleterious++;
            if (line[tcgaIndex].equals("True"))
                counts.tcga++;
            if (line[cosmicIndex].equals("True"))
                counts.cosmic++;
            String variant = line[variantIndex];
            if (variant.equals("damaging"))
                counts.damaging++;
            else if (variant.equals("other non-conserving"))
                counts.nonconserving++;
            else if (variant.equals("other conserving"))
                counts.conserving++;
            else if (variant.equals("silent"))
                counts.silent++;
        }

        // combine mutation data with other data:
        for (Map.Entry<String, HashMap<String, MutationCounts>> entry1 : cellLineCounts.entrySet()) {
            HashMap<String, Item> geneData = cellLineData.get(entry1.getKey());
            if (geneData == null) // no other data for this cell line - skip
                continue;
            for (Map.Entry<String, MutationCounts> entry2 : entry1.getValue().entrySet()) {
                Item data = geneData.get(entry2.getKey());
                if (data == null) // no other data for this gene - skip (TODO: or not?)
                    continue;
                MutationCounts counts = entry2.getValue();
                if (counts.deleterious > 0)
                    data.setAttribute("mutationCountDeleterious", String.valueOf(counts.deleterious));
                if (counts.tcga > 0)
                    data.setAttribute("mutationCountHotspotTCGA", String.valueOf(counts.tcga));
                if (counts.cosmic > 0)
                    data.setAttribute("mutationCountHotspotCOSMIC", String.valueOf(counts.cosmic));
                if (counts.damaging > 0)
                    data.setAttribute("mutationCountVariantDamaging", String.valueOf(counts.damaging));
                if (counts.nonconserving > 0)
                    data.setAttribute("mutationCountVariantOtherNonconserving",
                                      String.valueOf(counts.nonconserving));
                if (counts.conserving > 0)
                    data.setAttribute("mutationCountVariantOtherConserving", String.valueOf(counts.conserving));
                if (counts.silent > 0)
                    data.setAttribute("mutationCountVariantSilent", String.valueOf(counts.silent));
            }
        }
    }
}
