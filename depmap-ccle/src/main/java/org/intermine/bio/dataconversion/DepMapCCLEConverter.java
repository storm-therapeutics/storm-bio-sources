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
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
// @TODO: replace with javafx.util.Pair (requires additional dependency)?
import java.util.AbstractMap.SimpleEntry;

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
public class DepmapCCLEConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap/CCLE data";
    private static final String DATA_SOURCE_NAME = "DepMap";
    private static final String HUMAN_TAXON_ID = "9606";

    private static final String SAMPLE_INFO_FILE = "sample_info.csv";
    private static final String EXPRESSION_FILE = "CCLE_expression.csv";
    private static final String COPY_NUMBER_FILE = "CCLE_gene_cn.csv";
    private static final String CRISPR_FILE = "CRISPR_gene_effect.csv";

    private static final Logger LOG = Logger.getLogger(DepmapCCLEConverter.class);

    protected GeneLookup geneLookup;
    // mapping: cell line -> gene -> data ('DepMapCCLEData' item)
    protected Map<String, Map<String, Item>> cellLineData = new HashMap<String, HashMap<>>();
    protected Map<String, Item> cellLines = new HashMap<String, Item>();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapCcleMatrixConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        geneLookup = new GeneLookup();
        readSampleInfo(new File(dataDir, SAMPLE_INFO_FILE));
        readDataMatrix(new File(dataDir, EXPRESSION_FILE), "", "geneExpression");
        readDataMatrix(new File(dataDir, COPY_NUMBER_FILE), "", "copyNumber");
        readDataMatrix(new File(dataDir, CRISPR_FILE), "DepMap_ID", "crisprGeneEffect");
    }


    private void readSampleInfo(File inputPath) throws RuntimeException {
        if (!inputPath.exists()) {
            throw new RuntimeException("DepMap/CCLE sample info file not found: " + inputPath);
        }

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
                LOG.warning("Column name not found in DepMap/CCLE sample info file: " + entry.getKey());
            }
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }
            Item cellLine = createItem("cellLine");
            for (Map.Entry<String, Integer> entry : columnIndexes.entrySet()) {
                String name = entry.getKey();
                int index = entry.getValue();
                String value = line[index];
                name = StormOmicsMetadata.convertToAttributeName(name.replace('_', ' '));
                if (!value.isEmpty()) {
                    cellLine.setAttribute(name, value);
                }
            }
            cellLines.put(line[0], cellLine);

            for (Item cellLine : cellLines.values()) {
                if (cellLine.hasAttribute("parentDepmapId")) { // create proper reference to parent cell line
                    String parentId = cellLine.getAttribute("parentDepmapId").getValue();
                    Item parentLine = cellLines.get(parentId);
                    if (parentLine != null) {
                        cellLine.setReference("parentCellLine", parentLine);
                    }
                    else {
                        LOG.warning("Parent cell line not found: " + parentId);
                    }
                    cellLine.removeAttribute("parentDepmapId");
                }
                store(cellLine);
            }
        }
    }


    private void readDataMatrix(File inputPath, String headerStart, String outputAttribute) throws RuntimeException {
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
        String[] geneIDs() = new String[header.length - 1];
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
            if (id.startswith("ENSG")) {
                geneItem = geneLookup.getGene(null, id, symbol);
                geneIDs[i - 1] = geneItem.getAttribute("primaryIdentifier").getValue();
            }
            else {
                LOG.warning("Unexpected gene ID format: " + id + " - skipping");
            }
        }
        // check for duplicate gene IDs:
        Set<String> geneIDSet = new HashSet<String>();
        for (int i = 0; i < geneIDs.length; i++) {
            if (!geneIDSet.add(geneIDs[i])) {
                LOG.warning("Duplicate gene ID found: " + geneIDs[i] + " - removing");
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
                LOG.warning("Cell line not found: " + cellLine + " - skipping");
                continue;
            }
            HashMap<String, Item> geneData = cellLineData.get(cellLine);
            if (geneData == null) {
                geneData = new HashMap<String, Item>();
                cellLineData.put(cellLine, geneData);
            }

            for (int i = 1; i < line.length; i++) {
                String geneID = geneIDs[i - 1];
                if ((geneID != null) && (!line[i].isEmpty())) {
                    Item dataItem = geneData.get(geneID);
                    if (dataItem == null) {
                        dataItem = createItem("DepMapCCLEData");
                        geneData.put(geneID, dataItem);
                    }
                    dataItem.setAttribute(outputAttribute, line[i]);
                    // dataItem.setReference("gene", geneItems[i - 1]);
                    // dataItem.setReference("cellLine", cellLine);
                    // store(dataItem);
                }
            }
        }
    }
}
