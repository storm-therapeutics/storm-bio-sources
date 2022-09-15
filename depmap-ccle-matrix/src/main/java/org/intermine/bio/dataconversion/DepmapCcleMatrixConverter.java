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
 * Read data from DepMap/CCLE in simple matrix format
 * @author Hendrik Weisser
 */
public class DepmapCcleMatrixConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap/CCLE data";
    private static final String DATA_SOURCE_NAME = "DepMap";
    private static final String HUMAN_TAXON_ID = "9606";

    private static final Logger LOG = Logger.getLogger(DepmapCcleMatrixConverter.class);

    protected IdResolver resolver;
    protected String outputAttribute;
    protected String inputFile;
    protected String headerStart = "";
    protected boolean idsArePrimary = false;

    // genes corresponding to columns in the data matrix; null if gene not found:
    protected Item[] geneItems;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapCcleMatrixConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (resolver == null) {
            resolver = IdResolverService.getIdResolverByOrganism(HUMAN_TAXON_ID);
        }
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        if (inputFile == null) {
            throw new RuntimeException("Input filename not defined. Set property 'input.file' in project.xml.");
        }
        if (outputAttribute == null) {
            throw new RuntimeException("Output attribute (of class 'DepMapCCLEData') not defined. Set property 'output.attribute' in project.xml.");
        }
        FileReader reader = new FileReader(new File(dataDir, inputFile));
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // parse header with gene symbols and IDs:
        String[] header = lineIter.next();
        if (!header[0].equals(headerStart)) {
            throw new RuntimeException("Unexpected item at start of input file");
        }
        parseGenes(header); // fills 'geneItems'

        // parse rows with cell line data:
        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }

            Item cellLine = createItem("CellLine");
            cellLine.setAttribute("DepMapID", line[0]);
            store(cellLine);

            for (int i = 1; i < header.length; i++) {
                if ((geneItems[i - 1] != null) && (!line[i].isEmpty())) {
                    Item dataItem = createItem("DepMapCCLEData");
                    dataItem.setReference("gene", geneItems[i - 1]);
                    dataItem.setReference("cellLine", cellLine);
                    dataItem.setAttribute(outputAttribute, line[i]);
                    store(dataItem);
                }
            }
        }
    }


    private void parseGenes(String[] header) throws ObjectStoreException {
        String regExp;
        if (idsArePrimary) {
            regExp = "[^ ]+ \\(\\d+\\)"; // expected format: "symbol (ID)"
        }
        else {
            regExp = "([^ ]+ \\()?ENSG\\d+\\)?"; // as above, or just "ID"
        }
        geneItems = new Item[header.length - 1]; // first column is row names
        // mapping: gene primary id. -> column index, supported by symbol?
        // (to detect "collisions" of column names mapping to the same primary id.)
        Map<String, SimpleEntry<Integer, Boolean>> primaryIdMap = new HashMap<>();

        for (int i = 1; i < header.length; i++) {
            String colName = header[i];
            if (!colName.matches(regExp)) {
                throw new RuntimeException("Unexpected format for column name: " + colName);
            }
            // use the ID (NCBI or Ensembl) from the header to resolve the gene,
            // but also look up the symbol as a consistency check:
            String[] parts = colName.split(" ");
            String geneId;
            Set<String> resolvedFromSymbol = Collections.emptySet();
            if (parts.length == 2) {
                geneId = parts[1].substring(1, parts[1].length() - 1);
                String geneSymbol = parts[0];
                resolvedFromSymbol = resolver.resolveId(HUMAN_TAXON_ID, "gene", geneSymbol);
                if (resolvedFromSymbol.isEmpty()) { // no match
                    LOG.info("Resolving gene '" + colName + "' by symbol: not found");
                }
                else if (resolvedFromSymbol.size() > 1) { // multiple matches
                    LOG.info("Resolving gene '" + colName + "' by symbol: " + resolvedFromSymbol.size() + " matches");
                }
            }
            else {
                geneId = parts[0];
            }
            boolean symbolMatches = false; // do gene resolution by symbol and by ID agree?
            if (idsArePrimary) {
                // check if gene with this NCBI ID (number) exists:
                if (!resolver.isPrimaryIdentifier(HUMAN_TAXON_ID, "gene", geneId)) {
                    LOG.warn("Resolving gene '" + colName + "' by primary identifier: not found - skipping");
                    continue;
                }
                if (!resolvedFromSymbol.isEmpty()) {
                    if (resolvedFromSymbol.contains(geneId)) {
                        symbolMatches = true;
                    }
                    else {
                        LOG.info("Resolving gene '" + colName + "': primary identifier/symbol mismatch");
                    }
                }
            }
            else { // IDs are secondary identifiers
                Set<String> resolvedFromId = resolver.resolveId(HUMAN_TAXON_ID, "gene", geneId);
                if (resolvedFromId.isEmpty()) { // no match
                    LOG.warn("Resolving gene '" + colName + "' by secondary identifier: not found - skipping");
                    continue;
                }
                else if (resolvedFromId.size() == 1) {
                    geneId = resolvedFromId.iterator().next();
                    if (!resolvedFromSymbol.isEmpty()) {
                        if (resolvedFromSymbol.contains(geneId)) {
                            symbolMatches = true;
                        }
                        else {
                            LOG.info("Resolving gene '" + colName + "': secondary identifier/symbol mismatch");
                        }
                    }
                }
                else { // multiple matches
                    LOG.info("Resolving gene '" + colName + "' by secondary identifier: " +
                             resolvedFromId.size() + " matches");
                    if (!resolvedFromSymbol.isEmpty()) { // try matching by symbol and identifier
                        resolvedFromId.retainAll(resolvedFromSymbol); // intersect
                    }
                    if (resolvedFromId.size() != 1) {
                        LOG.warn("Resolving gene '" + colName + "' by secondary identifier and symbol: failed (not found or not unique) - skipping");
                        continue;
                    }
                    symbolMatches = true;
                    geneId = resolvedFromId.iterator().next();
                }
            }

            // complication: multiple columns can map to the same primary identifier
            // -> if available, one where both ID and symbol match takes precedence
            SimpleEntry<Integer, Boolean> entry = primaryIdMap.get(geneId);
            if (entry == null) { // first match to this primary identifier
                Item gene = createItem("Gene");
                gene.setAttribute("primaryIdentifier", geneId);
                store(gene);
                geneItems[i - 1] = gene;
                primaryIdMap.put(geneId, new SimpleEntry(i - 1, symbolMatches));
            }
            else { // collision
                int otherIndex = entry.getKey();
                if (!entry.getValue() && symbolMatches) { // new entry takes precedence
                    LOG.warn("Multiple matches to gene with primary identifier '" + geneId + "': " +
                             header[otherIndex + 1] + " - skipping, " + colName + " - keeping");
                    // swap entries in gene list:
                    geneItems[i - 1] = geneItems[otherIndex];
                    geneItems[otherIndex] = null;
                    primaryIdMap.replace(geneId, new SimpleEntry(i - 1, true));
                }
                else { // keep old entry, skip current column
                    LOG.warn("Multiple matches to gene with primary identifier '" + geneId + "': " +
                             header[otherIndex + 1] + " - keeping, " + colName + " - skipping");
                }
            }
        }
    }


    public void setInputFile(String inputFile) {
        this.inputFile = inputFile;
    }

    public void setHeaderStart(String headerStart) {
        this.headerStart = headerStart;
    }

    public void setIdsArePrimary(String idsArePrimary) {
        this.idsArePrimary = Boolean.parseBoolean(idsArePrimary);
    }

    public void setOutputAttribute(String outputAttribute) {
        this.outputAttribute = outputAttribute;
    }
}
