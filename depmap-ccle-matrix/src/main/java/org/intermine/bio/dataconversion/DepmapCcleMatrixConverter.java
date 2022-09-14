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
import java.util.Iterator;
import java.util.Set;

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
        String[] header = lineIter.next();
        if (!header[0].equals(headerStart)) {
            throw new RuntimeException("Unexpected item at start of input file");
        }

        Item[] geneItems = new Item[header.length - 1];

        int maxIter = header.length;
        for (int i = 1; i < maxIter; i++) {
            Item gene = parseGene(header[i]); // null if gene not found
            geneItems[i - 1] = gene;
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            if (line.length != header.length) {
                throw new RuntimeException("Unexpected number of items per line");
            }

            Item cellLine = createItem("CellLine");
            cellLine.setAttribute("DepMapID", line[0]);
            store(cellLine);

            for (int i = 1; i < maxIter; i++) {
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


    private Item parseGene(String headerItem) throws ObjectStoreException {
        String regExp1 = "[^ ]+ \\(\\d+\\)"; // expected format: "symbol (ID)"
        String regExp2 = "([^ ]+ \\()?ENSG\\d+\\)?"; // as above, or just ID
        if ((idsArePrimary && !headerItem.matches(regExp1)) ||
            (!idsArePrimary && !headerItem.matches(regExp2))) {
            throw new RuntimeException("Unexpected item format in header: " + headerItem);
        }
        // use the ID (NCBI or Ensembl) from the header to resolve the gene,
        // but also look up the symbol as a consistency check:
        String[] parts = headerItem.split(" ");
        String geneId;
        Set<String> resolvedFromSymbol = Collections.emptySet();
        if (parts.length == 2) {
            geneId = parts[1].substring(1, parts[1].length() - 1);
            String geneSymbol = parts[0];
            resolvedFromSymbol = resolver.resolveId(HUMAN_TAXON_ID, "gene", geneSymbol);
            if (resolvedFromSymbol.isEmpty()) { // no match
                LOG.info("Resolving gene '" + headerItem + "' by symbol: not found");
            }
            else if (resolvedFromSymbol.size() > 1) { // multiple matches
                LOG.info("Resolving gene '" + headerItem + "' by symbol: " + resolvedFromSymbol.size() + " matches");
            }
        }
        else {
            geneId = parts[0];
        }

        if (idsArePrimary) {
            // check if gene with this NCBI ID (number) exists:
            if (!resolver.isPrimaryIdentifier(HUMAN_TAXON_ID, "gene", geneId)) {
                LOG.warn("Resolving gene '" + headerItem + "' by primary identifier: not found - skipping");
                return null;
            }
            if (!resolvedFromSymbol.isEmpty() && !resolvedFromSymbol.contains(geneId)) {
                LOG.info("Resolving gene '" + headerItem + "': primary identifier/symbol mismatch");
            }
        }
        else {
            Set<String> resolvedFromId = resolver.resolveId(HUMAN_TAXON_ID, "gene", geneId);
            if (resolvedFromId.isEmpty()) { // no match
                LOG.warn("Resolving gene '" + headerItem + "' by secondary identifier: not found - skipping");
                return null;
            }
            else if (resolvedFromId.size() == 1) {
                geneId = resolvedFromId.iterator().next();
                if (!resolvedFromSymbol.isEmpty() && !resolvedFromSymbol.contains(geneId)) {
                    LOG.info("Resolving gene '" + headerItem + "': secondary identifier/symbol mismatch");
                }
            }
            else if (resolvedFromId.size() > 1) { // multiple matches
                LOG.info("Resolving gene '" + headerItem + "' by secondary identifier: " + resolvedFromId.size() + " matches");
                if (!resolvedFromSymbol.isEmpty()) { // try matching by symbol and identifier
                    resolvedFromId.retainAll(resolvedFromSymbol); // intersect
                }
                if (resolvedFromId.size() != 1) {
                    LOG.warn("Resolving gene '" + headerItem + "' by secondary identifier and symbol: failed (not found or not unique) - skipping");
                    return null;
                }
                geneId = resolvedFromId.iterator().next();
            }
        }

        Item gene = createItem("Gene");
        gene.setAttribute("primaryIdentifier", geneId);
        store(gene);
        return gene;
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
