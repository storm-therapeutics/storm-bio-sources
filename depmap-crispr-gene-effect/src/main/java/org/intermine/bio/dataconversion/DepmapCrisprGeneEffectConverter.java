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
import java.util.Iterator;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 * Read CRISPR gene effect data from DepMap
 * @author Hendrik Weisser
 */
public class DepmapCrisprGeneEffectConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "CRISPR (Chronos) gene effect scores";
    private static final String DATA_SOURCE_NAME = "DepMap";
    private static final String INPUT_FILE_NAME = "CRISPR_gene_effect.csv";

    private static final String HUMAN_TAXON_ID = "9606";

    private static final Logger LOG = Logger.getLogger(DepmapCrisprGeneEffectConverter.class);

    protected IdResolver resolver;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapCrisprGeneEffectConverter(ItemWriter writer, Model model) {
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
        FileReader reader = new FileReader(new File(dataDir, INPUT_FILE_NAME));
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        String[] header = lineIter.next();
        if (!header[0].equals("DepMap_ID")) {
            throw new RuntimeException("Unexpected item at start of input file");
        }

        Item[] geneItems = new Item[header.length - 1];

        int maxIter = header.length;
        for (int i = 1; i < maxIter; i++) {
            String[] parts = header[i].split(" ");
            if ((parts.length != 2) || (!parts[1].startsWith("(")) || (!parts[1].endsWith(")"))) {
                throw new RuntimeException("Unexpected item format in header: " + header[i]);
            }
            String geneSymbol = parts[0];
            String genePrimaryId = parts[1].substring(1, parts[1].length() - 1);

            // check if gene with this NCBI ID (number) exists:
            if (!resolver.isPrimaryIdentifier(HUMAN_TAXON_ID, "gene", genePrimaryId)) {
                LOG.warn("No such gene found (by primary ID): " + header[i] + " - skipping");
                continue;
            }
            // check if NCBI ID (number) and gene symbol match:
            Set<String> resolvedIds = resolver.resolveId(HUMAN_TAXON_ID, "gene", geneSymbol);
            if (resolvedIds.isEmpty()) { // no match
                LOG.warn("No such gene found (by symbol): " + header[i]);
                // @TODO: can we find out what the "other" symbol is (that matches the ID)?
            }
            else if (resolvedIds.size() > 1) { // multiple matches
                LOG.warn(resolvedIds.size() + " matches for gene symbol: " + header[i]);
            }
            else { // one match
                String primaryId = resolvedIds.iterator().next();
                if (!primaryId.equals(genePrimaryId)) {
                    LOG.warn("Primary ID/symbol mismatch for gene: " + header[i]);
                }
            }

            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", genePrimaryId);
            // @TODO: what happens if we set the symbol and it doesn't match (see case above)?
            store(gene);
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
                    Item effect = createItem("DepMapCRISPRGeneEffect");
                    effect.setReference("gene", geneItems[i - 1]);
                    effect.setReference("cellLine", cellLine);
                    effect.setAttribute("geneEffectScore", line[i]);
                    store(effect);
                }
            }
        }
    }
}
