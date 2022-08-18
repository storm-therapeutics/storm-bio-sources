package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

// import java.io.File;
// import java.io.FileReader;
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
public class DepmapCrisprGeneEffectConverter extends BioFileConverter // BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "CRISPR (Chronos) gene effect scores";
    private static final String DATA_SOURCE_NAME = "DepMap";

    // private static final String HUMAN_TAXON_ID = "9606";

    // private static final Logger LOG = Logger.getLogger(DepmapCrisprGeneEffectConverter.class);

    // protected IdResolver resolver;

    // private String[] geneIds;

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapCrisprGeneEffectConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    //     LOG.info("init base class");
    //     if (resolver == null) {
    //         resolver = IdResolverService.getIdResolverByOrganism(HUMAN_TAXON_ID);
    //     }
    //     LOG.info("init resolver");
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws Exception {
    // public void process(File dataDir) throws Exception {
        // LOG.info("start process");

        // FileReader reader = new FileReader(new File(dataDir, "CRISPR_gene_effect.csv"));
        // Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // LOG.info("parse file");
        // String[] header = lineIter.next();
        // if (!header[0].equals("DepMap_ID")) {
        //     throw new RuntimeException("Unexpected item at start of input file");
        // }

        // geneIds = new String[header.length - 1];

        // for (int i = 1; i < 11; i++) { // header.length; i++) {
        //     String[] parts = header[i].split(" ");
        //     if ((parts.length != 2) || (!parts[1].startsWith("(")) || (!parts[1].endsWith(")"))) {
        //         throw new RuntimeException("Unexpected item format in header: " + header[i]);
        //     }
        //     String geneSymbol = parts[0];
        //     String geneId = parts[1].substring(1, parts[1].length() - 1);

        //     // check if gene with this NCBI ID (number) exists:
        //     if (!resolver.isPrimaryIdentifier(HUMAN_TAXON_ID, "gene", geneId)) {
        //         LOG.info("No such gene found (by primary ID): " + header[i] + " - skipping");
        //         continue;
        //     }
        //     String resolvedId1 = resolver.resolveId(HUMAN_TAXON_ID, "gene", geneId).iterator().next();
        //     // check if NCBI ID (number) and gene symbol match:
        //     Set<String> resolvedIds = resolver.resolveId(HUMAN_TAXON_ID, "gene", geneSymbol);
        //     if (resolvedIds.isEmpty()) {
        //         LOG.info("No such gene found (by symbol): " + header[i] + " - skipping");
        //         continue;
        //     } else if (resolvedIds.size() > 1) {
        //         LOG.info(resolvedIds.size() + " matches for gene symbol: " + header[i]);
        //     }
        //     String resolvedId2 = resolvedIds.iterator().next();
        //     LOG.info(header[i] + ": " + resolvedId1 + "/" + resolvedId2);

            // if (!primaryId.equals(geneId)) {
            //     LOG.info("Primary ID/symbol mismatch for gene: " + header[i] + " - skipping");
            //     continue;
            // }
            // geneIds[i - 1] = resolvedId; // all is fine with this gene
        //  }

    //     while (lineIter.hasNext()) {
    //         String[] line = lineIter.next();
    //         if (line.length != header.length) {
    //             throw new RuntimeException("Unexpected number of items per line");
    //         }

    //         String cellLine = line[0];
    //         for (int i = 1; i < line.length; i++) {
    //             String score = line[i];
    //         }
    //     }

    }
}
