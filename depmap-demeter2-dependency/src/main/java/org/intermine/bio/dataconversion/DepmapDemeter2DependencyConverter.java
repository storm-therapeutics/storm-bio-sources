package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2018 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 *
 * @author
 */
public class DepmapDemeter2DependencyConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap RNAi Gene Dependency (DEMETER2)";
    private static final String DATA_SOURCE_NAME = "DepMap";
    private static final String D2_CSV_FILE = "D2_combined_gene_dep_scores.csv";
    private static final String TAXON_ID = "9606"; // Human Taxon ID

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(DepmapDemeter2DependencyConverter.class);

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();
    private Map<String, String> cellLinesMap = new HashMap<String, String>();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapDemeter2DependencyConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {
        FileReader reader = new FileReader(new File(dataDir, D2_CSV_FILE));
        Iterator<String[]> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header (cell line names):
        String[] firstLine = (String[]) lineIter.next();
        Item[] cellLineItems = new Item[firstLine.length - 1];
        for (int i = 1; i < firstLine.length; i++) {
            String formattedCL = firstLine[i].trim().replaceAll("\"", "");
            Item cl = createItem("CellLine");
            cl.setAttribute("CCLEname", formattedCL);
            store(cl);
            cellLineItems[i - 1] = cl;
        }

        while (lineIter.hasNext()) {
            String[] line = lineIter.next();
            // first column contains gene info:
            // format: "symbol (ID)", or "symbol1&symbol2 (ID1&ID2)" etc.
            String geneId = line[0].split(" ")[1]; // resolve by ID
            geneId = geneId.substring(1, geneId.length() - 1);
            if (geneId.contains("&")) // multiple genes - skip
                continue;
            if (!rslv.isPrimaryIdentifier(TAXON_ID, "gene", geneId)) {
                LOG.warn("Resolving gene '" + line[0] + "' by primary identifier: not found - skipping");
                continue;
            }
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", geneId);
            store(gene);
            // following columns contain dependency (gene effect) values:
            for (int i = 1; i < line.length; i++) {
                String dependencyValue = line[i];
                if (!dependencyValue.isEmpty() && !dependencyValue.equals("NA")) {
                    Item dataItem = createItem("DepMapCCLEData");
                    dataItem.setReference("gene", gene);
                    dataItem.setReference("cellLine", cellLineItems[i - 1]);
                    dataItem.setAttribute("rnaiGeneEffect", dependencyValue);
                    store(dataItem);
                }
            }
        }
    }
}
