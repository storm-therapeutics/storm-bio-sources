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

import java.util.*;


/**
 * 
 * @author
 */
public class DepmapCnvConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap Copy Number";
    private static final String DATA_SOURCE_NAME = "DepMap Public 19Q3";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String CN_CSV_FILE = "CCLE_gene_cn.csv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> cellLines = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(DepmapCnvConverter.class);

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapCnvConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {

        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processCopyNumber(new FileReader(files.get(CN_CSV_FILE)));

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processCopyNumber(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();
        ArrayList<String> genes = new ArrayList<String>();
        for(int i = 1; i < firstLine.length; i++) {
            String formattedGene = firstLine[i].split(" ")[0].trim();
            genes.add(formattedGene);
        }

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String cellLine = line[0];
            for(int i = 1; i < line.length; i++) {
                String cnvValue = line[i];
                String theGeneForThisItem = genes.get(i-1);
                Item CopyNumberItem;

                CopyNumberItem = createItem("DepMapCopyNumber");

                if(!cellLine.isEmpty()) {
                    CopyNumberItem.setReference("depMapID", getCellLine(cellLine));
                } else {
                    continue;
                }

                if(!theGeneForThisItem.isEmpty()) {
                    String geneId = getGeneId(theGeneForThisItem);

                    if (StringUtils.isEmpty(geneId)) {
                        continue;
                    }

                    CopyNumberItem.setReference("gene", geneId);
                } else {
                    continue;
                }

                if(!cnvValue.isEmpty() && StringUtils.isNumeric(cnvValue)) {
                    CopyNumberItem.setAttribute("value", cnvValue);
                } else {
                    continue;
                }

                store(CopyNumberItem);
                //cellLines.put(cellLine, CopyNumberItem.getIdentifier());
            }
        }
    }

    private String getGeneId(String primaryIdentifier) throws ObjectStoreException {
        String resolvedIdentifier = resolveGene(primaryIdentifier);
        if (StringUtils.isEmpty(resolvedIdentifier)) {
            return null;
        }
        String geneId = genes.get(resolvedIdentifier);
        if (geneId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", resolvedIdentifier);
            gene.setReference("organism", getOrganism(TAXON_ID));
            store(gene);
            geneId = gene.getIdentifier();
            genes.put(resolvedIdentifier, geneId);
        }
        return geneId;
    }

    private String resolveGene(String identifier) {
        String id = identifier;

        if (rslv != null && rslv.hasTaxon(TAXON_ID)) {
            int resCount = rslv.countResolutions(TAXON_ID, identifier);
            if (resCount != 1) {
                LOG.info("RESOLVER: failed to resolve gene to one identifier, ignoring gene: "
                        + identifier + " count: " + resCount + " Human identifier: "
                        + rslv.resolveId(TAXON_ID, identifier));
                return null;
            }
            id = rslv.resolveId(TAXON_ID, identifier).iterator().next();
        }
        return id;
    }

    public String getCellLine(String identifier) {
        String refId = cellLines.get(identifier);
        if (refId == null) {
            Item cl = createItem("CellLine");
            cl.setAttribute("DepMapID", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store cell line with DepMapID: " + identifier, e);
            }
            refId = cl.getIdentifier();
            cellLines.put(identifier, refId);
        }
        return refId;
    }
}
