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

import java.util.*;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.DataConverter;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Look up genes based on identifiers or symbols
 *
 * @author Hendrik Weisser
 */
public class GeneLookup
{
    private static final String HUMAN_TAXON_ID = "9606";
    private static final String MOUSE_TAXON_ID = "10090";

    // mapping: NCBI/Entrez (primary) ID -> gene (InterMine Item)
    private Map<String, Item> geneItems = new HashMap<String, Item>();
    // mapping: Ensembl ID/gene symbol combination -> gene (Item, also in 'geneItems') or null
    private Map<String, Item> resolutions = new HashMap<String, Item>();

    protected IdResolver resolver;
    protected String taxon_id;
    private DataConverter converter; // needed to create and store Items
    private static final Logger LOG = Logger.getLogger(GeneLookup.class);

    public GeneLookup(DataConverter converter) {
        this("human", converter);
    }

    public GeneLookup(String species, DataConverter converter) {
        if (species.equals("human")) {
            taxon_id = HUMAN_TAXON_ID;
        }
        else if (species.equals("mouse")) {
            taxon_id = MOUSE_TAXON_ID;
        }
        else {
            throw new IllegalArgumentException("Species not supported: " + species);
        }
        resolver = IdResolverService.getIdResolverByOrganism(taxon_id);
        this.converter = converter;
    }

    public void clear() {
        geneItems.clear();
        resolutions.clear();
    }

    public static boolean isValid(String s) {
        return (s != null) && !s.isEmpty() && !s.equals("NA");
    }

    /**
     * Write all gene Items in the lookup table to the database
     */
    public void storeAllGeneItems() throws ObjectStoreException {
        for (Item gene : geneItems.values()) {
            converter.store(gene);
        }
    }

    /**
     * Look up a gene by primary identifier (NCBI ID)
     *
     * If a gene with this primary identifier was looked up before, the corresponding stored Item is returned.
     * Otherwise a new Item will be created, added to the look-up table, and returned.
     */
    public Item getGene(String ncbiId) {
        // TODO: check if 'ncbiId' is a valid primary ID using IdResolver?
        Item gene = geneItems.get(ncbiId);
        if (gene == null) {
            gene = converter.createItem("Gene");
            gene.setAttribute("primaryIdentifier", ncbiId);
            geneItems.put(ncbiId, gene);
        }
        return gene;
    }

    /**
     * Look up a gene by identifier and/or symbol
     *
     * Return the Item for the gene, or 'null' if not found.
     * Create new entries in the lookup tables the first time each piece of information is looked up.
     */
    public Item getGene(String ncbiId, String ensemblId, String symbol) {
        // NCBI ID given? - if yes, use it directly as primary identifier for the gene:
        if (isValid(ncbiId)) {
            return getGene(ncbiId);
        }

        // no NCBI (primary) ID given? - look it up based on Ensembl ID and/or symbol:
        boolean ensemblValid = isValid(ensemblId);
        boolean symbolValid = isValid(symbol);
        if (!ensemblValid && !symbolValid)
            return null; // can't even write an informative error message in this case...

        String combined = "";
        if (ensemblValid) {
            ensemblId = ensemblId.split("\\.")[0]; // remove version number (if any)
            combined += ensemblId;
        }
        combined += "/";
        if (symbolValid) {
            symbol = symbol.split("\\.")[0]; // remove version number (if any)
            combined += symbol;
        }
        if (resolutions.containsKey(combined)) // resolved this combination of ID/symbol before?
            return resolutions.get(combined);

        List<String> given = new ArrayList<String>(); // what information was given (for error messages)
        Set<String> idsFromEnsembl = Collections.emptySet();
        Set<String> idsFromSymbol = Collections.emptySet();
        if (ensemblValid) {
            given.add("Ensembl ID '" + ensemblId + "'");
            idsFromEnsembl = resolver.resolveId(taxon_id, "gene", ensemblId);
        }
        if (symbolValid) {
            given.add("gene symbol '" + symbol + "'");
            idsFromSymbol = resolver.resolveId(taxon_id, "gene", symbol);
        }

        boolean possibleConflict = false; // conflicting results from Ensembl ID and symbol?
        Set<String> ids = idsFromEnsembl;
        if (ids.isEmpty()) {
            ids = idsFromSymbol;
        }
        else if (!idsFromSymbol.isEmpty()) {
            ids.retainAll(idsFromSymbol); // intersect sets
            possibleConflict = true;
        }

        if (ids.size() == 1) { // success!
            String primaryId = ids.iterator().next();
            Item gene = getGene(primaryId);
            resolutions.put(combined, gene); // save for next time
            return gene;
        }
        if (ids.isEmpty()) {
            if (possibleConflict)
                LOG.warn("Conflicting gene matches for " + String.join(" and ", given));
            else
                LOG.warn("No gene found for " + String.join(" and ", given));
        }
        else { // ids.size() > 1
            LOG.warn("Multiple gene matches for " + String.join(" and ", given));
        }

        resolutions.put(combined, null);
        return null;
    }
}
