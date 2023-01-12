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
    // store results of IDResolver queries:
    // mapping: Ensembl (secondary) ID -> NCBI/Entrez (primary) IDs, or 'null'
    private Map<String, Set<String>> geneIds = new HashMap<String, Set<String>>();
    // mapping: gene symbol -> NCBI/Entrez IDs, or 'null'
    private Map<String, Set<String>> geneSymbols = new HashMap<String, Set<String>>();

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
        geneIds.clear();
        geneSymbols.clear();
    }

    public static boolean isValid(String s) {
        return (s != null) && !s.isEmpty() && !s.equals("NA");
    }

    public Item getGene(String ncbiId) throws ObjectStoreException {
        Item gene = geneItems.get(ncbiId);
        if (gene == null) {
            gene = converter.createItem("Gene");
            gene.setAttribute("primaryIdentifier", ncbiId);
            converter.store(gene);
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
    public Item getGene(String ncbiId, String ensemblId, String symbol) throws ObjectStoreException {
        // NCBI ID given? - if yes, use it directly as primary identifier for the gene:
        if (isValid(ncbiId)) {
            return getGene(ncbiId);
        }

        // no NCBI (primary) ID given? - look it up based on Ensembl ID and/or symbol:
        List<String> given = new ArrayList<String>(); // what information was given (for error messages)
        Set<String> idsFromEnsembl = Collections.emptySet();
        Set<String> idsFromSymbol = Collections.emptySet();
        if (isValid(ensemblId)) {
            ensemblId = ensemblId.split("\\.")[0]; // remove version number (if any)
            given.add("Ensembl ID '" + ensemblId + "'");
            idsFromEnsembl = geneIds.get(ensemblId);
            if (idsFromEnsembl == null) { // not looked up previously - do it now
                idsFromEnsembl = resolver.resolveId(taxon_id, "gene", ensemblId);
                geneIds.put(ensemblId, idsFromEnsembl);
            }
        }
        if (isValid(symbol)) {
            symbol = symbol.split("\\.")[0]; // remove version number (if any)
            given.add("gene symbol '" + symbol + "'");
            idsFromSymbol = geneSymbols.get(symbol);
            if (idsFromSymbol == null) { // not looked up previously - do it now
                idsFromSymbol = resolver.resolveId(taxon_id, "gene", symbol);
                geneIds.put(symbol, idsFromSymbol);
            }
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
            return getGene(primaryId);
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

        return null;
    }
}
