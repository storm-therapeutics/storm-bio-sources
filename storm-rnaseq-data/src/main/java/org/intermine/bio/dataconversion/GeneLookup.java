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

    // mapping: Ensembl ID/gene symbol combination -> NCBI ID or 'null'
    private Map<String, String> resolutions = new HashMap<String, String>();

    protected IdResolver resolver;
    protected String taxon_id;
    private static final Logger LOG = Logger.getLogger(GeneLookup.class);

    public GeneLookup() {
        this("human");
    }

    public GeneLookup(String species) {
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
    }

    public void clear() {
        resolutions.clear();
    }

    public static boolean isValid(String s) {
        return (s != null) && !s.isEmpty() && !s.equals("NA");
    }

    /**
     * Look up a gene by identifier and/or symbol
     *
     * Return the NCBI (primary) ID for the gene, or 'null' if not found.
     * Create new entries in the lookup tables the first time each piece of information is looked up.
     */
    public String getGene(String ncbiId, String ensemblId, String symbol) throws ObjectStoreException {
        // NCBI ID given? - if yes, use it directly as primary identifier for the gene:
        if (isValid(ncbiId)) {
            return ncbiId;
        }

        // no NCBI (primary) ID given? - look it up based on Ensembl ID and/or symbol:
        boolean ensemblValid = isValid(ensemblId);
        boolean symbolValid = isValid(symbol);
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
            resolutions.put(combined, primaryId); // save for next time
            return primaryId;
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
