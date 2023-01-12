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

import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
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

    // mapping: NCBI/Entrez ID -> gene (InterMine Item)
    private Map<String, Item> geneItems = new HashMap<String, Item>();
    // mapping: Ensembl (secondary) ID -> NCBI/Entrez (primary) ID, or 'null'
    private Map<String, String> geneIDs = new HashMap<String, String>();
    // mapping: gene symbol -> NCBI/Entrez ID, or 'null'
    private Map<String, String> geneSymbols = new HashMap<String, String>();

    protected IdResolver resolver;
    protected String taxon_id;
    private static final Logger LOG = Logger.getLogger(GeneLookup.class);

    public GeneLookup() {
        GeneLookup("human");
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

    public clear() {
        geneItems.clear();
        geneIDs.clear();
        geneSymbols.clear();
    }

    public static boolean isValid(String s) {
        return (s != null) && !s.isEmpty() && !s.equals("NA");
    }

    /**
     * Look up a gene by identifier and/or symbol
     *
     * Return the Item for the gene, or 'null' if not found.
     * Create new entries in the lookup tables the first time each piece of information is looked up.
     */
    public Item getGene(String ncbiId, String ensemblId, String symbol) throws ObjectStoreException {
        boolean ensemblValid = isValid(ensemblId);
        boolean symbolValid = isValid(symbol);
        // remove version numbers (if any):
        if (ensemblValid)
            ensemblId = ensemblId.split("\\.")[0];
        if (symbolValid)
            symbol = symbol.split("\\.")[0];

        // NCBI ID given? - if yes, use it directly as primary identifier for the gene:
        if (isValid(ncbiId)) {
            Item gene = geneItems.get(ncbiId);
            if (gene == null) {
                gene = createItem("Gene");
                gene.setAttribute("primaryIdentifier", ncbiId);
                geneItems.put(ncbiId, gene);
            }
            // check/update other lookup tables if necessary:
            if (ensemblValid) {
                String ncbiId2 = geneIds.get(ensemblId);
                if (ncbiId2 == null) {
                    geneIds.put(ensemblId, ncbiId);
                }
                else if (!ncbiId2.equals(ncbiId)) {
                    LOG.warn("Conflicting primary IDs for Ensembl ID '" +
                             ensemblId + "': " + ncbiId2 + ", " + ncbiId);
                }
            }
            if (symbolValid) {
                String ncbiId2 = geneSymbols.get(symbol);
                if (ncbiId2 == null) {
                    geneSymbol.put(symbol, ncbiId);
                }
                else if (!ncbiId2.equals(ncbiId)) {
                    LOG.warn("Conflicting primary IDs for gene symbol '" +
                             symbol + "': " + ncbiId2 + ", " + ncbiId);
                }
            }
            return gene;
        }

        // no NCBI (primary) ID given? - look it up based on Ensembl ID and/or symbol:
        boolean ensemblExists = false;
        boolean symbolExists = false;
        String idFromEnsembl = null;
        String idFromSymbol = null;
        if (ensemblValid) {
            ensemblExists = geneIds.containsKey(ensemblId); // entry may exist with value 'null'
            if (ensemblExists) {
                idFromEnsembl = geneIds.get(ensemblId);
            }
        }
        if (symbolValid) {
            symbolExists = geneSymbols.containsKey(symbol); // entry may exist with value 'null'
            if (symbolExists) {
                idFromSymbol = geneSymbols.get(symbol);
            }
        }

        if ((idFromEnsembl != null) && (idFromSymbol != null)) {
            if (!idFromEnsembl.equals(idFromSymbol)) {
                LOG.warn("Conflicting primary IDs based on Ensembl ID '" + ensemblId +
                         "' and gene symbol '" + symbol + "': " + idFromEnsembl + ", " + idFromSymbol);
            }
            return geneItems.get(idFromEnsembl); // if in doubt, trust Ensembl ID over symbol
        }
        else if (idFromEnsembl != null) {
            if (symbolValid) { // update symbol lookup table (doesn't matter if entry was missing or 'null')
                geneSymbols.put(symbol, idFromEnsembl);
            }
            return geneItems.get(idFromEnsembl);
        }
        else if (idFromSymbol != null) {
            if (ensemblValid) { // update ID lookup table (doesn't matter if entry was missing or 'null')
                geneIDs.put(ensemblId, idFromSymbol);
            }
            return geneItems.get(idFromSymbol);
        }
        else if (ensemblExists && symbolExists) { // tried IDResolver before without success
            return null;
        }
        else { // information not previously seen - use IDResolver
            Set<String> resolvedFromEnsembl = null;
            Set<String> resolvedFromSymbol = null;
            if (ensemblValid && !ensemblExists) {
                resolvedFromEnsembl = resolver.resolveId(taxon_id, "gene", ensemblId);
                if (resolvedFromEnsembl.size() > 1) {
                    LOG.warn("Multiple gene matches for Ensembl ID '" + ensemblId + "': " +
                             ", ".join(resolvedFromEnsembl));
                }
            }
            if (symbolValid && !symbolExists) {
                resolvedFromSymbol = resolver.resolveId(taxon_id, "gene", symbol);
                if (resolvedFromSymbol.size() > 1) {
                    LOG.warn("Multiple gene matches for symbol '" + symbol + "': " +
                             ", ".join(resolvedFromSymbol));
                }
            }
            Set<String> resolved;
            if (!resolvedFromEnsembl.isEmpty()) {
                if (!resolvedFromSymbol.isEmpty()) {
                    resolvedFromEnsembl.retainAll(resolvedBySymbol); // intersect sets
                }
                resolved = resolvedFromEnsembl;
            }
            else if (!resolvedFromSymbol.isEmpty()) {
                resolved = resolvedFromSymbol;
            }
            if (resolved == null) { // no hits
                if (ensemblValid) {
                    LOG.warn("No gene found for Ensembl ID: " + ensemblId);
                    geneIds.put(ensemblId, null);
                }
                if (symbolValid) {
                    LOG.warn("No gene found for symbol: " + symbol);
                    geneSymbols.put(symbol, null);
                }
                return null;
            }
            else if (resolved.isEmpty()) { // conflict
                LOG.warn("No common gene match for Ensembl ID '" + ensemblId + "' and symbol '" + symbol "'");
                // @TODO: add 'null's to lookup tables in this case?
                return null;
            }
            else if (resolved.size() == 1) { // success
                String primaryId = resolvedFromEnsembl.iterator().next();
                Item gene = geneItems.get(primaryId);
                if (gene == null) {
                    gene = createItem("Gene");
                    gene.setAttribute("primaryIdentifier", primaryId);
                    geneItems.put(primaryId, gene);
                }
                return gene;
            }
            else { // multiple hits
                String message = "Failed to resolve multiple gene matches";
                if ()

                LOG.warn("Failed to resolve multiple gene matches for Ensembl ID " + ensemblId);
            }
        }
    }
}
