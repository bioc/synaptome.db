##### Use cases 7 Extract the PPIs for specific brain region and specie #####

#' Get list of all Brain regions in the database
#'
#' Get full list of brain regions considered in the database. Returns
#' table with regions and their respective IDs.
#'
#' @return data.frame with following columns:
#' * ID: Brain region internal ID
#' * Name: name of the region
#' * Description: extended description of the region
#' * InterlexID
#' * ParentID: ID of the containing brain region
#' @export
#' @md
#' @family {BrainRegion functions}
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' t <- getBrainRegions()
getBrainRegions <- function() {
    gns <- get_dbconn() %>%
        dplyr::tbl("BrainRegion") %>%
        collect()
    return(gns)
}

#' Get all genes for brain region for particular specie
#'
#' Get all genes annotated for specific brain region for specific specie.
#' Function lookups Brain Region ID and Specie Tax ID columns and returns
#' main information like internal Gene IDs, MGI ID, Human Entrez ID,
#' Human Gene Name, Mouse Entrez ID, Mouse Gene Name,
#' Localisation (presynaptic, postsynaptic, synaptosomal),
#' PaperPMID and BrainRegion.
#'
#' @param brainRegion ID for specific brain region
#' @param taxID specie ID
#'
#' @return data.frame with the following columns:
#' * GeneID
#' * Localisation
#' * MGI
#' * HumanEntrez
#' * MouseEntrez
#' * HumanName
#' * MouseName
#' * PMID
#' * Paper
#' * Year
#' * SpeciesTaxID
#' * BrainRegion
#' @export
#' @md
#' @family {BrainRegion functions}
#' @family {BrainRegion Gene functions}
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' gns <- getAllGenes4BrainRegion(brainRegion = "Striatum", taxID = 10090)
#' head(gns)
getAllGenes4BrainRegion <- function(brainRegion, taxID) {
    idsC <- get_dbconn() %>%
        dplyr::tbl("FullGeneFullPaperFullRegion") %>%
        dplyr::filter(BrainRegion == brainRegion & SpeciesTaxID == taxID) %>%
        dplyr::select(
            GeneID, Localisation, MGI, HumanEntrez, MouseEntrez,
            HumanName, MouseName, PaperPMID, Paper, Year,
            SpeciesTaxID, BrainRegion
        ) %>%
        dplyr::rename(PMID = PaperPMID) %>%
        collect()
    return(idsC)
}

#' Select genes from the list that found in brain region of
#' particular specie
#'
#' Select genes from your list annotated for specific brain region. Should
#' be used with [findGenesByEntrez()]  or
#' [findGenesByName()] functions to obtain
#' list of internal IDs for your list of genes. Function lookups the
#' Gene table for specified localisation and returns main gene information
#' like internal Gene IDs, MGI ID, Human Entrez ID, Human Gene Name,
#' Mouse Entrez ID, Mouse Gene Name, PaperPMID, Localisation and BrainRegion.
#'
#' @param ids internal IDs for list of genes
#' @param brainRegion ID for the brain region of interest
#' @param taxID taxon ID specie tax ID
#'
#' @return data.frame with the following columns:
#' * GeneID
#' * Localisation
#' * MGI
#' * HumanEntrez
#' * MouseEntrez
#' * HumanName
#' * MouseName
#' * PMID
#' * Paper
#' * Year
#' * SpeciesTaxID
#' * BrainRegion
#' @export
#' @md
#' @family {BrainRegion functions}
#' @family {BrainRegion Gene functions}
#' @export
#' @importFrom dplyr tbl select filter pull collect "%>%"
#' @examples
#' Genes <- getGenes4BrainRegion(c(1, 15, 156, 1500, 3000, 7000),
#'     brainRegion = "Striatum", taxID = 10090
#' ) # 5 rows
getGenes4BrainRegion <- function(ids, brainRegion, taxID) {
    idsC <- get_dbconn() %>%
        dplyr::tbl("FullGeneFullPaperFullRegion") %>%
        dplyr::filter(
            BrainRegion == brainRegion &
                SpeciesTaxID == taxID &
                GeneID %in% ids
        ) %>%
        dplyr::select(
            GeneID, Localisation, MGI, HumanEntrez, MouseEntrez,
            HumanName, MouseName, PaperPMID, Paper, Year,
            SpeciesTaxID, BrainRegion
        ) %>%
        dplyr::rename(PMID = PaperPMID) %>%
        collect()
    return(idsC)
}

#' Prepare induces network for compartment
#'
#' @param brainRegion region ID
#' @param taxID taxon ID
#' @param ids gene IDs
#'
#' @return  tbl_lazy
#' @importFrom dplyr tbl select filter pull collect
#' @keywords internal
getInducedPPI4BrainRegion <- function(ids, brainRegion, taxID) {
    cids <- getGenes4BrainRegion(ids, brainRegion, taxID)$GeneID
    aids <- getAllGenes4BrainRegion(brainRegion, taxID)$GeneID
    gns <- getPPIQuery() %>%
        dplyr::filter(A %in% cids | B %in% cids) %>%
        dplyr::filter(A %in% aids & B %in% aids)
    return(gns)
}

#' Prepare limited network for compartment
#'
#' @param brainRegion region ID
#' @param taxID taxon ID
#' @param ids gene IDs
#'
#' @return  tbl_lazy
#' @importFrom dplyr tbl select filter pull collect
#' @keywords internal
getLimitedPPI4BrainRegion <- function(ids, brainRegion, taxID) {
    cids <- getGenes4BrainRegion(ids, brainRegion, taxID)$GeneID
    gns <- getPPIQuery() %>% dplyr::filter(A %in% ids & B %in% ids)
    return(gns)
}
#' Prepare induced or limited network for brain region
#'
#' Prepare PPI network for genes from your list annotated
#' for specific brain region. Should be used with
#' [findGenesByEntrez()])
#' or [findGenesByName()]
#' functions to obtain list of internal IDs for your list of genes.
#' Could be used with getAllGenes4BrainRegion functions to obtain
#' all genes belonging to respective brain region. Function lookups
#' the PPI table for gene IDs from the list and returns “limited”
#' or “induced” interactors GeneIDs for the specified region.
#'
#' @param ids gene IDs
#' @param brainRegion region ID
#' @param taxID taxon ID
#' @param type type of the PPI network should be either `induced` (for
#'     all the PPIs for specific genes, including external genes) or
#'     `limited` (for PPIs between the genes specified in the query).
#'     Type could be shortened to recognizable minimum like 'ind'
#'     or 'lim'.
#'
#' @return data.frame with interactors internal GeneID in columns A and B
#' @family {PPI functions}
#' @family {BrainRegion functions}
#' @md
#' @export
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' # getting all genes for mouse Striatum
#' gns <- getAllGenes4BrainRegion(brainRegion = "Striatum", taxID = 10090)
#' head(gns)
#'
#' # getting full PPI network for postsynaptic compartment
#' ppi <- getPPIbyIDs4BrainRegion(
#'     gns$GeneID,
#'     brainRegion = "Striatum",
#'     taxID = 10090,
#'     type = "limited"
#' )
#' head(ppi)
getPPIbyIDs4BrainRegion <- function(ids, brainRegion, taxID,
    type = c("limited","induced")) {
    netType <- match.arg(type)
    gns <- switch(netType,
        induced = getInducedPPI4BrainRegion(ids, brainRegion, taxID),
        limited = getLimitedPPI4BrainRegion(ids, brainRegion, taxID)
    )
    df <- gns %>% collect()
    return(df)
}
