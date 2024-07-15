##### Use cases 6 Extract the PPIs for specific compartment #####

#' Get full list of compartments
#'
#' Get list of synaptic compartments considered in the database.
#' Returns table with compartments: “Presynaptic”, “Postsynaptic”,
#' “Synaptosome”, and their internal IDs.
#'
#' @return data.frame described above
#' @family {Compartment functions}
#' @export
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' c <- getCompartments()
getCompartments <- function() {
    gns <- get_dbconn() %>%
        dplyr::tbl("Localisation") %>%
        dplyr::collect()
    return(gns)
}

#' Extract all genes found in compartment
#'
#' Get all genes annotated for specific compartment. Function
#' returns main information like internal Gene IDs, MGI ID,
#' Human Entrez ID, Human Gene Name, Mouse Entrez ID,
#' Mouse Gene Name, Rat Entrez ID, Rat Gene Name.
#'
#' @param compartmentID ID for specific compartment
#'
#' @return data.frame described in [getGenesByID()]
#' @family {Gene functions}
#' @family {Compartment functions}
#' @export
#' @md
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' G <- getAllGenes4Compartment(compartmentID = 1) # 5560 rows
getAllGenes4Compartment <- function(compartmentID) {
    idsC <- get_dbconn() %>%
        dplyr::tbl("FullGenePaper") %>%
        dplyr::filter(LocalisationID == compartmentID) %>%
        dplyr::select(GeneID) %>%
        dplyr::pull(GeneID) %>%
        unique()

    return(getGenesByID(idsC))
}

#' Select genes from the list that found in compartment
#'
#' Select all genes from your list annotated for specific
#' compartment. Should be used with
#' [findGenesByEntrez()]
#' or [findGenesByName()]
#' functions to obtain list of internal
#' IDs for your list of genes. Function lookups the Gene table
#' for specified localisation and returns main gene information
#' like internal Gene IDs, MGI ID, Human Entrez ID, Human Gene Name,
#' Mouse Entrez ID, Mouse Gene Name, Rat Entrez ID, Rat Gene Name.
#'
#' @param ids Gene IDs
#' @param compartmentID compartment ID
#'
#' @return data.frame described in [getGenesByID()]
#' @family {Gene functions}
#' @family {Compartment functions}
#' @export
#' @md
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' Genes <- getGenes4Compartment(c(1, 15, 156, 1500, 3000, 7000),
#'     compartmentID = 1
#' )
getGenes4Compartment <- function(ids, compartmentID) {
    idsC <- get_dbconn() %>%
        dplyr::tbl("FullGenePaper") %>%
        dplyr::filter(LocalisationID == compartmentID & GeneID %in% ids) %>%
        dplyr::select(GeneID) %>%
        dplyr::pull(GeneID) %>%
        unique()
    return(getGenesByID(idsC))
}

#' Prepare induces network for compartment
#'
#' @param ids gene IDs
#' @param compartmentID compartment ID
#'
#' @return  tbl_lazy
#' @importFrom dplyr tbl select filter pull collect
#' @keywords internal
getInducedPPI4Compartment <- function(ids, compartmentID) {
    cids <- getGenes4Compartment(ids, compartmentID)$GeneID
    aids <- getAllGenes4Compartment(compartmentID)$GeneID
    gns <- getPPIQuery() %>%
        dplyr::filter(A %in% cids | B %in% cids) %>%
        dplyr::filter(A %in% aids & B %in% aids)
    return(gns)
}

#' Prepare limited network for compartment
#'
#' @param ids Gene IDs
#' @param compartmentID compartment ID
#'
#' @return tbl_lazy
#' @importFrom dplyr tbl select filter pull collect
#' @keywords internal
getLimitedPPI4Compartment <- function(ids, compartmentID) {
    cids <- getGenes4Compartment(ids, compartmentID)$GeneID
    gns <- getPPIQuery() %>%
        dplyr::filter(A %in% cids & B %in% cids)
    return(gns)
}

#' Prepare induced or limited network for compartment
#'
#' Prepare PPI network for genes from your list annotated for
#' specific compartment. Should be used with findGenesByEntrez
#' or findGenesByName functions to obtain list of internal IDs
#' for your list of genes. Could be used with getAllGenes4Compartment
#' functions to obtain all genes belonging to respective compartment.
#' Function lookups the PPI table for gene IDs from the list and
#' returns “limited” or “induced” interactors GeneIDs from the
#' same compartment.
#'
#' @param ids internal gene IDs
#' @param compartmentID ID for specific compartment
#' @param type type of the PPI network should be either `induced` (for
#'     all the PPIs for specific genes, including external genes) or
#'     `limited` (for PPIs between the genes specified in the query).
#'     Type could be shortened to recognizable minimum like 'ind'
#'     or 'lim'.
#'
#' @return data.frame with interactors internal GeneID in columns A and B
#' @family {PPI functions}
#' @family {Compartment functions}
#' @export
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' ppi <- getPPIbyIDs4Compartment(c(1, 15, 156, 1500, 3000, 7000),
#'     compartmentID = 1, type = "induced"
#' ) # 201 rows
getPPIbyIDs4Compartment <- function(ids, compartmentID,
    type = c("limited","induced")) {
    netType <- match.arg(type)
    gns <- switch(netType,
        induced = getInducedPPI4Compartment(ids, compartmentID),
        limited = getLimitedPPI4Compartment(ids, compartmentID)
    )
    df <- gns %>% collect()
    return(df)
}
