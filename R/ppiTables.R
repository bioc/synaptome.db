##### Use cases 5 Extract the PPIs for my list of genes #####
#' Prepare query for AllMutationAllPapers table
#'
#' @return tbl_lazy
#' @importFrom dplyr tbl select filter pull collect distinct
#' @keywords internal
getPPIQuery <- function() {
    #  gns<-get_dbconn() %>% dplyr::tbl("AllPpiAllPapers") %>%
    gns <- get_dbconn() %>%
        dplyr::tbl("PPI") %>%
        dplyr::select(
            # ID,#ppiID,
            A,
            B
        ) %>%
        distinct()
    # ,
    # method,
    # type,
    # taxID)
    return(gns)
}

#' Extract the PPIs for my list of genes defined by GeneID
#'
#' Get Protein-Protein interactions (PPIs) for the set of internal GeneIDs.
#' Function lookups for PPIs for specific GeneID and returns either “induced”
#' (all interaction for this GeneID) or “limited” (only interactions between
#' GeneIDs specified in the query) table of A and B interacting genes, where
#' A and B are respective GeneIDs.
#'
#' @param ids Gene IDs
#' @param type type of the PPI network should be either `induced` (for
#'     all the PPIs for specific genes, including external genes) or
#'     `limited` (for PPIs between the genes specified in the query).
#'     Type could be shortened to recognizable minimum like 'ind'
#'     or 'lim'.
#'
#' @return data.frame with interactors internal GeneID in columns A and B
#' @export
#'
#' @family {PPI functions}
#' @seealso [getPPIbyName()] and
#'     [getPPIbyEntrez()] to get
#'     PPI \code{data.frame},
#'     [getIGraphFromPPI()] to get igraph
#'     representation of the PPI \code{data.frame} and
#'     [getTableFromPPI()] to
#'     get interpretable
#'     table representation of the PPI \code{data.frame}.
#'
#' @md
#' @importFrom dplyr tbl select filter pull collect
#' @examples
#' t <- getPPIbyIDs(c(48, 585, 710), type = "limited") # (16 rows)
#' t <- getPPIbyIDs(c(48, 585, 710), type = "induced") # 306 rows
getPPIbyIDs <- function(ids, type = c("limited","induced")) {
    netType <- match.arg(type)
    gns <- switch(netType,
        induced = getPPIQuery() %>%
            dplyr::filter(A %in% ids | B %in% ids),
        limited = getPPIQuery() %>%
            dplyr::filter(A %in% ids & B %in% ids)
    )
    df <- gns %>% collect()
    return(df)
}

#' Extract the PPIs for my list of genes defined by Entrez IDs
#'
#' Get Protein-Protein interactions (PPIs) for the set of ID. Function
#' lookups for PPIs for the list of Entrez IDs and returns either “induced”
#' (all available interactions for these genes) or “limited” (only
#' interactions between genes specified in the query) table of
#' interacting genes A and B, where A and B are respective Internal IDs.
#'
#' @param entrez Entrez IDs
#' @param type type of the PPI network should be either `induced` (for
#'     all the PPIs for specific genes, including external genes) or
#'     `limited` (for PPIs between the genes specified in the query).
#'     Type could be shortened to recognizable minimum like 'ind'
#'     or 'lim'.
#'
#' @return data.frame with interactors internal GeneID in columns A and B
#' @export
#'
#' @family {PPI functions}
#' @seealso [getPPIbyName()] and
#'     [getPPIbyIDs()] to get
#'     PPI \code{data.frame},
#'     [getIGraphFromPPI()]
#'     to get igraph representation of the
#'     PPI \code{data.frame} and
#'     [getTableFromPPI()] to
#'     get interpretable
#'     table representation of the PPI \code{data.frame}.
#'
#' @md
#' @examples
#' t <- getPPIbyEntrez(c(1739, 1740, 1742, 1741), type = "ind")
getPPIbyEntrez <- function(entrez, type = c("limited","induced")) {
    ids <- getGeneIdByEntrez(entrez)
    df <- getPPIbyIDs(ids, type)
    return(df)
}

#' Extract the PPIs for my list of genes defined by Gene name
#'
#' Get Protein-Protein interactions (PPIs) for the set of gene names.
#' Function lookups for PPIs for the list of GeneIDs and returns either
#' “induced” (all interaction for this GeneID) or “limited” (only
#' interactions between GeneIDs specified in the query) table of
#' interacting genes A and B, where A and B are respective gene names.
#'
#' @param name Gene names
#' @param type type of the PPI network should be either `induced` (for
#'     all the PPIs for specific genes, including external genes) or
#'     `limited` (for PPIs between the genes specified in the query).
#'     Type could be shortened to recognizable minimum like 'ind'
#'     or 'lim'.
#'
#' @return data.frame with interactors internal GeneID in columns A and B
#' @export
#'
#' @family {PPI functions}
#' @seealso [getPPIbyEntrez()] and
#'     [getPPIbyIDs()] to get
#'     PPI \code{data.frame},
#'     [getIGraphFromPPI()] to get igraph
#'     representation of the PPI \code{data.frame} and
#'     [getTableFromPPI()] to
#'     get interpretable
#'     table representation of the PPI \code{data.frame}.
#'
#' @md
#' @examples
#' t <- getPPIbyName(
#'     c("CASK", "DLG4", "GRIN2A", "GRIN2B", "GRIN1"),
#'     type = "lim"
#' )
getPPIbyName <- function(name, type = c("limited","induced")) {
    ids <- getGeneIdByName(name)
    df <- getPPIbyIDs(ids, type)
    return(df)
}

#' Get table representation of the PPI.
#'
#' Combine information from PPI \code{data.frame} obtained with functions like
#' \code{\link{getPPIbyName}} or \code{\link{getPPIbyEntrez}} with information
#' about genes obtained from
#' \code{\link{getGenesByID}} to make  interpretable table representation.
#'
#' @param ppi PPI \code{data.frame} with columns A and B, obtained
#' from functions like \code{\link{getPPIbyName}}
#'
#' @return table with the following columns added for both
#' interactors (A and B):
#' \itemize{
#' \item{GeneID internal database ID}
#' \item{MGI MGI ID}
#' \item{HumanEntrez Human Entrez ID}
#' \item{MouseEntrez Mouse Entrez ID}
#' \item{HumanName Human gene name}
#' \item{MouseName Mouse gene name}
#' \item{RatEntrez Rat Entrez ID}
#' \item{Rat Name Rat gene name}
#' }

#' @export
#' @importFrom dplyr inner_join
#'
#' @family {PPI functions}
#' @family {DF functions}
#' @seealso [getPPIbyName()],
#'     [getPPIbyEntrez()]
#'     and [getPPIbyIDs()] to get
#'     PPI \code{data.frame},
#'     [getIGraphFromPPI()] to
#'     get igraph representation of the PPI \code{data.frame}.
#'
#' @md
#' @examples
#' tbl <- getTableFromPPI(getPPIbyIDs(c(48, 585, 710), type = "limited"))
getTableFromPPI <- function(ppi) {
    aTbl <- getGenesByID(ppi$A)
    names(aTbl) <- paste0(names(aTbl), ".A")
    bTbl <- getGenesByID(ppi$B)
    names(bTbl) <- paste0(names(bTbl), ".B")
    res <- ppi %>%
        dplyr::inner_join(aTbl, by = c("A" = "GeneID.A")) %>%
        dplyr::inner_join(bTbl, by = c("B" = "GeneID.B"))
    return(res)
}

#' Get Igraph representation of PPI
#'
#' Combine information from PPI \code{data.frame} obtained with functions like
#' \code{\link{getPPIbyName}} or \code{\link{getPPIbyEntrez}} with information
#' about genes obtained from
#' \code{\link{getGenesByID}} to make  interpretable undirected PPI graph in
#' \code{\link{igraph}} format. In this format network could be further
#' analysed and visualized by algorithms in \code{\link{igraph}} package.
#'
#' @param ppi PPI \code{data.frame} with columns A and B, obtaioed
#' from functions like \code{\link{getPPIbyName}}
#'
#' @return \code{\link{igraph}} object with specified PPI network.
#'
#' @family {PPI functions}
#' @family {Igraph functions}
#' @seealso [getPPIbyName()],
#'     [getPPIbyEntrez()]
#'     and [getPPIbyIDs()] to get
#'     PPI \code{data.frame},
#'     [getTableFromPPI()]
#'     to get interpretable table representation of
#'     the PPI \code{data.frame}.
#'
#' @export
#' @importFrom igraph graph_from_data_frame
#'
#' @md
#' @examples
#' library(igraph)
#' g <- getIGraphFromPPI(
#'     getPPIbyIDs(c(48, 129, 975, 4422, 5715, 5835), type = "lim")
#' )
#' plot(g, vertex.label = V(g)$RatName, vertex.size = 35)
getIGraphFromPPI <- function(ppi) {
    nids <- unique(c(ppi$A, ppi$B))
    nodes <- getGenesByID(nids)
    g <- graph_from_data_frame(ppi, directed = FALSE, vertices = nodes)
}
