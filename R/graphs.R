#' Utility function to create network from
#' \code{\link[synaptome.db]{synaptome.db}} data
#'
#' @param entrez vector of EntrezIDs for network vertices
#' @param LCC if TRUE only largest connected component is returned
#' @param simplify if TRUE loops and multiple edges will be removed
#'
#' @return network defined by the gene table
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' geneTable<-getAllGenes4Compartment(cid)
#' gg<-graphFromSynaptomeByEntrez(geneTable$HumanEntrez)
graphFromSynaptomeByEntrez<-function(entrez,LCC=TRUE,simplify=TRUE){
    geneTable<-findGenesByEntrez(entrez)
    gg<-graphFromSynaptomeGeneTable(geneTable,LCC=LCC,simplify=simplify)
    return(gg)
}

#' Utility function to create network from
#' \code{\link{synaptome.db}} data
#'
#' @param geneTable data.frame described in
#'        \code{\link{getGenesByID}}
#' @param LCC if TRUE only largest connected component is returned
#' @param simplify if TRUE loops and multiple edges will be removed
#' @param type type of the PPI network should be either `induced` (for
#'     all the PPIs for specific genes, including external genes) or
#'     `limited` (for PPIs between the genes specified in the query).
#'     Type could be shortened to recognizable minimum like 'ind'
#'     or 'lim'.
#'
#' @return network defined by the gene table
#' @export
#' @seealso [getPPIbyIDs()]
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic', getCompartments()$Name)
#' geneTable<-getAllGenes4Compartment(cid)
#' gg<-graphFromSynaptomeGeneTable(geneTable)
graphFromSynaptomeGeneTable<-function(geneTable,LCC=FALSE,simplify=FALSE,
                                      type = c("induced", "limited")){
    p<-getPPIbyIDs(geneTable$GeneID, type = type)
    aidx<-match(p$A, geneTable$GeneID)
    bidx<-match(p$B, geneTable$GeneID)
    # TODO: uncomment once BioNAR>= 1.3.7 in Bioconductor
    gg<-buildNetwork(data.frame(A=geneTable$HumanEntrez[aidx],
                                B=geneTable$HumanEntrez[bidx]),
                    LCC=LCC,simplify=simplify)
    return(gg)
}

#' @importFrom igraph simplify
buildNetwork<-function(ff, kw=NA,LCC=TRUE,simplify=TRUE){
    #--- build raw graph
    GG <- graph_from_data_frame(ff[, seq_len(2)], directed=FALSE)
    if(simplify){
        GG <- simplify(GG, remove.multiple=TRUE, remove.loops=TRUE)
    }
    if(LCC){
        bnr <- system.file(package = "BioNAR")
        if(nzchar(bnr)){
            GG <- BioNAR::findLCC(GG)
        }else{
        warning("LCC calculations not implemented yet.\n",
                "Install BioNAR and use BioNAR::findLCC.\n")
        }
    }
    return(GG)
}
