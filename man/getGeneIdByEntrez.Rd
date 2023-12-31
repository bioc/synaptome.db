% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneInfo.R
\name{getGeneIdByEntrez}
\alias{getGeneIdByEntrez}
\title{Get list of GeneIDs corresponding to provided Entrez IDs.}
\usage{
getGeneIdByEntrez(entrez)
}
\arguments{
\item{entrez}{\code{vector} of Entres IDs. Function accepts both
integers and characters.}
}
\value{
\code{vector} of GeneID values.
}
\description{
Get internal GeneID values for set of gene Entrez IDs. Function
lookups for provided values in Human Entrez ID, Mouse Entrez ID and
Rat Entrez ID columns and returns obtained GeneIDs.
}
\examples{
t <- synaptome.db:::getGeneIdByEntrez(c(1742, 1741, 1739, 1740))
}
\keyword{internal}
