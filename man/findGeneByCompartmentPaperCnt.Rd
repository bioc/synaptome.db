% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneInfo.R
\name{findGeneByCompartmentPaperCnt}
\alias{findGeneByCompartmentPaperCnt}
\title{Get gene table of frequently found genes within compartments}
\usage{
findGeneByCompartmentPaperCnt(cnt = 1)
}
\arguments{
\item{cnt}{minimal number of times mentioned gene}
}
\value{
\code{data.frame} with 9 columns: 8 specified in
        \code{\link{getGenesByID}} and \code{Npmid} column for the paper
        count.
}
\description{
Get gene table and paper count for genes mentioned \code{cnt}
or more times in different compartment-paper pairs.
}
\examples{
cntT <- findGeneByPaperCnt(47)
head(cntT)
}
\seealso{
getGenesByID

Other {Gene functions}: 
\code{\link{findGeneByPaperCnt}()},
\code{\link{findGeneByPapers}()},
\code{\link{findGenesByEntrez}()},
\code{\link{findGenesByName}()},
\code{\link{getAllGenes4Compartment}()},
\code{\link{getGenes4Compartment}()},
\code{\link{getGenesByID}()}
}
\concept{{Gene functions}}
