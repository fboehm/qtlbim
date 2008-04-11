\name{summary.qb.scanone}
\alias{summary.qb.scanone}
\alias{summary.qb.scantwo}
\alias{print.summary.qb.scantwo}
\title{Summary of qb.scanone or qb.scantwo object.}
\description{
Summary of a qb.scanone or qb.scantwo object.
}
\usage{
\method{summary}{qb.scanone}(object, chr, threshold = 0,
  sort = "no", n.qtl = 0.05, \dots)
\method{summary}{qb.scantwo}(object, chr, threshold = 0,
  sort = "no", which.pos = "upper", min.iter,
  refine = FALSE, width = 10, smooth = 3, n.qtl = 0.05,
  weight = c("sqrt","count","none"), \dots)
\method{print}{summary.qb.scantwo}(x, digits = 3, \dots)
}

\arguments{
  \item{object}{A \code{\link{qb.scanone}} or \code{\link{qb.scantwo}} object.}
  \item{x}{An object of class \code{qb.scantwo}.}
  \item{chr}{Chromosomes to include in summary (must be integers for now).}
  \item{threshold}{Threshold(s) for inclusion in summary (see below).}
  \item{sort}{Sort by selected column of \code{object} ("no" indicates
    sort by \code{chromosome}).}
  \item{which.pos}{Base position estimate on this summary for maximal
    statistics such as \code{LOD}.}
  \item{min.iter}{Minimum number of iterations included at each position
    (default gleaned from \code{object}).}
  \item{refine}{Refine estimates if \code{TRUE}.}
  \item{width}{Window width for refinement.}
  \item{smooth}{Degree of nearest neighbor smoothing to determine maxima.}
  \item{n.qtl}{Minimum number of estimated QTL per chromosome or
    chromosome pair.}
  \item{weight}{Weights to use for nearest neighbor
    smoothing. \code{sqrt} is square root of count per locus. Used only
    if \code{smooth} > 0.}
  \item{digits}{Significant digits to round with print.}
  \item{\dots}{Additional arguments for multiple linked loci (see details).}
}

\value{
Matrix with chromosome \code{chr}, estimated position \code{pos}
(or chromosome pairs\code{chr1} and \code{chr2} and two columns for
\code{pos1} and \code{pos2} in the case of \code{summary.qb.scantwo}) and
means or modes of each column of \code{object}. Means are weighted by
a smooth average of the number of MCMC sample iterations.
}
\details{
These summary method report estimates by
chromosome (or chromosome pair) at the maximum poster. Threshold can be
used to condense summary to a subset of chromosomes (or chromosome
pairs). Threshold is a 
vector with names corresponding to a subset of column names of
\code{object}. Positive threshold values select chromosomes where that
column average is above given value; negative threshold values select
chromosomes with mean value within that value of the maximum across
chromosomes. Thresholding is inclusive rather than exclusive.

It can be helpful to use summary.qb.scanone as an initial screen of
chromosomes worth a further look. Since marginal summaries can include
effects of multiple QTL and epistasis. Subsets based on 1-D scans can be
used for 2-D subsequent screens. See \code{\link{qb.demo}} for an
example.
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{
  \code{\link{qb.scanone}}, \code{\link{plot.qb.scanone}},
  \code{\link{qb.split.chr}}}
\examples{
data(qbExample)

temp <- qb.scanone(qbExample)
summary(temp, threshold = c(sum=15), sort = "sum")

temp <- qb.scantwo(qbExample)
summary(temp, threshold = c(upper=3), sort = "upper")
}
\keyword{utilities}
