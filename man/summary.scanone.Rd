\name{summary.qb.scanone}
\alias{summary.qb.scanone}
\alias{summary.qb.scantwo}
\title{Summary of qb.scanone or qb.scantwo object.}
\description{
Summary of a qb.scanone object.
}
\usage{
\method{summary}{qb.scanone}(object, chr, threshold = 0,
  sort = "no", digits = 3, smooth = 3, n.qtl = 0.05,
  min.iter, verbose = FALSE, \dots)
\method{summary}{qb.scantwo}(object, chr, threshold = 0,
  sort = "no", digits = 3, which.pos = "upper", min.iter,
  refine = FALSE, width = 10, smooth = 3, n.qtl = 0.05,
  verbose = FALSE, \dots)
}

\arguments{
  \item{object}{A \code{\link{qb.scanone}} object.}
  \item{chr}{Chromosomes to include in summary (must be integers for now).}
  \item{threshold}{Threshold(s) for inclusion in summary (see below).}
  \item{sort}{Sort by selected column of \code{object} ("no" indicates
    sort by \code{chromosome}).}
  \item{digits}{Significant digits to round summary.}
  \item{which.pos}{Base position estimate on this summary for maximal
    statistics such as \code{LOD}.}
  \item{min.iter}{Minimum number of iterations included at each position
    (default gleaned from \code{object}).}
  \item{refine}{Refine estimates if \code{TRUE}.}
  \item{width}{Window width for refinement.}
  \item{smooth}{Degree of nearest neighbor smoothing to determine maxima.}
  \item{n.qtl}{Minimum number of estimated QTL per chromosome or
    chromosome pair.}
  \item{verbose}{Show meta-information if \code{TRUE}.}
  \item{\dots}{Not used.}
}

\value{
Matrix with chromosome \code{chr}, estimated position \code{pos}
(or chromosome pairs\code{chr1} and \code{chr2} and two columns for
\code{pos1} and \code{pos2} in the case of \code{summary.qb.scantwo}) and
means or modes of each column of \code{object}. Means are weighted by number of
MCMC sample iterations.
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
used for 2-D subsequent screens. See \code{demo(qb.qb.scan.tour)} for an
example.
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link{qb.scanone}}, \code{\link{plot.qb.scanone}}}
\examples{
example(qb.scanone)
summary(temp, threshold = c(sum=5), sort = "sum")

example(qb.scantwo)
summary(temp, threshold = c(upper=3), sort = "upper")
}
\keyword{utilities}
