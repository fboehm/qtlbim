\name{subset.qb}
\alias{subset.qb}
\title{Subsetting Bayesian interval mapping data}
\description{
Subset Bayesian interval mapping iterations on number of QTL and/or
chromosome pattern of QTL, using exact match or inclusive subsetting.
}
\usage{
subset.qb(x, nqtl=1, pattern=NULL, exact=FALSE, chr,
  region, offset = TRUE, restrict.pair = TRUE, \dots)
}
\arguments{
  \item{x}{ object of class \code{qb} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{chr}{subset of chromosomes to plot (numerical indices or
    chromosome names)}
  \item{region}{list containing \code{chr}, and \code{start} and \code{end}
    positions, for regions to include}
  \item{offset}{indicates that \code{start} and \code{end} are in cM
    position if \code{TRUE}; otherwise they are in distance from first marker}
  \item{restrict.pair}{Restrict \code{chr} selection to
    epistatic pairs with both in selected subset.}
  \item{\dots}{additional arguments to \code{\link[base]{subset}}}
}
\details{
Subset to include only iterations with at least \code{nqtl} number of
QTL and at least the \code{pattern} across chromosomes. \code{pattern}
  is a vector of chromosome identices, with repeats for
  multiple linked QTL on a chromosome. If \code{exact=FALSE}, then all
  iterations with at least the given \code{pattern} and \code{nqtl} are
  included. \code{nqtl} will be reset to \code{length(pattern)} if it is
  smaller than this value. Note that \code{pattern} should be number
  codes corresponding to those used in the \code{x} object. At
  present, chromosome names are not allowed. Further subsets to only
  include QTL from these iterations that are on chromosomes \code{chr}.
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link[qtl]{read.cross}}}

\examples{
\dontshow{
qb.load(cross, qbExample)
}
## Subset to chr 1,3,5.
qbSubset <- subset(qbExample, chr = c(1,3,5))

## Subset to chr 1,3,5 but restrict to only pairs with 1,3,5.
qbSubset <- subset(qbExample, chr = c(1,3,5), restrict.pair = TRUE)

## Subset to chr 1,3,5, and to within 10 cM of QTL on chr 3,5.
qbSubset <- subset(qbExample, chr = c(1,3,5),
  region = data.frame(chr = c(3,5), start = c(2,5), end = c(22,25)))
}
\keyword{utilities}
