\name{qb.slicetwo}
\alias{qb.slicetwo}
\alias{summary.qb.slicetwo}
\alias{plot.qb.slicetwo}
\alias{print.qb.slicetwo}
\title{Slices for epistatic pairs.}
\description{
These routines refine QTL positions for epistatic pairs and show plots
to reveal the nature of epistasis.
}
\usage{
qb.slicetwo(qbObject, chr, pos, type.scan = "2logBF", width = 10, \dots)
\method{summary}{qb.slicetwo}(object, \dots)
\method{print}{qb.slicetwo}(x, \dots)
\method{plot}{qb.slicetwo}(x, byrow = TRUE, figs, auto.par = TRUE, col, lty, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.slicetwo}.}
  \item{x}{Object of class \code{qb.slicetwo}.}
  \item{chr}{Chromosome vector.}
  \item{pos}{Position vector corresponding to \code{chr}.}
  \item{type.scan}{Type of profile scan; see \code{\link{qb.scanone}}.}
  \item{width}{Width of slice.}
  \item{byrow}{Arrange plots by row (for slides) if \code{TRUE}.}
  \item{figs}{Plot only selected figures. Full set of \code{c("profile",
      "effects", "cellmean", "effectplot")} is default.}
  \item{auto.par}{Automatic setting of plot parameters for multiple
    plots if \code{TRUE}.}
  \item{col,lty}{colors and line types for cellmean and interaction plots}
  \item{\dots}{Extra plot options.}
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link{qb.scantwo}}, \code{\link{qb.sliceone}}}
\examples{
data(qbExample)

temp <- qb.slicetwo(qbExample, chr = c(1,2), pos = c(45,12))
summary(temp)
plot(temp)
}
\keyword{hplot}
