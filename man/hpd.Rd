\name{qb.hpdone}
\alias{qb.hpdone}
\alias{summary.qb.hpdone}
\alias{print.qb.hpdone}
\alias{plot.qb.hpdone}
\title{Highest probability density (HPD) region.}
\description{
Determine HPD region across genome, including position of posterior mode.
}
\usage{
qb.hpdone(qbObject, level = 0.5, profile = "2logBF",
  effects = "cellmean", scan = "sum", chr, smooth = 3, \dots)
\method{summary}{qb.hpdone}(object, chr, digits = 3, \dots)
\method{print}{qb.hpdone}(x, \dots)
\method{plot}{qb.hpdone}(x, chr, smooth = 3, \dots)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.hpdone}.}
  \item{x}{Object of class \code{qb.hpdone}.}
  \item{level}{Value between 0 and 1 of HPD coverage.}
  \item{scan}{Elements to scan; usually one of \code{"sum"},
    \code{"mean"}, \code{"epistasis"}, \code{"GxE"}.}
  \item{smooth}{Degree of smoothing.}
  \item{chr}{Chromosomes to include; default determined by HPD region.}
  \item{effects}{Effects are \code{"cellmean"} for means by genotype;
    \code{"estimate"} for estimates of Cockerham main effects.}
  \item{profile}{Objective profile for plot; default is \code{"2logBF"};
    other choices found in option \code{type} for
    \code{\link{qb.scanone}}.}
  \item{digits}{Number of digits for \code{\link[base]{round}}.}
  \item{\dots}{Extra parameters passed along to plot.}
}
\details{
Determine 100*\code{level} percent HPD region. Subset chromosomes based
on HPD region. Create genome scans for \code{profile} and \code{effects}.
}
\value{
  \code{qb.hpdone} is a list with a \code{hpd.region} summary matrix and
  \code{\link{qb.scanone}} objects for the \code{profile} and
  \code{effects}. A summary of a \code{qb.hpdone} object yields a matrix
  with columns for
  \item{chr}{chromosome number}
  \item{n.qtl}{estimated number of QTL on chromosome}
  \item{pos}{estimated position of QTL}
  \item{lo.nn\%}{lower nn\% HPD limit}
  \item{hi.nn\%}{upper nn\% HPD limit}
  \item{profile}{Peak of profile, identifed by the profile type.}
  \item{effects}{Columns for the effects, appropriately labeled.}
}
\references{http://www.ssg.uab.edu/qtlbim}
\author{Brian S. Yandell}
\seealso{\code{\link{qb.scanone}}}
\examples{
\dontshow{
qb.load(cross, qbExample)


}
temp <- qb.hpdone(qbExample)
summary(temp)
plot(temp)
}
\keyword{hplot}
