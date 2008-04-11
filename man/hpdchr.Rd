\name{qb.hpdchr}
\alias{qb.hpdchr}
\title{Find highest probability density (HPD) region by Chromosome.}
\description{
Find area under posterior by chromosome for genome-wide HPD region.
}
\usage{
qb.hpdchr(qbObject, level = 0.5, height, hpd, chr, smooth = 3)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{level}{Value between 0 and 1 of HPD coverage. Ignored if
    \code{height} is supplied.}
  \item{height}{Height of posterior corresponding to
    \code{level}. Determined from \code{level} and \code{hpd} by default.}
  \item{hpd}{Object of class \code{qb.hpdone}. Created internally be default.}
  \item{chr}{Chromosomes to include; default determined by genome-wide HPD region.}
  \item{smooth}{Degree of smoothing.}
}
\details{
Determine 100*\code{level} percent HPD coverage by chromosome for a
given \code{height} threshold.
}
\value{
  A List with the following elements:
  \item{hpd.height}{Height of HPD threshold. Name of \code{hpd.height}
    is the \code{level}.}
  \item{chr.posterior}{Posterior probability by chromosome as percent.}
}
\references{http://www.qtlbim.org}
\author{Brian S. Yandell}
\seealso{\code{\link{qb.scanone}}}
\examples{
data(qbExample)

qb.hpdchr(qbExample)
}
\keyword{utilities}
