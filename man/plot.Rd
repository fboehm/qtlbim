\name{plot.qb}
\alias{plot.qb}
\alias{print.qb}
\alias{summary.qb}
\title{Diagnostics plots for Bayesian interval mapping}
\description{
  Diagnostic plots highlight putative QTL loci and effects as well as
  providing graphical model assessment tools.
}
\usage{
\method{plot}{qb}(x, ask = dev.interactive(), verbose = TRUE, \dots)
\method{print}{qb}(x, \dots)
\method{summary}{qb}(object, cutoff = 1, \dots)
}
\arguments{
  \item{x}{An object of class \code{qb}.}
  \item{object}{An object of class \code{qb}.}
  \item{verbose}{Verbose summaries if \code{TRUE}.}
  \item{ask}{Ask before each plot if \code{TRUE}.}
  \item{cutoff}{Cutoff passed to \code{\link{qb.BayesFactor}}.}
  \item{\dots}{graphical parameters can be given as arguments to \code{\link[base]{plot}}}
}
\details{
  This generic \code{\link[base]{plot}} routine takes an object of class
  \code{qb} created by \code{\link{qb.mcmc}} and produces plots via
  calls to several other plot routines. The generic
  \code{\link[base]{summary}} produces a summary, while the generic
  \code{\link[base]{print}} passes through to \code{summary}.
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link[base]{plot}}, \code{\link[base]{print}},
  \code{\link[base]{summary}},
  \code{\link{qb.mcmc}},
  \code{\link{qb.coda}},
  \code{\link{qb.loci}},
  \code{\link{qb.BayesFactor}},
  \code{\link{qb.hpdone}},
  \code{\link{qb.epistasis}},
  \code{\link{qb.diag}} }

\examples{
\dontshow{
qb.load(cross, qbExample)
}
plot(qbExample)
summary(qbExample)
}
\keyword{hplot}
\keyword{models}
