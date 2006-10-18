\name{plot.qb.multloci}
\alias{qb.multloci}
\alias{plot.qb.multloci}
\alias{print.qb.multloci}
\alias{summary.qb.multloci}
\title{Summaries of multiple loci on a chromosome.}
\description{
Summaries and detailed scatterplot showing all loci found in MCMC
samples for a chromosome.
}
\usage{
qb.multloci(qbObject, chr)
\method{plot}{qb.multloci}(x, amount = 0.75, cex, \dots )
\method{print}{qb.multloci}(x, \dots )
\method{summary}{qb.multloci}(object, \dots )
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.multloci}.}
  \item{x}{Object of class \code{qb.multloci}.}
  \item{chr}{Identifier for one chromosome.}
  \item{amount}{Amount to \code{\link[base]{jitter}} points.}
  \item{cex}{Character expansion of plot symbols.}
  \item{\dots}{Parameters to methods.}
}
\details{
Find multiple loci in MCMC samples for chromosome \code{chr}. Produce
scatter plot, density plots and histogram of counts with generic
\code{plot} or show numerical \code{summary}. The plot provides position
detail complementary to \code{\link{plot.qb.scantwo}}.
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.qb}}, \code{\link{qb.scantwo}}}

\examples{
\dontshow{
qb.load(cross, qbExample)
}
temp <- qb.multloci(qbExample, 1)
summary(temp)
plot(temp)
}
\keyword{models}
