\name{plot.qb.pairloci}
\alias{qb.pairloci}
\alias{plot.qb.pairloci}
\alias{print.qb.pairloci}
\alias{summary.qb.pairloci}
\title{Summaries of epistatic pairs of loci.}
\description{
Summaries and detailed scatterplot showing all MCMC samples for
epistatic pairs for selected chromosomes.
}
\usage{
qb.pairloci(qbObject, chr)
## S3 methods.
plot.qb.pairloci(x, main, cex = 0.75, \dots )
print.qb.pairloci(x, \dots )
summary.qb.pairloci(object, \dots )
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.pairloci}.}
  \item{x}{Object of class \code{qb.pairloci}.}
  \item{chr}{Identifiers for one or two chromosomes.}
  \item{main}{Main title for plot.}
  \item{cex}{Character expansion of plot symbols.}
  \item{\dots}{Parameters to methods.}
}
\details{
Find pairs of loci in MCMC samples. Produce scatter plot with generic
\code{plot} or show numerical \code{summary}. The plot provides position
detail complementary to \code{\link{qb.multloci}} and
\code{\link{qb.scantwo}}.
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.qb}}, \code{\link{qb.scantwo}},
    \code{\link{qb.multloci}}}

\examples{
\dontshow{
qb.load(cross, qbExample)
  example(qb.mcmc)

}
temp <- qb.pairloci(qbExample, 1)
summary(temp)
plot(temp)
}
\keyword{models}
