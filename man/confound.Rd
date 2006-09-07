\name{qb.confound}
\alias{qb.confound}
\alias{plot.qb.confound}
\title{Examine confounding of covariate with pseudomarkers.}
\description{
Covariates used in gene mapping may be correlated with covariates. These
routines examine the pattern of confounding.
}
\usage{
qb.confound(qbObject, covar = 1)
plot.qb.confound(x, ylim, main, ...)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{x}{Object of class \code{qb.confound}.}
  \item{covar}{Index to covariate}
  \item{ylim}{Limits for \code{y} (vertical) plotting axis.}
  \item{main}{Title for plot.}
  \item{\dots}{Additional parameters passed alone.}
}
\details{
This examines possible confounding between a covariate and pseudomarkers
across the genome. Confouding, evidenced by large correlation with a
marker, would raise suspicions about mapping in a genomic region,
unless of course the covariate is a marker in that region. Blue curves
are correlation with additive effect; red curves are correlation with
dominance effect. Dashed lines at 5 percent significance limits.
}
\value{
\code{qb.confound} returns a matrix with columns for:
  \item{coradd}{Correlation with additive pseudomarker effect.}
  \item{cordom}{Correlation with dominance pseudomarker effect (if F2).}
  \item{chr}{Chromosome identifier.}
The object inherits from \code{\link[qtl]{scanone}} objects.
}
\references{http://www.soph.uab.edu/qtlbim}
\author{Brian S. Yandell}
\seealso{\code{\link{qb.mcmc}}}
\examples{
\dontshow{
qb.load(cross, qbExample)


}
temp <- qb.confound(qbExample)
plot(temp)
}
\keyword{hplot}
