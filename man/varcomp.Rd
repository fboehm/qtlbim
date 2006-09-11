\name{qb.varcomp}
\alias{qb.varcomp}
\alias{summary.qb.varcomp}
\alias{print.qb.varcomp}
\alias{plot.qb.varcomp}
\title{Variance components for Bayesian multiple QTL}
\description{
These routines extract and summarize variance components for Bayesian
multiple QTL. Variance components are averaged over genome
loci. Covariates and GxE may be included.
}
\usage{
qb.varcomp(qbObject, scan, aggregate = TRUE)
summary.qb.varcomp(object, ...)
print.qb.varcomp(x, ...)
plot.qb.varcomp(x, log = TRUE, percent = 5, cex, ...)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.varcomp}.}
  \item{x}{Object of class \code{qb.varcomp}.}
  \item{scan}{Aggregated terms to include in created object (see below).}
  \item{aggregate}{Sum over individual components of aggregated terms if \code{TRUE}.}
  \item{log}{Use log10 of variances in plot if \code{TRUE}.}
  \item{percent}{Percentile between 0 and 100 for summaries.}
  \item{cex}{Character expansion for plot symbols. Default shrinks with
    number of MCMC iterations.}
  \item{\dots}{Arguments to pass along.}
}
\details{
Variance components are organized as "main" ("add" and "dom"),
"epistasis" ("aa", etc.), "fixcov" (for all fixed covariate terms),
"rancov" (random covariates), and "GxE" (genotype by environment,
including additive and dominance terms). Any subset of these may be chosen.
}
\value{
\code{qb.varcomp} creates a matrix with columns of samples for the
variance components. Each row represents an MCMC iteration. Values are
averaged over loci. 
}
\references{http://www.ssg.uab.edu/qtlbim}
\author{Brian S. Yandell}
\seealso{\code{\link{qb.mcmc}}}
\examples{
\dontshow{
qb.load(cross, qbExample)


}
temp <- qb.varcomp(qbExample)
summary(temp)
plot(temp)
}
\keyword{hplot}
