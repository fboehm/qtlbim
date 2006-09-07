\name{qb.meancomp}
\alias{qb.meancomp}
\alias{summary.qb.meancomp}
\alias{print.qb.meancomp}
\alias{plot.qb.meancomp}
\title{Examine grand mean and covariate MCMC samples.}
\description{
Examine grand mean and covariate Monte Carlo samples to glean estimates
of data center and importance of covariates.
}
\usage{
qb.meancomp(qbObject)
summary.qb.meancomp(object, percent = 5, ...)
print.qb.meancomp(x, ...)
plot.qb.meancomp(x, covar, percent = 5, cex, ...)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.meancomp}.}
  \item{x}{Object of class \code{qb.meancomp}.}
  \item{percent}{Percentile between 0 and 100 for summaries.}
  \item{covar}{Sequence of covariate identifiers for plot.}
  \item{cex}{Character expansion for plot symbols. Default shrinks with
    number of MCMC iterations.}
  \item{...}{Extra parameters passed along.}
}
\details{
Grand mean is adjusted to mean level of covariates. Diagonal of
scatterplot matrix includes density plot. Setting \code{covar = 0}
yields a density plot for the grand mean alone.
}
\value{
\code{qb.meancomp} is a matrix with columns for the grand mean and for
each fixed covariate. Summaries show mean and upper and lower percentiles.
}
\references{http://www.soph.uab.edu/qtlbim}
\author{Brian S. Yandell}
\seealso{\code{\link{qb.mcmc}}}
\examples{
\dontshow{
qb.load(cross, qbExample)


}
temp <- qb.meancomp(qbExample)
summary(temp)
plot(temp)
}
\keyword{hplot}
