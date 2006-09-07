\name{qb.covar}
\alias{qb.covar}
\alias{summary.qb.covar}
\alias{print.qb.covar}
\alias{plot.qb.covar}
\title{Examine GxE effect of covariates on main genetic effects.}
\description{
Compare main effects with GxE effects to address correlation of estimates.
}
\usage{
qb.covar(qbObject, element = "add", covar = 1, chr, ...)
summary.qb.covar(object, percent = 5, digits = 3, ...)
print.qb.covar(x, ...)
plot.qb.covar(x, percent = 5, cex, include.zero = TRUE, ...)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.covar}.}
  \item{x}{Object of class \code{qb.covar}.}
  \item{element}{Main effect to examine ("add" or "dom").}
  \item{covar}{Index to covariates used in MCMC samples.}
  \item{chr}{Subset of chromosomes as integer vector.}
  \item{percent}{Percentile (0 to 100) for summaries.}
  \item{digits}{Number of significant digits to print.}
  \item{cex}{Character expansion for plots (default decreases with MCMC
    sample size).}
  \item{include.zero}{Include zero values in plot when \code{TRUE}.}
  \item{\dots}{Arguments passed through to inherited routines.}
}
\details{
The diagonal dark green line of points on plots by chromosome indicate
adjustment for covariates that have not been centered. Main effects are
generally less correlated with GxE when covariates are first centered to
have mean zero.
}
\value{
Objects of class \code{qb.covar} have three columns: main effect, GxE
effect and chromosome. Summary objects have eight columns, three for main
effect and GxE (mean, lower and upper percentile), followed by
correlation and p-value. Summaries are done by chromosome.
}
\references{http://www.soph.uab.edu/qtlbim}
\author{Brian S. Yandell}
\seealso{\code{\link{qb.mcmc}}}
\examples{
\dontshow{
qb.load(cross, qbExample)


}
temp <- qb.covar(qbExample)
summary(temp)
plot(temp)
}
\keyword{hplot}
