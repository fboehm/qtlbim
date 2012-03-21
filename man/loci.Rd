\name{plot.qb.loci}
\alias{qb.loci}
\alias{plot.qb.loci}
\alias{print.qb.loci}
\alias{summary.qb.loci}
\title{Jittered plot of Bayesian QTL loci samples by chromosome}
\description{
Each point is one locus from the Bayesian QTL estimates, plotted
vertically by chromosome, jittered to give a sense of density.
Separate colored vertical bands by loci element.
}
\usage{
qb.loci(qbObject, loci = c("main", "epistasis", "GxE"), covar, \dots)
\method{plot}{qb.loci}(x, loci, labels = FALSE, amount = 0.35, cex, col, \dots)
\method{print}{qb.loci}(x, \dots)
\method{summary}{qb.loci}(object, digit = 1, \dots)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{loci}{Character string identifying MCMC sample elements; may include
    \code{"main"}, \code{"epistasis"}, \code{"GxE"} and \code{"all"}.}
  \item{covar}{Fixed covariate(s) for \code{"GxE"} loci; default is all
    fixed covariates involved in GxE interactions.}
  \item{x}{Object of class \code{qb.loci}.}
  \item{object}{Object of class \code{qb.loci}.}
  \item{labels}{Include marker labels if TRUE.}
  \item{amount}{Amount of \code{\link[base]{jitter}} (between 0 and .45)}
  \item{cex}{Character expansion (may be invisible if too
    small--default based on number of MCMC samples).}
  \item{col}{Character string with colors named by \code{loci}; also
    includes color for \code{marker} lines.}
  \item{digit}{Number of digits for roundoff of loci quantiles.}
  \item{\dots}{Graphical parameters can be given as arguments to
    \code{\link[graphics]{plot}}. Not used in \code{qb.loci}.} 
}
\details{
  Focuses attention on chromosome lengths and concentration of QTL loci
  estimates. Horizontal lines at markers. Separate bands by \code{loci}
  for each chromosome. Adjust \code{amount} and \code{cex} to modify
  look of plot. Most useful when looking at multiple chromosomes.
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link[base]{jitter}},\code{\link{subset.qb}}}

\examples{
data(qbExample)

temp <- qb.loci(qbExample)
plot(temp)
summary(temp)
temp <- qb.loci(qbExample, "all")
plot(temp)
summary(temp)
}
\keyword{models}
