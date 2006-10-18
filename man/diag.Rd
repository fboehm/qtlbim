\name{plot.qb.diag}
\alias{qb.diag}
\alias{plot.qb.diag}
\alias{print.qb.diag}
\alias{summary.qb.diag}
\title{Marginal and model-conditional summaries of Bayesian interval
  mapping diagnostics}
\description{
A density histogram is drawn for model-averaged summary diagnostics such
as LOD, variance, or heritability.
}
\usage{
qb.diag(qbObject, items= c("mean","envvar","var","herit"), \dots)
\method{plot}{qb.diag}(x, \dots )
\method{print}{qb.diag}(x, \dots )
\method{summary}{qb.diag}(object, digits = 5, \dots )
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{object}{Object of class \code{qb.diag}.}
  \item{x}{Object of class \code{qb.diag}.}
  \item{items}{Diagnostics to be summarized; must be name of a column in
    \code{element}.}
  \item{digits}{Number of significant digits.}
  \item{\dots}{Parameters to methods. Not used for \code{qb.diag}.}
}
\details{
  Model-averaged density is smooth kernel estimate similar to ordinary
  histogram. A \code{\link[base]{boxplot}} (without outliers) is
  overlaid for comparison   with conditional boxplots. Conditional
  boxplots by number of QTL may   show   indication of model bias for
  small number of QTL. This and  \code{\link{qb.BayesFactor}} can help
  suggest the minimal model. Diagnostic   items that make sense to plot
  are \code{"LOD"}, \code{"envvar"}   (environmental variance),
  \code{"herit"} (heritability), \code{"mean"} (grand mean),
  \code{"addvar"} (variance of \code{add}), \code{"domvar"} (variance of
  \code{add}). Marginal and conditional medians are printed.
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.qb}}, \code{\link[stats]{density}},
  \code{\link[graphics]{boxplot}}, \code{\link{qb.BayesFactor}}}

\examples{
\dontshow{
qb.load(cross, qbExample)
}
temp <- qb.diag(qbExample)
summary(temp)
plot(temp)
}
\keyword{models}
