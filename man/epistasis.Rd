\name{plot.qb.epistasis}
\alias{qb.epistasis}
\alias{qb.intcov}
\alias{plot.qb.epistasis}
\alias{print.qb.epistasis}
\alias{summary.qb.epistasis}
\title{Density Plots for Models Showing Epistasis and GxE Interactions.}
\description{
Produces density plots of models showing epsistasis (\code{qb.epistasis})
or GxE interactions (\code{qb.intcov}). The vertical axis
shows magnitude of effect, horizontal axis shows chromosomes in
epistatic pairs or covariate by chromosome.  Paralell plots are produced
for each of the entries in the \code{effects} parameter.}
\usage{
qb.epistasis(qbObject, effects = c("aa", "ad", "da", "dd"),
  cutoff = 1, maxpair = 5, pairs, \dots)
qb.intcov(qbObject, covar, effects = c("add","dom"),
  cutoff = 1, nmax = 5, cov.chr, \dots)
\method{plot}{qb.epistasis}(x, effects, cex = 0.5, main, \dots )
\method{print}{qb.epistasis}(x, \dots )
\method{summary}{qb.epistasis}(object, \dots )
}
\arguments{
  \item{qbObject}{An object of class \code{qb}.}
  \item{object}{Object of class \code{qb.epistasis}.}
  \item{x}{Object of class \code{qb.epistasis}.}
  \item{cutoff}{The \code{cutoff} parameter for number of epistatic pairs.}
  \item{maxpair}{Maximum number of epistatic pairs shown.}
  \item{pairs}{A character vector of chromosome pairs to examing for epistatic
    pairs.  Chromosome names are separated by a dot.}
  \item{covar}{Covariate(s) to include; default is \code{seq(nfixcov)}
    where \code{nfixcov} is taken from \code{\link{qb.data}}.}
  \item{nmax}{Maximum number of covariate chromosomes shown.}
  \item{cov.chr}{A character vector of covariate by chromosome pairs to
    examing for GxE effects.  Covariate names and chromosome names are
    separated by a dot.}
  \item{effects}{Character string of model effects.}
  \item{cex}{Horizontal expansion factor for characters in the plot.
    See \code{\link[graphics]{par}}.}
  \item{main}{Main titles for plots; default is \code{effects}.}
  \item{\dots}{Arguments passed to generic \code{\link[graphics]{plot}}.}
}

\value{
  Returns a table of counts of epistatic pairs with counts above the
  cutoff value.
}
\references{\url{http://www.qtlbim.org}}
\author{Brian S. Yandell, yandell@stat.wisc.edu }

\examples{
\dontshow{
qb.load(cross, qbExample)

}
temp <- qb.epistasis(qbExample)
summary(temp)
plot(temp)
temp <- qb.intcov(qbExample)
summary(temp)
plot(temp)
}

\keyword{hplot}
\keyword{models}
