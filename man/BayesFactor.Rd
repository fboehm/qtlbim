\name{qb.BayesFactor}
\alias{qb.BayesFactor}
\alias{plot.qb.BayesFactor}
\alias{print.qb.BayesFactor}
\alias{summary.qb.BayesFactor}
\title{Bayesian model selection via Bayes factors.}
\description{
Model-averaged posteriors and Bayes factors computed for number and
pattern of QTL, chromosomes and pairs of chromosomes showing epistasis.
}
\usage{
qb.BayesFactor(qbObject, items = c("nqtl","pattern","chrom","pairs"),
  cutoff.pattern = 0.2, cutoff.pairs = 1, nmax = 15)
## S3 methods.
plot.qb.BayesFactor(x, \dots)
summary.qb.BayesFactor(object, sort = TRUE, digits = 3, \dots)
print.qb.BayesFactor(x, \dots)
}
\arguments{
  \item{qbObject}{An object of class \code{qb}.}
  \item{object}{Object of class \code{qb.BayesFactor}.}
  \item{x}{Object of class \code{qb.BayesFactor}.}
  \item{items}{Items to include in model selection assessment.}
  \item{cutoff.pattern}{Percent cutoff for pattern inclusion in model
    selection.} 
  \item{cutoff.pairs}{Percent cutoff for epistatic pair inclusion in
    model selection.}
  \item{nmax}{Maximum number of model terms included per item (for
    \code{items} "pattern" and "pairs" only).}
  \item{sort}{Sort by Bayes factor if \code{TRUE}.}
  \item{digits}{Number of significant digits for summary.}
  \item{\dots}{Additional arguments passed to generic plot, summary or print.}
}
\details{
  \code{qb.BayesFactor} creates model selection results for selected
  items. These are based on marginal posteriors and priors, averaged
  over all other model parameters. The posterior may be influenced by
  prior, while Bayes factors are empirically less sensitive for QTL
  model selection.  The Bayes factors are computed relative to the
  smallest term for each item, using the ratios of
  \code{posterior/prior}. Any pair of model terms can be compared as the
  ratio of their Bayes factors. The \code{items} evaluated are:
  \item{nqtl}{Number of QTLs.}
  \item{pattern}{Pattern of QTL across chromosomes. Identifiers are
    comma-separated chromosome numbers, with asterisk \code{n*} for
    multiple QTL per chromosome.}
  \item{chrom}{Chromosome.}
  \item{pairs}{Epistatic pairs of chromosomes.}
}
\value{
  List with \code{items}, each containing:
  \item{posterior}{Posterior frequency of MCMC samples.}
  \item{prior}{Prior frequency.}
  \item{bf}{Rank-ordered Bayes factors relative to smallest value.}
  \item{bfse}{Approximate standard error for bf computed using binomial
    variance of MCMC samples.}
}
\references{ \url{http://www.soph.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.qb}}, \code{\link{qb.mcmc}}}

\examples{
\dontshow{
qb.load(cross, qbExample)
  example(qb.mcmc)

}
temp <- qb.BayesFactor(qbExample)
summary(temp)
plot(temp)
}
\keyword{models}
