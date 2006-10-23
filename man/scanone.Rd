\name{qb.scanone}
\alias{qb.scanone}
\title{Genome Scan for Main Loci Involved in Phenotypic Trait}
\description{
This method extracts iteration diagnostics and mainloci from the
\code{qb} object and returns a data frame (of class \code{qb.scanone})
containing information on environmental variance, explained variance
components, non-epistatic variance components.
}
\usage{
qb.scanone(qbObject, epistasis = TRUE, scan, type, covar, adjust.covar,
  chr, sum.scan = "yes", min.iter = 1,
  aggregate = TRUE, half = FALSE, verbose = FALSE)
}
\arguments{
  \item{qbObject}{An object of class \code{qb}.}
  \item{epistasis}{If \code{TRUE} then information about
    epistasis is included.}
  \item{scan}{Vector of diagnostics to scan (see below).}
  \item{type}{Type of scan; default is "heritability" (see below).}
  \item{covar}{Covariate(s) to include; default is \code{seq(nfixcov)}
    where \code{nfixcov} is taken from \code{\link{qb.data}}. Set to 0
    to exclude any covariates.}
  \item{adjust.covar}{Adjustments to covariates. Default is
    \code{NA}, which adjusts by covariate mean values. Values are
    assumed to be in order of fixed covariates.}
  \item{chr}{Chromosomes to subset on if not \code{NULL}.}
  \item{sum.scan}{Sum over \code{scan} diagnostics if "yes" or "only";
    only report \code{sum} if "only".}
  \item{min.iter}{Include only samples at loci if minimum number of
    iterations is at least \code{min.iter}; default is to include all
    (\code{min.iter} = 1).}
  \item{aggregate}{Aggregate effects into main, epis, gbye if \code{TRUE}.}
  \item{half}{Cut epistatic effects in half if \code{TRUE}.}
  \item{verbose}{ Give verbose feedback if \code{TRUE}.}
}

\value{
Returns an object of class \code{qb.scanone} (a data frame) containing
effects selected according to \code{type} and \code{scan}.
}
\details{
The \code{type} specifies what type of scan is performed. Scan produces
marginal estimates of diagnostics at each potential loci across the
genome. That is, values are adjusted for other possible QTL simply by
taking the marginal average over MCMC samples. Choices of \code{type} are
"heritability", "LPD", "LR", "deviance", "detection", "variance",
"estimate", "cellmean", "count", "log10", "posterior", "logposterior"
(i.e. log10(posterior)), "BF", "2logBF" (i.e. 2*ln(BF)), and "nqtl"
(number of linked QTL). Default is "LPD".

Type "heritability" is actually R-squared at this point, not the theoretical
heritability.
Types "LPD", "LR" and "deviance" are all proportional to each other in
the usual sense; "LPD" is computed to agree with \code{lod} from
\code{\link[qtl]{scanone}} if models were restricted to one QTL and missing
genotypes are imputed. Detection is the marginal posterior probability of
detectio of a QTL at a locus.
Types "variance" and "estimate" yield, respectively, the marginal variance
components and the marginal parameter estimates at each loci. Type
"cellmean" gives marginal estimates for A, H, B genotypes (these are
single character codes for AA, AB, BB, respectively). The remaining
count types provide diagnostics. Types "count" and "log10" report on
number of MCMC samples in raw or logged scale. Type 
"posterior" ("logposterior") yields the marginal (log) posterior
probability. Type "BF" ("2logBF") gives the marginal
Bayes factor per loci; both are proportional to "count". Type "nqtl"
gives the average number of linked loci, which can be useful in sorting
out multiple linked loci.

The \code{scan} specifies the model effects to include for all types
except the counts. Aggregated effects (default except for type
"cellmean") are "main", "epistasis" and "GxE" (genotype by
environment). Individual model effects can be requested 
as "add", "dom", "aa", "ad", "da", "dd". In addition, GxE terms, if
present are included automatically if \code{covar} is not 0. For type
"estimate", main effects for "add" and "dom" are adjusted for any
covariate GxE effects. The \code{sum.scan} is used for all types but the
counts to get a summary across \code{scan} effects.
}

\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link{summary.qb.scanone}}, \code{\link{plot.qb.scanone}}}
\examples{
\dontshow{
qb.load(cross, qbExample)

}
temp <- qb.scanone(qbExample)
summary(temp)
plot(temp)
}
\keyword{models}
\keyword{hplot}
