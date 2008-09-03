\name{qb.close}
\alias{qb.close}
\alias{plot.qb.close}
\alias{summary.qb.close}
\title{Measures closeness of genetic architectures to target.}
\description{
Boxplots and summaries of how close MCMC samples of genetic
architectures are to target architecture.
}
\usage{
qb.close(qbObject, target = NULL, epistasis = TRUE, signed = FALSE,
  score.type = c("sq.atten","attenuation","variance","recombination",
  "distance"), \dots)
\method{plot}{qb.close}(x, category = c("pattern", "nqtl"), xlab,
  cutoff, sort.pattern = c("percent","score"), \dots)
\method{summary}{qb.close}(object, cutoff, digits = 0,
  show = "score", \dots )
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{x,object}{Object of class \code{qb.close}.}
  \item{target}{Target architure as data frame with columns \code{chrom}
    and \code{locus}. Extracted from \code{summary} if \code{target} is
    a \code{\link{qb.scanone}} object. If \code{target} = \code{NULL},
    \code{score.type} is set to \code{"variance"} to be compared with
    null model.}
  \item{epistasis}{Include epistasis in patterns if \code{TRUE}.}
  \item{signed}{Sign score. Most useful to examine single chromosome.}
  \item{score.type}{Type of score to use as distance.}
  \item{category}{Boxplots indexed by \code{nqtl} or \code{pattern}.}
  \item{xlab}{Label for X axis (default taken from \code{x} object).}
  \item{cutoff}{Percent cutoff for pattern inclusion in model
    selection. Default is 0.25 (0.5) if \code{epistasis} is \code{TRUE}
    (\code{FALSE}).}
  \item{sort.pattern}{If \code{type} = \code{"pattern"}, sort by
    \code{percent} posterior or by median of \code{score}.}
  \item{digits}{Number of digits displayed for locus.}
  \item{show}{Character string with name from \code{object} to show.}
  \item{\dots}{Parameters to methods.}
}
\details{
Closeness for each loci is measured as \code{1-2r}, with \code{r} the
recombination rate. Thus unlinked loci have measure 0. Loci between
a MCMC sample architecture and the target architecture on the same
chromosome are matched by closest distance in cM between subsets of the
the same length (if target has 2 QTL on chr 3 and sample has 3, consider
all pairs from sample to find closest pair in 2-D). Measure per sample
is sum across all loci. A quick way to generate a target is to use
\code{\link{qb.BestPattern}}. 

The \code{score.type} is \code{"recombination"} = r, the recombination
rate; \code{"attenuation"} = 1-2r; \code{"sq.atten"} = squared
attenuation, \code{"distance"} in cM, or genetic \code{"variance"}.
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{qb.BestPattern}}, \code{\link{qb.hpdone}}}

\examples{
data(qbExample)

## Here target is known for simulated data.
target <- data.frame(chrom = c(1,1,2,3),
  locus = c(15,45,12,15))
temp <- qb.close(qbExample, target)
summary(temp)
plot(temp)
plot(temp, category = "nqtl")

## Or use target from qb.BestPattern
target <- summary(qb.BestPattern(qbExample))$best
qb.close(qbExample, target)
}
\keyword{models}
