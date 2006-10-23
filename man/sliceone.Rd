\name{qb.sliceone}
\alias{qb.sliceone}
\alias{print.qb.sliceone}
\alias{summary.qb.sliceone}
\alias{plot.qb.sliceone}
\title{Genome Slice to detect QTL for Phenotypic Trait}
\description{
This method extracts iteration diagnostics and mainloci from the
\code{qb} object and returns a data frame (of class \code{qb.sliceone}).
Generic summary and plot can be used for display.
}
\usage{
qb.sliceone(qbObject, slice, epistasis = TRUE, scan, type, covar,
  adjust.covar, chr, sum.scan = "yes", min.iter = 1,
  aggregate = TRUE, verbose = FALSE)
\method{summary}{qb.sliceone}(object, chr, \dots)
\method{print}{qb.sliceone}(x, \dots)
\method{plot}{qb.sliceone}(x, \dots, scan, auto.par = TRUE)
}
\arguments{
  \item{qbObject}{An object of class \code{qb}.}
  \item{object}{Object of class \code{qb.sliceone}.}
  \item{x}{Object of class \code{qb.sliceone}.}
  \item{slice}{Chromosomes to slice upon.}
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
  \item{verbose}{ Give verbose feedback if \code{TRUE}.}
  \item{auto.par}{Automatic setting of plot parameters for multiple
    plots if \code{TRUE}.}
  \item{\dots}{Arguments to be passed along.}
}

\value{
\code{qb.sliceone} returns an object of class \code{qb.sliceone} (a data frame) containing
effects selected according to \code{type} and \code{scan}.
}
\details{
All arguments except \code{slice} agree with \code{\link{qb.scanone}}. The
slice specifies a chromosome upon which to slice, yielding a 1-D scan of
what might be seen on a 2-D scan using \code{\link{qb.scantwo}}. One
advantage of \code{qb.sliceone} is that you can get 2-QTL cell means for the
slice chromosome and the scanned chromosomes.

The summary invokes \code{\link{summary.qb.scanone}} to summarize slice by
chromosome. The plot will by default give separate plots for each slice
genotype and use \code{\link{plot.qb.scanone}} to scan the chromosomes. If
\code{scan} is specified for \code{plot.qb.sliceone}, then those elements
will be plotted. For instance, \code{plot(x,scan="slice")} will plot the
running average locus on the slice chromosome with respect to the other
chromosomes.
}

\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link{summary.qb.scanone}}, \code{\link{plot.qb.scanone}}}
\examples{
\dontshow{
qb.load(cross, qbExample)


}
## Get profile of heritability.
temp <- qb.sliceone(qbExample, slice = 1, chr = 2:3)
summary(temp)
plot(temp)

## Get profile of cell means.
temp <- qb.sliceone(qbExample, slice = 1, chr = 2:3, type = "cellmean")
summary(temp)
plot(temp)
}
\keyword{models}
\keyword{hplot}
