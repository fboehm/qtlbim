\name{qb.arch}
\alias{qb.arch}
\alias{summary.qb.arch}
\alias{print.qb.arch}
\alias{qb.arch.default}
\alias{qb.arch.step.fitqtl}
\title{Build genetic architecture with chromosomes, positions and
  epistatic pairs.}
\description{
These routines work in conjunction with \code{\link{qb.hpdone}},
\code{\link{qb.scantwo}}, \code{\link{qb.sliceone}} and
\code{\link{step.fitqtl}} to infer the number, pattern and position of
QTL from MCMC samples.
}
\usage{
qb.arch(object, \dots)
qb.arch.default(object, chr, pos, tolerance = 10, \dots)
qb.arch.step.fitqtl(object, main, epistasis, \dots)
\method{summary}{qb.arch}(object, \dots)
\method{print}{qb.arch}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object for appropriate method: summary of object of
    class \code{\link{qb.scantwo}} for default; object of class
    \code{\link{step.fitqtl}}.}
  \item{x}{Object of class \code{qb.arch}.}
  \item{chr}{Vector of chromosome numbers.}
  \item{pos}{Vector of positions on chromosomes (much be same length as \code{chr}).}
  \item{tolerance}{Minimum distance for two QTL to be considered
    distinct.}
  \item{main}{Vector of chromosome identifiers with only main effects.}
  \item{epistasis}{Data frame with a 2-element vector of chromosome
    identifiers for each epistatic pair.}
  \item{\dots}{Not used here.}
}
\details{
Extract architecture in terms of chromosomes and positions of main QTL
and identifiers of epistatic pairs of QTL. The \code{\link{step.fitqtl}}
approach is used to compare an automatic fit to a user-defined set of
main chromosomes and epistatic pairs.
}
\value{
  \item{qtl}{Data frame with main QTL as \code{chr} and \code{pos}.}
  \item{by.num}{Data frame with epistatic pairs indexed by chromosome
    number, labeled \code{qtla} and \code{qtlb}.}
  \item{by.chr}{List with elements \code{chr} and \code{pos} showing
    epistatic pairs. These elements are data frames with chromosomes and
    positions for each epistatic pair: rows are QTL number, columns are
    \code{qtla} and \code{qtlb}.} 
  \item{by.set}{List of connected sets of epistatic chromosomes.}
}
\references{http://www.ssg.uab.edu/qtlbim}
\author{Brian S. Yandell}
\seealso{\code{\link{step.fitqtl}}, \code{\link{qb.sweave}}}

\examples{
\dontshow{
qb.load(cross, qbExample)
}
## Run qb.scantwo and get summary to use in qb.arch
temp <- summary(qb.scantwo(qbExample, type = "2logBF"),
  threshold = c(upper = 10))
## qb.arch default use. 
cross.arch <- qb.arch(temp, chr = c(1,1,2,3), pos = c(15,45,12,15))
cross.arch
}
\keyword{utilities}
