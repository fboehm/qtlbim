\name{plot.qb.mainmodes}
\alias{qb.epimodes}
\alias{qb.mainmodes}
\alias{qb.split.chr}
\alias{plot.qb.mainmodes}
\alias{summary.qb.mainmodes}
\title{Summaries of QTL modes by chromosome.}
\description{
Determine number of QTL per chromosome and estimate peaks and valleys.
}
\usage{
qb.split.chr(qbObject, split, \dots)
qb.epimodes(qbObject, cutoff = 1, nqtl, n.iter, pairloci, \dots)
qb.mainmodes(qbObject, cutoff = 25, nqtl, n.iter, mainloci, \dots)
\method{summary}{qb.mainmodes}(object, digits = 4, \dots)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{split}{List of split locations; names of list must correspond to
    chromosome names in \code{qbObject}. Default set by call to
    \code{qb.mainmodes} with addition arguments in \code{\dots}.}
  \item{object}{Object of class \code{qb.mainmodes}.}
  \item{cutoff}{Cutoff for neglible number of QTL.}
  \item{nqtl}{Vector of number of QTL per chromosome.}
  \item{n.iter}{Number of iterations (for internal use only).}
  \item{mainloci}{Object containing \code{mainloci} data (for internal
    use only).}
  \item{pairloci}{Object containing \code{pairloci} data (for internal
    use only).}
  \item{digits}{Number of significant digits.}
  \item{\dots}{Parameters to \code{\link{qb.mainmodes}} or to methods.}
}
\details{
Cut off histogram of number of QTL per chromosome such that cumulative
percent above number is less than \code{cutoff}. Once \code{nqtl} is
determined or provided in call, divide MCMC samples using linear
discriminant analysis (\code{\link[MASS]{lda}}) and find peak locations per
class. Use these peak locations to find locations of valleys between
peaks. These valleys are used to divide MCMC samples into separate QTL
for analysis. Currently this is used by \code{\link{summary.qb.scanone}}
and \code{\link{qb.multloci}}.

\code{qb.split.chr} sets up a \code{split} for chromosomes
for at valleys between inferred multiple QTL as an
attribute of a returned \code{qbObject}. This is done by default when
\code{qbObject} is created, and is stored as an attribute available as
\code{qb.get(qbObject, "split.chr")}.
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{qb.multloci}}, \code{\link{summary.qb.scanone}}}

\examples{
data(qbExample)

temp <- qb.mainmodes(qbExample)
summary(temp)
}
\keyword{models}
