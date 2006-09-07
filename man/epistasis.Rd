\name{plot.qb.epistasis}
\alias{qb.epistasis}
\alias{plot.qb.epistasis}
\alias{print.qb.epistasis}
\alias{summary.qb.epistasis}
\title{Density Plots for Models Showing Epistasis}
\description{
  Produces density plots of models showing epsistasis. The vertical axis
shows degree of epstatic effect, horizontal axis shows identities of
pairs.  Paralell plots are produced for each of the entries in the
\code{effects} parameter.}
\usage{
qb.epistasis(qbObject, effects = c("aa", "ad", "da", "dd"),
  cutoff = 1, maxpair = 5, pairs, \dots)
## S3 methods.
plot.qb.epistasis(x, effects = c("aa", "ad", "da", "dd"), cex = 0.5, \dots )
print.qb.epistasis(x, \dots )
summary.qb.epistasis(object, \dots )
}
\arguments{
  \item{qbObject}{An object of class \code{qb}.}
  \item{object}{Object of class \code{qb.epistasis}.}
  \item{x}{Object of class \code{qb.epistasis}.}
  \item{cutoff}{The \code{cutoff} parameter for number of epistatic pairs.}
  \item{maxpair}{Maximum number of epistatic pairs shown.}
  \item{pairs}{A list of chromosome pairs to examing for epistatic
    pairs.  The default uses chromosome names separated by a dot.}
  \item{effects}{Paralell plots are produced for each of the entries in the
\code{effects}.}
  \item{cex}{Horizontal expansion factor for characters in the plot.
    See \code{\link[graphics]{par}}.}
  \item{\dots}{Arguments passed to generic \code{\link[graphics]{plot}}.}
}

\value{
  Returns a table of counts of epistatic pairs with counts above the
  cutoff value.
}
\references{\url{http://www.soph.uab.edu/qtlbim}}
\author{Brian S. Yandell, yandell@stat.wisc.edu }

\examples{
\dontshow{
qb.load(cross, qbExample)


}
temp <- qb.epistasis(qbExample)
summary(temp)
plot(temp)
}

\keyword{hplot}
\keyword{models}
