\name{plot.qb.loci}
\alias{qb.loci}
\alias{plot.qb.loci}
\title{Jittered plot of Bayesian QTL loci samples by chromosome}
\description{
Each point is one locus from the Bayesian QTL estimates, plotted
vertically by chromosome, jittered to give a sense of density.
}
\usage{
qb.loci(qbObject, \dots)
plot.qb.loci(x, labels=FALSE, amount=0.35, cex, \dots )
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{x}{Object of class \code{qb.loci}.}
  \item{labels}{Include marker labels if TRUE.}
  \item{amount}{Amount of \code{\link[graphics]{jitter}} (between 0 and .45)}
  \item{cex}{Character expansion (may be invisible if too
    small--default based on number of MCMC samples).}
  \item{\dots}{Graphical parameters can be given as arguments to
    \code{\link[base]{plot}}. Not used in \code{qb.loci}.} 
}
\details{
  Focuses attention on chromosome lengths and concentration of QTL loci
  estimates. Horizontal blue lines at markers if \code{cross}
  included. Adjust \code{amount} and \code{cex} to modify look of 
  plot. Most useful when looking at multiple chromosomes.
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link[base]{jitter}},\code{\link{subset.qb}}}

\examples{
\dontshow{
qb.load(cross, qbExample)
  example(qb.mcmc)

}
temp <- qb.loci(qbExample)
plot(temp)
}
\keyword{models}
