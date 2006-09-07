\name{fisch}
\alias{fisch}
\non_function{}
\title{Eight QTL Stephens and Fisch simulated data}
\description{
Contains genotypes and phenotypes for data simulated using model in
Stephens and Fisch (1998) but with 90 percent heritability.
}
\usage{
data(fisch)
}
\format{
  \code{fisch} is f2 (see \code{\link[qtl]{read.cross}} for format).
}
\seealso{\code{\link[qtl]{read.cross}}, \code{\link{plot.qb}}, \code{\link{qb.mcmc}}}
\author{ Brian S. Yandell, \url{mailto:yandell@stat.wisc.edu}}
\source{
Patrick J. Gaffney (\url{mailto:paga@lubrizol.com}), Lubrizol Corp.
}
\references{
\url{http://www.soph.uab.edu/qtlbim}
}
\examples{
data(fisch)
summary(fisch)
}
\keyword{datasets}
