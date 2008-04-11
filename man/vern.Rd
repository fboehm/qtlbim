\name{vern}
\alias{vern}
\title{Eight week vernalization data for Brassica napus}
\description{
Contains genotypes and phenotypes for 8-week vernalization study used in
Satagopan et al. (1996).
}
\usage{
data(vern)
}
\format{
See \code{\link[qtl]{read.cross}} for format of
\code{vern}.
}
\seealso{\code{\link[qtl]{read.cross}}, \code{\link{plot.qb}},
\code{\link{qb.mcmc}}}
\source{
Thomas C. Osborn (\url{mailto:tcosborn@facstaff.wisc.edu}), Department of Agronomy,
UW-Madison. 
}
\references{
  \item \url{http://www.stat.wisc.edu/~yandell/qtl/data/osborn/Bnapus}
  \item Ferreira ME, Satagopan J, Yandell BS, Williams PH, Osborn TC
  (1995) Mapping loci controlling vernalization requirement and
  flowering time in Brassica napus. Theor Appl Genet
  90: 727-732. [original source and analysis]
  \item Kole C, Thorman CE, Karlsson BH, Palta JP, Gaffney P, Yandell
  BS, Osborn TC (2001)
  Comparative mapping of loci controlling winter survival and
  related traits in oilseed Brassica rapa and B. napus.
  Molecular Breeding 1: 329-339. [refined map and reanalysis]
  \item JM Satagopan, BS Yandell, MA Newton and TC Osborn (1996)
  Markov chain Monte Carlo approach to detect polygene loci
  for complex traits. Genetics 144: 805-816.
  \url{http://www.genetics.org/cgi/content/abstract/144/2/805}
  [first MCMC for experimental crosses; analysis of B. napus N2=LG9; see
  \code{\link{vern}} data]
}
\examples{
data(vern)
summary(vern)
}
\keyword{datasets}
