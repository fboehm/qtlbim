\name{Bnapus}
\alias{Bnapus}
\title{Cross structure for complete Brassica napus data}
\description{
Contains genotypes and phenotypes for Brassica napus study, including 0-
4- and 8-week vernalization, survival, and 19 chromosomes.
}
\usage{data(Bnapus)}
\format{
See \code{read.cross} in \code{library(qtl)} for format.
}
\seealso{\code{\link[qtl]{read.cross}},\code{\link{plot.qb}}}
\source{
Thomas C. Osborn (\url{mailto:tcosborn@facstaff.wisc.edu}), Department of Agronomy,
UW-Madison (\url{http://agronomy.wisc.edu}).
}
\details{
Traits included are percent winter survival for 1992-3, 1993-4,
1994-5, 1997-8, and 1999-2000 (surv92, surv93, surv94, surv97, surv99,
respectively), and days to flowering after no vernalization (flower0),
4 weeks vernalization (flower4) or 8 weeks vernalization
(flower8). Percents are of plants alive in the Fall, taken from the
middle of rows (totals unavailable). Days to flowering after
transplant are averages over four replicates from a RCB design (values
by block unavailable). First column has the trait name. The remaining
columns identify individual DH line (302-455). 

Marker genotype data for Major x Stellar double haploid (DH)
population. Double haploids have the same relation of recombination to
distance as backcrosses as there is just one meiosis tracked (in
F1). However, DH are homozygous at every locus (usually mapped as RI0
lines). Marker genotypes are coded as M = Major, S = Stellar, - =
missing. Data columns are
\item chrom = <em>B. napus</em> chromosome (N1, N2, etc.)
\item order = along chromosome
\item cM d= istance from proximal end
\item marker = marker name:
 E = AFLP;
 *ec, *tg, *wg = RFLP;
 *xxx = other markers from Arabidopsis:
Lem, eru1, eru2, fad3, isoDia, isoIdh, isoPgi, isoLap, pr2, slg6, Aca1, cor15
\item
Remaining columns are for individual DH lines (identifier 302-455).
}
\references{
  \item \url{http://www.stat.wisc.edu/~yandell/qtl/data/osborn/Bnapus}
  \item Ferreira ME, Satagopan J, Yandell BS, Williams PH, Osborn TC
  (1995) Mapping loci controlling vernalization requirement and
  flowering time in Brassica napus. Theor Appl Genet
  90: 727-732. [original source and analysis]
  \item JM Satagopan, BS Yandell, MA Newton and TC Osborn (1996)
  Markov chain Monte Carlo approach to detect polygene loci
  for complex traits. Genetics 144: 805-816.
  \url{http://www.genetics.org/cgi/content/abstract/144/2/805}
  [first MCMC for experimental crosses; analysis of B. napus N2=LG9; see
  \code{\link{vern}} data]
  \item Kole C, Thorman CE, Karlsson BH, Palta JP, Gaffney P, Yandell
  BS, Osborn TC (2001)
  Comparative mapping of loci controlling winter survival and
  related traits in oilseed Brassica rapa and B. napus.
  Molecular Breeding 1: 329-339. [refined map and reanalysis]
}
\examples{
data(Bnapus)
summary(Bnapus)
}
\keyword{datasets}
