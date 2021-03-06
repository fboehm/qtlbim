\name{qb.genoprob}
\alias{qb.genoprob}
\title{Grid point and genotype probability computation method}
\description{
  This function is used to compute putative QTL positions and genotypic probabilities
  at these positions. The genotypic probabilities for missing marker genotypes are also
  computed.
}
\usage{
qb.genoprob(cross, map.function, step = 2,
  tolerance = 1e-6, stepwidth = "variable", \dots)
} 
\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.}
 
  \item{map.function}{Indicates what map function to use when converting
    genetic distances to recombination fractions. See \code{\link[qtl]{calc.genoprob}}.}
  \item{step}{Distance (in cM) between positions at which putative QTL
    positions and their genotype probabilities are calculated. However,
    specifying \code{step = 0} would assume marker positions as putative
    QTL locations and genotypic probabilities would be calculated only
    for markers with missing genotype.}
  \item{tolerance}{Minimum separation of markers enforced by
    \code{\link[qtl]{jittermap}}.}
  \item{stepwidth}{Type of stepwidth for
    \code{\link[qtl]{calc.genoprob}}; \code{"variable"} works best with
    MCMC sampling.}
  \item{\dots}{Extra arguments to pass to \code{\link[qtl]{calc.genoprob}}.}
}
\value{
  \code{qb.genoprob} first ensures marker separation is at least
  \code{tolerance}, and then computes genotype probabilities at
  pseudomarkers spaced approximately \code{step} units apart using
  \code{\link[qtl]{calc.genoprob}}. See \code{\link[qtl]{calc.genoprob}}
  for value.
}
\references{ \url{http://www.qtlbim.org} }
\author{Dr. Nengjun Yi, et al.,  nyi@ms.ssg.uab.edu}

\seealso{\code{\link[qtl]{jittermap}}, \code{\link[qtl]{calc.genoprob}}.}

\examples{
\dontshow{
if(!exists("cross"))
  example(qb.sim.cross)

}
## calculate grids and genotypic probabilites
cross <- qb.genoprob(cross,map.func="haldane",step=2)
}
\keyword{models}

