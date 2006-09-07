\name{qb.model}
\alias{qb.model}
\title{Set up interacting QTL model for qb.mcmc}
\description{
 This function sets up a genome-wide interacting QTL model by specifying global
 constraints on models and priors on unknowns.
}
\usage{
qb.model(cross, epistasis = TRUE, main.nqtl = 3,
  mean.nqtl = main.nqtl + 3, max.nqtl = NULL, interval = NULL,
  chr.nqtl = NULL,intcov = c(0), depen = FALSE,
  prop = c(0.5, 0.1, 0.05), \dots)	     
  
}
\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{epistasis}{indicates if epistasis is included in the model: TRUE or FALSE}
  \item{main.nqtl}{prior expected number of main effect QTLs. }
  \item{mean.nqtl}{prior expected number for all QTLs on all chromosomes including 
    QTLs with main effects, epistatic effects and gene-environment interactions.}
  \item{max.nqtl}{maximum number of QTLs allowed in the model. Default
    is \eqn{l+3\sqrt{l}} where \eqn{l} is \code{main.qtl} for
    non-epistatic model and \code{mean.qtl} for epistatic model. }
  \item{interval}{minimum distance between any two flanking QTLs for all chromosomes.
    Default is the average distance between markers in each chromosome.}
  \item{chr.nqtl}{list of the maximum number of QTLs allowed to be detected 
    on each chromosome. Default is the length of the chromosome
    divided by \code{interval}.}
  \item{intcov}{list of the fixed covariates with which gene-environment interaction 
    will be considered.}
  \item{depen}{=TRUE will use dependent prior for indicator variables of epistatic effects.}
  \item{prop}{prior inclusion probabilities for epistatic effects in three different scenarios: 
    when none,one and both the main effects of the two interacting QTL are 
    included in the model. Note that the sum of the probabilities need not
    be equal to 1 and \code{prop} should be specified only when \code{depen=TRUE}.
  }
  \item{\dots}{Not used.}
}
\value{
  \item{qtl_envi}{ Indicates if there is an interaction between the QTLs and
                   environmental variables: TRUE or FALSE.}
}
\details{
This function defines the model for Bayesian QTL mapping using \code{\link{qb.mcmc}}.
This model considers two-way interaction as the highest level of both gene-gene and
gene-environment interactions.
}
\references{
  \url{http://www.soph.uab.edu/qtlbim}
  }


\note{ This function returns a list and hence should have a differenct name from
that of the \code{cross} object.
}


\author{Dr. Nengjun Yi, et al.,  nyi@ms.soph.uab.edu}

\seealso{\code{\link{qb.data}}, 
  \code{\link{qb.genoprob}},   \code{\link{qb.mcmc}}}


\examples{
\dontshow{
if(!exists("cross"))
  example(qb.genoprob)

}
qbModel <- qb.model(cross, chr.nqtl = rep(3,nchr(cross)),
  prop = c(0.5, 1, 0.5), intcov = 1)  
}
\keyword{models}
