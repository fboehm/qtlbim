\name{qb.data}
\alias{qb.data}
\title{Prepares data for qb.mcmc}
\description{
This function selects trait(s) and covariates from a \code{cross} object to build
a model (\code{\link{qb.model}}) for MCMC (\code{\link{qb.mcmc}}).  
}
\usage{
qb.data( cross, pheno.col = 1, trait = c("normal","binary","ordinal"), 
  fixcov = c(0), rancov = c(0), boxcox = FALSE, standardize = FALSE, \dots) 
}
\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{the column number for the phenotype used by \code{model}. 
                    Currently, only one phenotype can be analyzed at a time.}
  \item{trait}{Type of the quatitative trait or dependent variable: 
                "normal" or "binary" or "ordinal".}
  \item{fixcov}{list of fixed covariates. The column number(s) in \code{cross$pheno}
                which is(are) considered as fixed covariates.}
  \item{rancov}{list of random covariates.The column number(s) in \code{cross$pheno}
                which is(are) considered as random covariates.}
  \item{boxcox}{ Indicates whether to use a Boxcox transformation for the 
                  dependent variable or not: TRUE or FALSE. Note: trait has to be "normal"
                  and all phenotypic values have to be positive for using this option.}
  \item{standardize}{Indicates whether to standardize the dependent variable or not: 
    TRUE or FALSE. Note: trait has to be "normal" to use this option.}
  \item{\dots}{Extra terms not used.}
}
\value{
  \item{yvalue}{vector of the values of the dependent variable.}
  \item{ncategory}{number of category type if it is non-normal data.}
  \item{envi}{environment effect: TRUE or FALSE.}
  \item{nfixcov}{number of fixed covariates.}
  \item{nrancov}{number of random covariates.}
  \item{fixcoef}{values of the fixed covariate(s) for all individuals.}
  \item{rancoef}{values of the random covariate(s) for all individuals.}
  \item{nran}{number of categories defining the random covariate.} 
  \item{lamda}{value of lamda, the transformation parameter for the \code{boxcox} transformation.}
}
\details{    This function picks the relevant part of the data from the \code{cross}
object and prepares data for \code{\link{qb.model}} and \code{\link{qb.mcmc}}
It can also standardize or transform continuous data if specified.
}
\references{
  \url{http://www.qtlbim.org}}


\note{ This function returns a list and hence should have a differenct name from
that of the \code{cross} object.
}


\author{Dr. Nengjun Yi, et al.,  nyi@ms.ssg.uab.edu}


\seealso{\code{\link{qb.genoprob}}, 
            \code{\link{qb.model}},   \code{\link{qb.mcmc}}}


\examples{
\dontshow{
if(!exists("cross"))
  example(qb.genoprob)

}
qbData <- qb.data(cross, pheno.col = 3, rancov = 2, fixcov = 1)
}
\keyword{data}
