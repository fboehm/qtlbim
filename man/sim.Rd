\name{qb.sim.cross}
\alias{sim.data}
\alias{qb.sim.cross}
\alias{summary.qb.sim}
\title{Simulates QTL related data for an F2 or BC cross.}
\description{This function is used to simulate genotypic, phenotypic and
  covariate data for BC and F2 populations. The underlying genetic model
  is Cockerham's model and data for both continuous (normally
  distributed only) and ordinal traits can be generated.
}
\usage{
qb.sim.cross(len = rep(100,20), n.mar = 11, eq.spacing = TRUE,
  n.ind = 400, type = c("f2","bc"), missing.geno = 0.0,
  missing.pheno = 0.0, ordinal = c(0.5,0.5),
  qtl.pos = NULL, qtl.main = NULL, qtl.epis = NULL,
  covariate = NULL, gbye = NULL, seed = NULL )
\method{summary}{qb.sim}(object, \dots)
} 
\arguments{
  \item{len}{ defines the length (in cM) of each chromosome and number of
    chromosomes.  Thus \code{len = c(80,90,44)} would represent a model with three
     chromosomes of lengths 80, 90, and 44 respectively.}
  \item{n.mar}{The number of markers per chromosome.
                                  This can be specified as a single number or as a vector.
                                  If a single number is specified, all the chromosomes will have the same
                                  number of markers.  If n.mark is a vector then it must have
                                  the same number of entries as there are chromosomes.  For
                                  example, if n.mar = c(10,11,9) then we have a three chromosome
                                  model in which the first chromosome has 10 markers, the second
                                  has 11 and the third has 9.a vector specifying the number of markers per chromosome.}
  \item{eq.spacing}{if TRUE, markers will be equally spaced. Default is TRUE. 
                    If FALSE, markers are generated uniformly over the chromosome.}
  \item{n.ind}{specifies the number of individuals.}
  \item{type}{indicates whether to simulate an intercross ("f2") or a backcross ("bc").}
  \item{missing.geno}{the frequency of missing genotypes.}
  \item{missing.pheno}{the frequency of missing phenotypic values.}
  \item{ordinal}{define the probabilities of each ordinal category and the number
                  of elements in the vector will determine the number of categories.
                  The elements must be positive and the should sum up to 1}

  \item{qtl.pos}{
                This parameter specifies the positions
		of qtl as a matrix with
                dimensions (number of qtl) x 2.  Note
		that the row dimension is the number
		of qtl and is not the number of
		chromosomes.  Each row identifies a qtl, the
                first column entries represent the chromosome's index, the second
                column entries represent the location on the chromosome
		of the qtl.   The (row) order in which qtl are listed in this parameter
               is the index by which they are identified later on in the parameters
               \code{qtl.main} and \code{qtl.epi}.}
  \item{qtl.main}{The parameter \code{qtl.main}
    is a matrix specifying the main effects of QTLs.  The first
    column gives the qtl-index (the row index of the qtl in the \code{qtl.pos}
    parameter.), the second and third column gives the additive and dominance effects,
    respectively. There are two or three columns depending on \code{type} being "bc" or "f2".}
  \item{qtl.epis}{
    It is a matrix specifying epistatic effects. There are 3 or 6 columns 
    depending on \code{type} being "bc" or "f2". Each row gives an
    epistatic pair.  The first entry in a row gives the first qtl index,
    the second entry represents the index of the second qtl. The other
    entries give the value of the epistatic effects (additive-additive,
    additive-dominance, dominance-additive and dominance-dominance) 
    of the two qtls.  The indices used to represent the qtl are the row
    indices of the \code{qtl.pos} matrix which correspond to the first
    and second qtl in each epistatic pair.}
    
  \item{covariate}{
  A vector of two elements, the first being the true value of the coffecient 
  for the fixed covariate and the second the true value for the standard deviation 
  of the random covariate.
  } 
  \item{gbye}{
    A matrix specifying the interaction between the fixed covariate and QTL main
    effect. The first column is the index of the QTL, the other column(s) is(are)
    the value(s) of interaction(s).
  }
  \item{seed}{Set pseudo-random number seed with
    \code{\link[base]{set.seed}} if not \code{NULL}.}
  \item{object}{An object of class \code{qb.sim}, typically the
    \code{qtl} element of a \code{cross} object created by
    \code{qb.sim.cross}.
  }
  \item{\dots}{Not used here.}
  
}
\value{
  \code{qb.sim.cross} will returns an object of class \code{cross}.
  See \code{\link[qtl]{read.cross}} for details. In addition, a component
  \code{qtl} of class \code{qb.sim} is added which is a list of atmost 6
  components depending on the options specified.
   
          \item{geno}{is a matrix of true QTL genotypes for every
                                     individual and each locus. The genotypes are
                                     defined following \code{\link[qtl]{read.cross}}.}
                                     
          \item{pos}{is a matrix of true QTL position. Same as \code{qtl.pos}.}
          
          \item{herit.main}{is a matrix of the heritability of main effects.
                                          \code{nrow($qtl$herit.main)}=no. of QTLs and 
                                          \code{ncol($qtl$herit.main)}=2 or 3 depending on the type
                                          of genetic cross ("bc" or "f2"). The first column being the
                                          QTL index and the others being additive and dominant heritability
                                          respectively.}
          \item{herit.epis}{is a matrix of the heritability of epistatic effects.
                                          \code{nrow($qtl$herit.epis)}=no. of QTLs pairs interacting and 
                                          \code{ncol($qtl$herit.main)}=3 or 6 depending on the type
                                          of genetic cross ("bc" or "f2"). The first column being the
                                          QTL index and the others being additive-additive, 
                                          additive-dominant, dominant-additive and dominant-dominant
                                          heritability respectively.}
          
          \item{herit.cov}{is a vector of length 2 containing the heritability
                                           of the fixed and random covariate.}
          
          \item{herit.gbye}{is a matrix of heritability of GxE interactions.
                                          \code{nrow($qtl$herit.gbye)}= no. of GxE interactions and
                                          \code{ncol($qtl$herit.gbye)}= 2 or 3 depending on the type
                                          of genetic cross ("bc" or "f2"). The first column being the
                                          GxE index and the others being additive and dominant GxE 
                                          interaction heritability.}
          
          
  }
\details{
The most important difference of this simulation function from others
 is that it computes phenotype values with full genetic model. i.e. both additive, 
 dominance, and espistatic effects are considerred. Furthermore, environmental 
 effects and gene-environment interactions can be included in the model to 
 simulate phenotypes. 

The outputted genotypes for markers and qtls will be coded as 1 and 2 for BC and 
1,2, and 3 for F2. Missing data will be coded as NA.}

\references{
  \url{http://www.ssg.uab.edu/qtlbim}
  }
\author{Dr. Nengjun Yi, et al.,  nyi@ms.ssg.uab.edu}

\seealso{\code{\link{qb.genoprob}}, \code{\link{qb.data}} 
  \code{\link{qb.model}}, \code{\link{qb.mcmc}},
  \code{\link[qtl]{sim.cross}}
}


\examples{
## Simulate F2 cross.
cross <- qb.sim.cross(len = rep(100,20), n.mar = 11, eq.spacing =FALSE,
  n.ind = 500, type = "f2", ordinal = c(0.3,0.3,0.2,0.2),
  missing.geno = 0.03, missing.pheno = 0.03,
  qtl.pos = rbind(qtl.1=c(chr=1,pos=15),qtl.2=c(1,45),qtl.3=c(3,12),
    qtl.4=c(5,15),qtl.5=c(7,15),qtl.6=c(10,15),qtl.7=c(12,35),qtl.8=c(19,15)),
  qtl.main = rbind(main.1=c(qtl=1,add=0.5,dom=0),main.2=c(2,0,0.7),
    main3=c(3,-0.5,0),main4=c(4,0.5,-0.5)),
  qtl.epis = rbind(epis1=c(qtl.a=4,qtl.b=5,aa=-0.7,ad=0,da=0,dd=0),
    epis2=c(6,8,0,1.2,0,0)),
  covariate = c(fix.cov=0.5,ran.cov=0.07),
  gbye = rbind(GxE.1=c(qtl=7,add=0.8,dom=0)) ) 

## Examine simulation information.
summary(cross$qtl)
}
\keyword{datagen}

