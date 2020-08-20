\name{qb.mcmc}
\alias{qb.mcmc}
\title{Bayesian Mutiple Interacting QTL mapping using MCMC}
\description{
  A computationally efficient MCMC algorithm using the Gibbs sampler or 
  Metropolis-Hastings algorithm is used to produce posterior samples for 
  QTL mapping.   
}
\usage{
qb.mcmc(cross, data, model, mydir = ".", n.iter = 3000, n.thin = 20,
  n.burnin = 0.01*n.iter*n.thin,
  genoupdate = TRUE, seed = 0, verbose = TRUE, \dots)   
}
\arguments{
  \item{cross}{An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.}
  \item{data}{List returned by calling the function \code{\link{qb.data}} .}
  \item{model}{List returned by calling the function \code{\link{qb.model}}.}
  \item{mydir}{A directory to save output from \code{qb.mcmc} in several \file{*.dat} files.
                A directory is created using the trait name and the system time and date.
                If no directory is specified, the default directory is the current working 
                directory.}
  \item{n.iter}{number of iterations to be saved in \code{mydir}, the default being
                3000. Note that, \code{n.iter} is not the total number of         iterations performed but the number iterations saved or considered as posterior
                samples for future analysis. The actual number of iterations would be \code{n.burnin + n.iter*n.thin}                                                
                }
  \item{n.thin}{the thinning number which must be a positive number (default=40)}
  \item{n.burnin}{the initial burn-in period, i.e number of iterations to discard
                  at the beginning of the MCMC run default being 0.01*n.iter*n.thin.}
  \item{genoupdate}{=TRUE will update QTL genotypes and =FALSE will not do so and use
                      the expected value of the QTL genotypes.}
  \item{seed}{Specifies the seed for the random number generator. Using the same seed
              for two runs of the \code{qb.mcmc} function will generate the exact same
              output. The \code{seed} needs to be an integer. The default value for seed
              is the system time.}
  \item{verbose}{=TRUE will force periodic output of the number of MCMC iterations saved.
                 The location of the output directory where results are stored and
                  the time taken for the MCMC run will also be displayed to the user.}
  \item{\dots}{Paramters passed to \code{\link{qb.data}} or
    \code{\link{qb.model}} if \code{data} or \code{model}, respectively,
    is not provided.}
}
   

\value{
  Returns a list of class \code{qb}, including: 
  \item{args}{Arguments passed to \code{qb.mcmc}, \code{\link{qb.data}} and
    \code{\link{qb.model}}. An additional element for \code{subset} may be
    added by \code{\link{subset.qb}} if called.}
  \item{cross.object}{A \code{\link[qtl]{clean}}-ed version of the
    original cross object, but only with phenotypes used by
    \code{qb.mcmc}.}
  \item{mcmc.samples}{A list containing the MCMC samples for each
    phenotype. There are always data frames for \code{iterdiag} and
    \code{mainloci}, with optional data frames for \code{pairloci},
    \code{covariates} and \code{gbye}.}
}

\details{
  A composite model space approach to develop a Bayesian model selection framework for identifying
  interacting QTL for complex traits in experimental crosses from two inbred lines. By placing a liberal
  constraint on the upper bound of the number of detectable QTL we restrict attention to models of 
  fixed dimension. Either Gibbs sampler or Metroplis-Hastings algorithm can be applied to sample
  from the posterior distribution.            
  
  The following data frames in the \code{mcmc.samples} element of the
  \code{qb} object contain the MCMC samples. They are used by many other
  routines.
  
  The iterations data frame \code{iterdiag} has \code{n.iter} rows and
  5 major columns:
  \code{niter} = iteration number;
  \code{nqtl} =  number of putative QTLs included;
  \code{mean} =  overall mean;
  \code{envvar} = residual variance;
  \code{var} = total genetic variance.
  Depending on the type of cross, presence of covariates and epistatic
  effects there would be more columns in the following order:
  \code{varadd} = variance of all additive effects;
  \code{vardom} = variance of all dominant effects;
  \code{varaa} = variance of all additive-additive interactions;
  \code{varad} = variance of all additive-dominant interactions;
  \code{varda} = variance of all dominant-additive interactions;
  \code{vardd} = variance of all dominant-dominant interactions.
  Values for variance of environment-additive interaction,
  variance of environment-dominant interaction, and variance of
  environment effect have names that encode the covariate.

  Covariates are in data frame \code{covariates}, with \code{n.iter} rows and
  L+M(\code{length(fixcov)+length(rancov)}) columns:
  \code{L columns} : Coefficient of the fixed effect.
  \code{M columns} : Variance of the random effect.
  If an ordinal trait is analyzed, the cutoff points for the threshold
  model are also included in additional columns. There would be C-3
  bounded threshold values for an ordinal phenotype with C categories.
                             
  The mainloci data frame has N rows (N=sum of number of QTLs
  detected in \code{n.iter} iterations) and 6-8 columns:
  \code{niter} = iteration number;
  \code{nqtl} =  number of putative QTLs included;
  \code{chrom} = chromosome number;
  \code{locus} = locus in cM;
  \code{add} = additive effect;
  \code{dom} = dominance effect (if included);
  \code{varadd} = variance of additive effect;
  \code{vardom} = variance of dominant effect (if included).
                            
  The pairloci data frame has N rows (N=sum of number
  of pairs of QTLs with epistatic effect detected) and 8-14 columns:
  \code{niter} = iteration number;
  \code{n.epis} =  number of epistatic pairs included;
  \code{chrom1} = first chromosome number;
  \code{locus1} = first locus in cM;
  \code{chrom2} = second chromosome number;
  \code{locus2} = second locus in cM;
  \code{aa} = additive-additive effect;
  \code{ad} = additive-dominant effect (if included);
  \code{da} = dominant-additive effect (if included);
  \code{dd} = dominant-dominant effect (if included);
  \code{varaa} = variance of additive-additive interaction;
  \code{varad} = variance of additive-dominant interaction (if included);
  \code{varda} = variance of dominant-additive interaction (if included);
  \code{vardd} = variance of dominant-dominant interaction (if included).
                        
  The gbye (Gene by Environment) data frame has 7-9 columns:
  \code{niter} = iteration number;
  \code{n.gbye} = number of GxE terms included;
  \code{covar} = fixed covariate identifier;
  \code{chrom} = chromosome number;
  \code{locus} = locus in cM;
  \code{add} = additive effect;
  \code{dom} = dominance effect (if included);
  \code{varadd} = variance of additive effect;
  \code{vardom} = variance of dominant effect (if included).

  The deviance data frame has 1 column with the posterior deviance.
  There is one deviance value for each iteration, or \code{n.iter} values.
  The last value is the deviance calculated at the posterior means, known as Dhat.
}
   
   
\references{
  \url{http://www.qtlbim.org}}
\author{Nengjun Yi, nyi@ms.ssg.uab.edu}

\seealso{\code{\link{qb.sim.cross}},   \code{\link{qb.data}},   
  \code{\link{qb.model}},   \code{\link{qb.mcmc}}}


\examples{
\dontrun{
example(qb.sim.cross)

## Calculate grids and genotypic probabilites.
cross <- qb.genoprob(cross, step=2)

## Create MCMC samples
## First line as qb.data options; second line has qb.model options.
qbExample <- qb.mcmc(cross, pheno.col = 3, rancov = 2, fixcov = 1,
  chr.nqtl = rep(3,nchr(cross)), intcov = 1, interval = rep(10,3),
  n.iter = 1000, n.thin = 20)
}
}
\keyword{models}
\keyword{regression}
