#####################################################################
##
## $Id: mcmc.R,v 1.4.2.1 2006/09/07 21:15:39 byandell Exp $
##
##     Copyright (C) 2006 Nengjun Yi and Tapan Mehta
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

##############################################################################

qb.mcmc <- function(cross,
                     data = qb.data(cross, ...),
                     model = qb.model(cross, ...),
                     mydir = ".", 
                     n.iter = 3000, n.thin = 40, n.burnin = 0.01*n.iter*n.thin,
                     algorithm = c("M-H","Gibbs"), genoupdate = TRUE,
                     seed = 0, verbose = TRUE,
                     ... )
{
  if(class(cross)[2] != "cross")
     stop("The first input variable is not an object of class cross",call.= FALSE)

  algorithm <- algorithm[1]
  
  n.ind = nind(cross)           # number of individuals
  n.chr = nchr(cross)           # number of chromsomes
  n.gen = 2                     # number of genotypes
  if(class(cross)[1]=="f2") n.gen = 3

  cross.name <- deparse(substitute(cross))
  
  if(is.null(cross$geno[[1]]$prob)) {
    warning("First running qb.genoprob with default step size",
            call. = FALSE, immediate. = TRUE)
    cross <- qb.genoprob(cross)
  }

  ## loci on the genome
  loci <- pull.loci(cross)
  n.loci <- sapply(loci, length) # number of loci on each chr
  loci <- unlist(loci)

  ## Check for mismatch in loci and prob.
  tmp <- unlist(lapply(cross$geno, function(x) dim(x$prob)[2]))
  if(any(tmp != n.loci))
    stop("Please first run qb.genoprob")

  ## probabilities at all loci
  prob <- unlist(lapply(cross$geno, function(x) x$prob))

  if( algorithm == "M-H" ) algorithm = 0
  else algorithm = 1

  ## data 
  pheno = data$pheno.col
  yvalue = data$yvalue
  if( data$trait == "normal" ) traittype = 1
  if( data$trait == "binary" ) traittype = 2
  if( data$trait == "ordinal" ) traittype = 3
  ncategory = data$ncategory
  nrancov =  data$nrancov
  nfixcov =  data$nfixcov                               
  envi = data$envi                                                      
  fixcoef = data$fixcoef
  rancoef = data$rancoef
  nran = data$nran	

  ## model
  epis = model$epistasis
  emainqtl = model$main.nqtl
  eqtl = model$mean.nqtl
  mnqtl = model$max.nqtl
  interval = model$interval
  chrnqtl = model$chr.nqtl
  qtl_envi = model$qtl_envi                                   
  intcov = model$intcov
  depen = model$depen
  prop = model$prop

  if(verbose) {
    cat("Bayesian MCMC run in progress. The current saved iterations: \n")     
    start.walltime = Sys.time()
  }

  allTraits = colnames(cross$pheno)
  output = output.dir( qbDir = mydir, traitName = allTraits[pheno] )	# set up a directory to save outputs

  z<-.C("R_AnalysisEngine",

        as.integer(n.ind),           # number of individuals
        as.integer(n.chr),           # number of chromosomes
        as.integer(n.gen),           # number of genotypes
        as.integer(n.loci),          # numbers of loci on chromosomes
        as.double(loci),             # vector of loci
        as.double(prob),             # vector of probabilities

	as.double(yvalue),           # phenotypic value
	as.integer(traittype),       # type of traits (1: normal, 2: binary, 3: ordinal)
        as.integer(ncategory),       # number of categories for binary and ordinal traits    	
        as.integer(n.iter),          # number of iterations
	as.integer(n.thin),          # number of thin
        as.integer(n.burnin),        # number of burnin
        as.integer(algorithm),       # MCMC method (1: GIBBS, 0: MH sampler)
	as.integer(genoupdate),      # 1: update QTL genotypes; 0: don't update QTL genotypes
        as.integer(epis),            # epistatic model(epis=1) or nonepistatic model (epis=0)
	as.integer(emainqtl),        # prior number of main-effect QTL
	as.integer(eqtl),            # prior number of all QTL
	as.integer(mnqtl),           # upper bound of detectable QTL
	as.double(interval),         # distance of flanking QTL
	as.integer(chrnqtl),         # maximum number of QTL on each chromosome
        as.integer(envi),            # envi=1: include covariates, envi=0: covariates
        as.integer(qtl_envi),        # qtl_envi=1: include g-by-e, qtl_envi=0: exclude g-by-e 
	as.integer(nrancov),         # number of random covariates
	as.integer(nfixcov),         # number of fixed covariates
	as.integer(intcov),          # indicating which fixed covariates are included to interact QTL 
        as.double(rancoef),          # random covariates
        as.double(fixcoef),          # fixed covariates
        as.integer(nran),            # numbers of levels for random covariates 
	as.integer(depen),           # depen=1: use dependent prior for epistatic effects
	as.double(prop),             # prior prob for three types of epistatic effects
        as.integer(seed),	     # Seed specification for the pseudo-random number generator srand
        as.integer(verbose),
        PACKAGE="qtlbim")

  if(verbose) {
    cat("MCMC sample has been saved to: ")
    cat(output)
    cat(".\n")
  }

  ## create an object qb
  qb = c(cross.name = cross.name,
    cross = class(cross)[1],
    output.dir = output, 
    n.iter = n.iter,
    n.thin = n.thin,
    n.burnin = n.burnin,
    algorithm = algorithm, 
    genoupdate = genoupdate,
    step = attr(cross$geno[[1]]$prob, "step"),
    seed = seed,
    verbose = verbose,
    data,
    model ) 

  if(verbose) {
    stop.walltime = Sys.time()
    walltime = difftime(stop.walltime, start.walltime, units="min")
    cat("Bayesian MCMC took ")
    cat(as.character(round(walltime, digits=2)))
    cat(" minutes. \n")
  }

  qb = qb.reorder( qb )   # transfer mcmc output format for graphics

  gc()

  class(qb) = "qb"

  qb
}

###############################################################################################################

qb.data <- function( cross, pheno.col = 1, trait = c("normal","binary","ordinal"), 
                      fixcov = c(0),rancov = c(0), 
                      boxcox = FALSE, standardize = FALSE, ... )                   
{
  if(class(cross)[2] != "cross")
     stop("The first input variable is not an object of class cross",call.= FALSE)
  

  ## Make sure trait has only one value.
  trait <- trait[1]	
  
  yvalue = cross$pheno[,pheno.col]

  lamda = NULL
  if( boxcox & ( length(yvalue[yvalue>0])!=length(yvalue) | trait!="normal" ) )
    stop("The boxcox transformation cannot be used for this data")
  if(boxcox & trait=="normal") {
     require("MASS")
     if( length(yvalue[yvalue>0])==length(yvalue) ) {   
         BC = boxcox(yvalue ~ 1)
         lamda = BC$x[ which.max(BC$y) ] 
         if(lamda != 0) yvalue = (yvalue^lamda - 1)/lamda
         else
             yvalue = log10(yvalue)
     }
  }
  if(standardize & trait=="normal")
    yvalue = (yvalue - mean(yvalue, na.rm=T))/sd(yvalue, na.rm=T)
  yvalue[is.na(yvalue)] = 999

  if(trait=="normal") ncategory = 0	
  if(trait!="normal") {
     yvalue[yvalue==999] = NA
     yvalue = factor(yvalue)
     levels(yvalue) = c( 1:length(unique(yvalue)) )
     yvalue = as.numeric(yvalue) - 1        # recode the phenotype to 0,1,2,...
     yvalue[is.na(yvalue)] = 999
     ncategory = max(yvalue[yvalue!=999]) + 1  # the number of categories of ordinal traits 
  }   

  nrancov = 0; nfixcov = 0                                
  if(rancov[1]!=0) nrancov = length(rancov)   # number of covariates can be calculated from rancov and fixcov
  if(fixcov[1]!=0) nfixcov = length(fixcov)
  envi = FALSE                                        
  if((nrancov+nfixcov)!=0) envi = TRUE              

               
  fixcoef = as.matrix(cross$pheno[,fixcov]) # input fixed covariates   
  fixcoef[is.na(fixcoef)] = 999

  rancoef = as.matrix(cross$pheno[,rancov])  # input random covariates 
  rancoef[rancoef==999] = NA
  if(nrancov!=0) {
    for(i in 1:nrancov)          # recode the random covariates to 0,1,2,...
    {
       rancoef[ ,i] = factor(rancoef[ ,i])
       levels( rancoef[ ,i] ) = c(1:length(unique(rancoef[ ,i])))
       rancoef[ ,i] = as.numeric(rancoef[ ,i])-1
    }
    rancoef[is.na(rancoef)] = 999
  }

  nran = rep(0, nrancov)
  if(nrancov!=0) {
    for(i in 1:nrancov)       # calculate the number of each random effects
    {
       x = rancoef[ ,i]
       x = x[x!=999]
       nran[i] = max(x) + 1  
    } 
  }

  data = list( pheno.col=pheno.col, yvalue=yvalue, trait=trait, ncategory=ncategory, 
               envi=envi, nfixcov=nfixcov, nrancov=nrancov, fixcoef=fixcoef, rancoef=rancoef, nran=nran, 
               boxcox=boxcox, lamda=lamda, standardize=standardize,
               covar = c(fixcov, rancov) ) 
  gc()
  data
}

#############################################################################################################

qb.model <- function( cross, epistasis = TRUE, 
                       main.nqtl = 3, mean.nqtl = main.nqtl + 3, max.nqtl = NULL, 
                       interval = NULL, chr.nqtl = NULL, 
                       intcov = c(0), depen = FALSE, prop = c(0.5, 0.1, 0.05), 
                       ... )	   
                     
{
  if(class(cross)[2] != "cross")
     stop("The first input variable is not an object of class cross",call.= FALSE)
	

  if(epistasis == FALSE) mean.nqtl = main.nqtl
  
  if(is.null(max.nqtl)) {
       if(epistasis == FALSE) max.nqtl = main.nqtl+3*main.nqtl^0.5 
       else max.nqtl = mean.nqtl+3*mean.nqtl^0.5
  } 
  max.nqtl = as.integer(max.nqtl) 
  
  if(is.null(interval)) interval = unlist( lapply(cross$geno, function(x) mean(diff(x$map))) ) 
  if(length(interval)<nchr(cross) & !is.null(interval))
     stop("You should specify interval for each chromosome") 
  
  if(length(chr.nqtl)<nchr(cross) & !is.null(chr.nqtl))
     stop("You should specify chr.nqtl for each chromosome") 
  if(is.null(chr.nqtl)) {
     chr.nqtl = unlist( lapply(cross$geno, function(x) diff(range(x$map))) )
     for(i in 1:nchr(cross)) {
        if(interval[i] != 0) chr.nqtl[i] = chr.nqtl[i]/interval[i]
     }
  }
  chr.nqtl = as.integer(chr.nqtl)

  for(i in 1:nchr(cross)) {
     if( interval[i] > diff(range(cross$geno[[i]]$map)) ) interval[i] = diff(range(cross$geno[[i]]$map))
     if(interval[i] != 0) {
     	d = diff(range(cross$geno[[i]]$map))/interval[i]   
     	if( d < chr.nqtl[i] ) chr.nqtl[i] = as.integer(d)
     }
     if(chr.nqtl[i] > length(cross$geno[[1]]$map)) chr.nqtl[i] = length(cross$geno[[1]]$map)  
     if(chr.nqtl[i] < 0) chr.nqtl[i] = 0      
  }
  if(max.nqtl > sum(chr.nqtl)) max.nqtl = sum(chr.nqtl)

  qtl_envi = FALSE                                   
  if(sum(intcov)!=0) qtl_envi = TRUE 

  model = list( epistasis=epistasis, main.nqtl=main.nqtl, mean.nqtl=mean.nqtl, max.nqtl=max.nqtl,                                                 
                interval=interval, chr.nqtl=chr.nqtl, qtl_envi=qtl_envi, intcov=intcov, 
                depen=depen, prop=prop )

  gc()
  model
}

#######################################################################################################

output.dir <- function( qbDir = getwd(), traitName = "trait1" )
{
 
   if(!file.exists(qbDir)){
      stop("This output directory does not exist.", call. = FALSE)
   }

#   dt = substring(date(),5,19)
#   dt = gsub(" ","-",dt)
   fullDt = date()
   dtMinusDay = substr(fullDt,5,nchar(fullDt))
   dt = gsub(":","",dtMinusDay)
   dt = gsub(" ","-",dt)
   dt = substr(dt,1,nchar(dt)-5)
   traitName = paste(traitName,dt,sep="_")

   traitDir<-file.path(qbDir,traitName)
  
   if(!dir.create(traitDir)){
      stop("This trait specific directory already exists.", call. = FALSE)
   }

   iterFile = file.path(traitDir,"iterdiag.dat")
   covFile = file.path(traitDir,"covariates.dat")
   mainFile = file.path(traitDir,"mainloci.dat")
   pairFile = file.path(traitDir,"pairloci.dat")
   gbyeFile = file.path(traitDir,"gbye.dat")
     
   z <- .C("R_OutputManager",
           as.character(iterFile),
           as.character(covFile),
           as.character(mainFile),
           as.character(pairFile),
           as.character(gbyeFile),
           PACKAGE="qtlbim")
  
   return(traitDir) 
}

###########################################################################################################

