#####################################################################
##
## $Id: mcmc.R,v 1.8 2006/12/15 19:05:40 dshriner Exp $
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
                     n.iter = 3000, n.thin = 20, n.burnin = 0.01*n.iter*n.thin,
                     genoupdate = TRUE,
                     seed = 0, verbose = TRUE,
                     ... )
{
  if(class(cross)[2] != "cross")
     stop("The first input variable is not an object of class cross",call.= FALSE)

  algorithm <- "M-H"

  cross.name <- deparse(substitute(cross))

  ## Save cross object if it is a transient object.
  tmp <- make.names(cross.name)
  is.transient.cross <- tmp != cross.name
  if(is.transient.cross) {
    ## Should make sure the names is unique.
    tmp <- make.dir.name(tmp, "")
    assign(tmp, cross, pos = 1)
    cross.name <- tmp
  }
  
  if(is.null(cross$geno[[1]]$prob)) {
    warning("First running qb.genoprob with default step size",
            call. = FALSE, immediate. = TRUE)
    cross <- qb.genoprob(cross)
  }
  
  n.ind = nind(cross)                   # number of individuals
  n.chr = nchr(cross)                   # number of chromsomes
  n.gen = dim(cross$geno[[1]]$prob)[3]  # number of genotypes
##############################
## For Multiple traits -- SB
  multiple.trait = 0        #logical if TRUE multiple trait analysis is done
  if(data$multiple.trait){ 
  multiple.trait=1
  if(tolower(data$multiple.type)=="traditional") multiple.trait=1 #redundant but being careful
  if(tolower(data$multiple.type)=="sur") multiple.trait=2 
   } 
##########################
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

  ## data 
  pheno = data$pheno.col
  n.pheno=length(pheno) ## For Multiple traits -- SB
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

  if(is.null(data$censor)) {
     censor.lo = rep(-10000000,n.ind)
     censor.hi = rep(10000000,n.ind)
  }
  if(!is.null(data$censor)) {
     censor.lo = data$censor[,1]
     censor.lo[censor.lo==-Inf] = -10000000
     censor.lo[is.na(censor.lo)] = 999
     censor.hi = data$censor[,2]
     censor.hi[censor.hi==Inf] = 10000000
     censor.hi[is.na(censor.hi)] = 999
  }	

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
  contrast = model$contrast
##############################
## For Multiple traits -- SB

  if(multiple.trait==2) diff.loc=model$diff.loc else diff.loc=FALSE

if(multiple.trait)
{ ## Initial Value for variance-cov matrix (Done here since calculation of inverse
  ##  is required
   Y <- yvalue
   Y[Y==999]=NA
  sigma.initial <- var(Y,na.rm=TRUE,use="complete.obs")
  sigma.initial <- solve(sigma.initial)
}  
##########################

  totaliter = n.iter*n.thin + n.burnin

  ## Check that intcov is of length nfixcov.
  intcov <- check.intcov(intcov, nfixcov)
  
  if(verbose) {
    cat(paste("qb.mcmc is running",format(totaliter,big.mark=","),"iterations. The current iterations are saved: \n",sep=" "))     
    start.walltime = Sys.time()
  }

  allTraits = colnames(cross$pheno)
  
########### Single trait Code #################################################

if(!multiple.trait)
{
output = output.dir( qbDir = mydir, traitName = allTraits[pheno] )	# set up a directory to save outputs

  z<-.C("RSingleTraitMCMCSetup",

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
        as.integer(contrast),        # 1: Cockerham model, 0: estimate genotypic values
        as.double(censor.lo),
        as.double(censor.hi),
    
        as.integer(seed),	     # Seed specification for the pseudo-random number generator srand
        as.integer(verbose),
        PACKAGE="qtlbim")

  if(verbose) {
    cat("MCMC sample has been saved to: ")
    cat(output)
    cat(".\n")
  }

  ## calculate pD and DIC for model comparison
  deviances = read.table( file = paste(output,"/deviance.dat",sep="") )
  pD1 = mean(deviances[1:n.iter,1],na.rm=TRUE) - deviances[n.iter+1,1]
  pD2 = 0.5*var(deviances[1:n.iter,1],na.rm=TRUE)
  DIC1 = mean(deviances[1:n.iter,1],na.rm=TRUE) + pD1
  DIC2 = mean(deviances[1:n.iter,1],na.rm=TRUE) + pD2

  ## create an object qb
  qb = c(cross.name = cross.name,
    cross = class(cross)[1],
    output.dir = output, 
    n.iter = n.iter,
    n.thin = n.thin,
    n.burnin = n.burnin,
    algorithm = algorithm, 
    genoupdate = genoupdate,
    pD1 = pD1,
    pD2 = pD2,
    DIC1 = DIC1,
    DIC2 = DIC2,
    step = attr(cross$geno[[1]]$prob, "step"),
    seed = seed,
    verbose = verbose,
    data,
    model ) 


  
}

########### Multiple trait Code #################################################

if(multiple.trait)
{
   yvalue = as.matrix(yvalue)
   dimnames(yvalue) <- NULL
   qtlloc=1
  if( algorithm == "M-H" ) algorithm = 0  else algorithm = 1

output = output.dir( qbDir = mydir, traitName = allTraits[pheno][1] )	# set up a directory to save outputs

  z<-.C("RMultipleTraitsMCMCSetup",

  as.integer(n.ind),           # number of individuals
  as.integer(n.chr),           # number of chromosomes
  as.integer(n.gen),           # number of genotypes
  as.integer(n.pheno),         # number of phenotypes to analyzed jointly
  as.integer(n.loci),          # numbers of loci on chromosomes
  as.double(loci),             # vector of loci
  as.double(prob),             # vector of probabilities

	as.vector(yvalue,mode="double"),                        # phenotypic value
  as.integer(multiple.trait),   # logical if =1 mutiple trait analysis performed
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
  as.integer(diff.loc),
  as.integer(qtlloc),
  PACKAGE="qtlbim")

  if(verbose) {
    cat("MCMC sample has been saved to: ")
    cat(output)
    cat(".\n")
  }

  ## create an object qb
  qb = list(cross.name = cross.name,
    cross = class(cross)[1],
    output.dir = output, 
    n.iter = n.iter,
    n.pheno=n.pheno,
    n.thin = n.thin,
    n.burnin = n.burnin,
    algorithm = algorithm, 
    genoupdate = genoupdate,
    step = attr(cross$geno[[1]]$prob, "step"),
    seed = seed,
    verbose = verbose)

    qb=c(qb,data,model) 

#   for(ph in 1:n.pheno) qb = qb.reorder( qb, pheno = ph )   # transfer mcmc output format for graphics

} 

  ## Assign qb.genoprob attributes to args list if not there.
  defaults <- qb.genoprob.defaults(cross)
  for(i in names(defaults))
    if(is.null(qb[[i]]))
       qb[[i]] <- defaults[[i]]

  if(verbose) {
    stop.walltime = Sys.time()
    walltime = difftime(stop.walltime, start.walltime, units="min")
    cat("Bayesian MCMC took ")
    cat(as.character(round(walltime, digits=2)))
    cat(" minutes. \n")
    qb=c(qb,time=walltime) ## SB add for simulations
    }

  qb = qb.reorder( qb )   # transfer mcmc output format for graphics
  gc()

  class(qb) = "qb"

  ## Reorganize as new qb object. For now do this at the end. Later do as created.
  qb <- qb.legacy(qb, remove = TRUE)

  ## Remove temporary cross if it was from a transient object.
  if(is.transient.cross)
    remove(list = cross.name, pos = 1)
  
  qb
}

qb.model <- function( cross, epistasis = TRUE, 
                       main.nqtl = 3, mean.nqtl = main.nqtl + 3, max.nqtl = NULL, 
                       interval = NULL, chr.nqtl = NULL, 
                       intcov = c(0), depen = FALSE, prop = c(0.5, 0.1, 0.05),
                       contrast = TRUE, 
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
  n.chr = nchr(cross)                   # number of chromsomes
  if(length(interval) == 1)
    interval <- rep(interval, n.chr)
  if(length(interval)<n.chr & !is.null(interval))
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
     if(chr.nqtl[i] > length(cross$geno[[i]]$map)) chr.nqtl[i] = length(cross$geno[[i]]$map)  
     if(chr.nqtl[i] < 0) chr.nqtl[i] = 0      
  }
  if(max.nqtl > sum(chr.nqtl)) max.nqtl = sum(chr.nqtl)

  qtl_envi = FALSE                                   
  if(sum(intcov)!=0) qtl_envi = TRUE 

  model = list( epistasis=epistasis, main.nqtl=main.nqtl, mean.nqtl=mean.nqtl, max.nqtl=max.nqtl,                                                 
                interval=interval, chr.nqtl=chr.nqtl, qtl_envi=qtl_envi, intcov=intcov, 
                depen=depen, prop=prop, contrast=contrast,diff.loc=TRUE )

  gc()
  model
}

#######################################################################################################
make.dir.name <- function(traitName, dtsep = "-")
{
  fullDt <- date()
  dtMinusDay <- substr(fullDt, 5, nchar(fullDt))
  dt <- gsub(":", "", dtMinusDay)
  dt <- gsub(" ", dtsep, dt)
  dt <- substr(dt, 1, nchar(dt) - 5)
  paste(traitName, dt, sep = "_")
}
#######################################################################################################

output.dir <- function( qbDir = getwd(), traitName = "trait1" )
{
 
   if(!file.exists(qbDir)){
      stop("This output directory does not exist.", call. = FALSE)
   }

   traitName <- make.dir.name(traitName)

   traitDir<-file.path(qbDir,traitName)
  
   if(!dir.create(traitDir)){
      stop("This trait specific directory already exists.", call. = FALSE)
   }

   iterFile = file.path(traitDir,"iterdiag.dat")
   covFile = file.path(traitDir,"covariates.dat")
   mainFile = file.path(traitDir,"mainloci.dat")
   pairFile = file.path(traitDir,"pairloci.dat")
   gbyeFile = file.path(traitDir,"gbye.dat")
   devFile = file.path(traitDir,"deviance.dat")
   sigmaFile =file.path(traitDir,"sigma.dat") 
     
   z <- .C("ROutputManager",
           as.character(iterFile),
           as.character(covFile),
           as.character(mainFile),
           as.character(pairFile),
           as.character(gbyeFile),
           as.character(devFile),
           as.character(sigmaFile),
           PACKAGE="qtlbim")
  
   return(traitDir) 
}

