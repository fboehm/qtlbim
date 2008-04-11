#####################################################################
##
##
## Copyright (C) 2006 Nengjun Yi and Tapan Mehta
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
qb.data <- function( cross, pheno.col = 1, trait = c("normal","binary","ordinal"),
                      censor = NULL, 
                      fixcov = c(0),rancov = c(0), 
                      boxcox = FALSE, standardize = FALSE, ... )                   
{
  trait <- match.arg(trait)
  
  qb.valid.phenoData(cross,pheno.col,trait,fixcov,rancov,censor)
  
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
  ## Change missing to 999.
  yvalue[is.na(yvalue)] = 999
  ## Change infinite to 999
  yvalue[yvalue == Inf | yvalue == -Inf] = 999

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

  if(nfixcov!=0){              
  fixcoef = as.data.frame(cross$pheno[,fixcov]) # input fixed covariates   
  ## Change infinite to 999
  fixcoef[fixcoef == Inf | fixcoef == -Inf] = NA
    for(i in 1:nfixcov)          # recode the random covariates to 0,1,2,...
      if(is.factor(fixcoef[[i]]))
      {
       levels( fixcoef[[i]] ) = c(1:length(levels(fixcoef[[i]])))
       fixcoef[[i]] = as.numeric(fixcoef[[i]])-1
      }
    fixcoef[is.na(fixcoef)] = 999
    fixcoef = as.matrix(fixcoef)   
  } else fixcoef = as.matrix(cross$pheno[,fixcov]) # input fixed covariates   
 
  rancoef = as.matrix(cross$pheno[,rancov])  # input random covariates 
  ## Change infinite to 999
  rancoef[rancoef == Inf | rancoef == -Inf] = 999
  if(nrancov!=0) {
    rancoef[rancoef==999] = NA
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

  if(trait=="binary"|trait=="ordinal") censor = NULL
  if( !is.null(censor) ) {
    if ( length(dim(censor))!=2 | dim(censor)[1]!=nind(cross) | dim(censor)[2]!=2 )
       stop("censor should be n x 2 matrix",call.= FALSE)
  }

  data = list( pheno.col=pheno.col, yvalue=yvalue, trait=trait, censor=censor, ncategory=ncategory, 
               envi=envi, nfixcov=nfixcov, nrancov=nrancov, fixcoef=fixcoef, rancoef=rancoef, nran=nran, 
               boxcox=boxcox, lamda=lamda, standardize=standardize,
               covar = c(fixcov, rancov) ) 
  gc()
  data
}





## function to check the input arguments
qb.valid.phenoData<-function(cross,pheno.col,trait,fixcov,rancov,censor){
    if(class(cross)[2] != "cross")
     stop("The first input variable is not an object of class cross",call.= FALSE)
    if(is.na(as.integer(pheno.col)))
     stop("Phenotype column should be an integer",call.=FALSE)
    if(length(na.omit(as.integer(fixcov)))<length(fixcov))
     stop("Fixed covariate column id should be a valid integer",call.=FALSE)
    if(length(na.omit(as.integer(rancov)))<length(rancov))
     stop("Random covariate column id should be a valid integer",call.=FALSE)
    if(length(trait)!= 1 | !(trait %in% c("normal","binary","ordinal")))
     stop("Check the type of trait",call.=FALSE)
    if( !is.null(censor) ) {
    if (length(dim(censor))!=2 | dim(censor)[1]!=nind(cross) | dim(censor)[2]!=2 )
       stop("censor should be n x 2 matrix",call.= FALSE)
    }
}

