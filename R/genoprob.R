####################################################################################

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

####################################################################################


### NOTE: Depends on R/qtl 1.03.
qb.genoprob <- function(cross, map.function = map.functions, step = 2,
                        tolerance = 1e-6, stepwidth = "variable", ...)
{
  ## This uses calc.genoprob directly.
  ## Call sequence maintained for legacy user code.
  map.functions <- unlist(as.list(formals(calc.genoprob)$map.function)[-1])
  map.function <- match.arg(tolower(map.function[1]), map.function)

  ## First make sure QTL are not too close.
  if(any(sapply(pull.map(cross), function(x)
                ifelse(length(x) > 1, max(diff(x)) < tolerance, FALSE))))
    cross <- jittermap(cross, tolerance)
  
  ## Call R/qtl routine calc.genoprob for calculations.
  cross <- calc.genoprob(cross, step, map.function = map.function,
                         stepwidth = stepwidth,
                         ...)

  ## Add tolerance as attribute.
  attr(cross$geno[[1]]$prob, "tolerance") <- tolerance
  
  return(cross)
}

qb.humanGenoProb <- function(cross,fill.missing=TRUE,cc.col=NULL){
    n.chr = nchr(cross)
    n.ind = nind(cross)
    n.mar = nmar(cross)
    grid = cross$geno
    for(i in 1:n.chr) {
      grid[[i]] = cross$geno[[i]]$map
      n.gen = max(cross$geno[[i]]$data[,1],na.rm=T) # number of genotypes
      cross$geno[[i]]$prob = array(1, c(n.ind,n.mar[i],n.gen))
      for(k in 1:n.gen) cross$geno[[i]]$prob[,,k] = (cross$geno[[i]]$data==k)*1
      colnames(cross$geno[[i]]$prob) = names(cross$geno[[i]]$map)
      if(n.gen==3) dimnames(cross$geno[[i]]$prob)[[3]] = c("AA","AB","BB")
      else {
        if(class(cross$geno[[i]])=="A") dimnames(cross$geno[[i]]$prob)[[3]] = c("AA","AB")
        else dimnames(cross$geno[[i]]$prob)[[3]] = c("g1","g2")
      }
      attr(cross$geno[[i]]$prob,"map") = cross$geno[[i]]$map
      attr(cross$geno[[i]]$prob,"error.prob") = 1e-04
      attr(cross$geno[[i]]$prob,"step") = 0
      attr(cross$geno[[i]]$prob,"off.end") = 0
      attr(cross$geno[[i]]$prob,"map.function") = "haldane"
      attr(cross$geno[[i]]$prob,"stepwidth") = "variable"
      attr(cross$geno[[i]]$prob,"tolerance") = 1e-06
      if(fill.missing) {
        if(is.null(cc.col)|length(unique(na.omit(cross$pheno[,cc.col])))!=2) {
          w = apply(cross$geno[[i]]$data,2,table)
          wsum = apply(w,2,sum)
          for(k in 1:n.gen) {
            z = rep(w[k,]/wsum,each=n.ind)
            na.index = which(is.na(cross$geno[[i]]$prob[,,k]))
            cross$geno[[i]]$prob[,,k][na.index] = z[na.index] 
          }
        }
        if(!is.null(cc.col)&length(unique(cross$pheno[,cc.col]))==2) {
          cc = cross$pheno[,cc.col]
          cc.type = unique(cc)
          cross1 = subset(cross, ind=(cc==cc.type[1]))
          w1 = apply(cross1$geno[[i]]$data,2,table)
          w1sum = apply(w1,2,sum)
          cross2 = subset(cross, ind=(cc==cc.type[2]))
          w2 = apply(cross2$geno[[i]]$data,2,table)
          w2sum = apply(w2,2,sum)
          for(k in 1:n.gen) {
            z1 = rep(w1[k,]/w1sum,each=n.ind)
            z2 = rep(w2[k,]/w2sum,each=n.ind)
            ccrep = rep(cc,n.mar[i])
            na.index1 = which(is.na(cross$geno[[i]]$prob[,,k])&ccrep==cc.type[1])
            na.index2 = which(is.na(cross$geno[[i]]$prob[,,k])&ccrep==cc.type[2])
            cross$geno[[i]]$prob[,,k][na.index1] = z1[na.index1] 
            cross$geno[[i]]$prob[,,k][na.index2] = z2[na.index2] 
          }
        }
      }
    }
    cross
}

calcPriorGenoProb<-function(genoVec){
     genoProb<-table(genoVec)/sum(table(genoVec))
     genoProb
}
