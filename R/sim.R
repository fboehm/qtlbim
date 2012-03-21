####################################################################################
##     Copyright (C) 2006 Samprit Banerjee and Nengjun Yi
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


sim.data <- function(...) qb.sim.cross(...)

qb.sim.cross <- function(len = rep(100,20), n.mar = 11, eq.spacing = TRUE,
                         n.ind = 400, type=c("f2","bc","riself","risib"), 
                         missing.geno = 0.0, missing.pheno = 0.0, ordinal = c(0.5,0.5),
                         qtl.pos = NULL, qtl.main = NULL, qtl.epis = NULL,
                         covariate = NULL, gbye = NULL,seed=NULL) 
{  
    type <- tolower(type[1])

if(!is.null(seed)) set.seed(seed)
multiple.trait <- FALSE;n.pheno <- 1;sigma<-diag(n.pheno);if(!multiple.trait & n.pheno > 1) print("Simulating independent phenotypes")

if(n.pheno==1)
   if(multiple.trait) warning("Ignoring option multiple.trait=TRUE")

if(!is.null(qtl.pos)){
    if(is.vector(qtl.pos)) qtl.pos = t(as.matrix(qtl.pos))
}

if(!is.null(gbye)){
if(is.vector(gbye)) gbye = t(as.matrix(gbye))
}

if(!is.null(qtl.main)){
if(is.vector(qtl.main)) qtl.main = t(as.matrix(qtl.main))
}
if(!is.null(qtl.epis)){
if(is.vector(qtl.epis)) qtl.epis = t(as.matrix(qtl.epis))
}
                      
# generate marker map

    if(length(len) != length(n.mar) && length(len) != 1 && length(n.mar) != 1) 
       stop("Lengths of vectors len and n.mar do not conform.")
    if(length(len) == 1) 
       len <- rep(len, length(n.mar))
    if(length(n.mar) == 1) 
      n.mar <- rep(n.mar, length(len))
    n.chr <- length(n.mar)
    map <- vector("list", n.chr)
    names(map) <- as.character(1:n.chr)
    
    
    
    for(i in 1:n.chr) {
       map[[i]] <-c(0,len[i])
       if(n.mar[[i]]>2) {
          if(!eq.spacing) map[[i]] <- sort(c(map[[i]],runif(n.mar[i]-2,0,len[i])))
          else map[[i]] <- seq(0,len[i],length=n.mar[i])
       }     
       names(map[[i]]) <- paste("C",names(map)[i],"M",1:n.mar[i],sep = "")
    }   

# insert QTL into marker map
         
    n.qtl <- 0
    if(!is.null(qtl.pos)) n.qtl <- nrow(qtl.pos) 

    if(n.qtl>0) {
      dimnames(qtl.pos) <- list(paste("qtl", seq(n.qtl), sep = "."),
                                c("chr", "pos"))
       for(i in 1:n.qtl) {
          temp <- map[[qtl.pos[i, 1]]]
          j <- max((seq(along = temp))[temp < qtl.pos[i, 2]])
          temp <- c(temp[1:j], qtl.pos[i, 2], temp[(j + 1):length(temp)])
          names(temp)[j + 1] <- paste("QTL", i, sep = "")
          map[[qtl.pos[i, 1]]] <- temp
      }
    }
    
 
# generate marker genotypes and QTL genotypes
   
    geno <- vector("list", n.chr)
    names(geno) <- names(map)
    n.mar <- sapply(map,length)     # number of markers and QTL
    mar.names <- lapply(map,names)
    
    for (i in 1:n.chr) {
        data <- matrix(nrow=n.ind,ncol=n.mar[i])
        dimnames(data) <- list(NULL,mar.names[[i]])
        d <- diff(map[[i]])
        r <- 0.5*(1-exp(-2*d/100))  # we only consider Haldane distance
        rbar <- 1-r
        if(type=="bc") {
           data[ ,1] <- sample(1:2, n.ind, replace = TRUE)
           if(n.mar[i] > 1) {
              for (j in 1:(n.mar[i]-1)) {
                  rec <- sample(0:1, n.ind, replace = TRUE, prob = c(1-r[j],r[j]))
                  data[rec==0,j+1] <- data[rec==0,j]
                  data[rec==1,j+1] <- 3-data[rec==1,j]
              }
            }
        }
        if(type=="f2") {
           data[, 1] <- sample(1:3, n.ind, replace = TRUE, prob=c(1,2,1))
           if(n.mar[i] > 1) {
               for(j in 1:(n.mar[i] - 1)) {
                   data[data[, j] == 1, j + 1] <- sample(1:3, sum(data[, j]==1), replace = TRUE, 
                        prob=c(rbar[j]*rbar[j], 2*r[j]*rbar[j], r[j]*r[j]))
                   data[data[, j] == 2, j + 1] <- sample(1:3, sum(data[, j]==2), replace = TRUE, 
                        prob = c(r[j]*rbar[j], rbar[j]*rbar[j] + r[j]*r[j], r[j]*rbar[j]))
                   data[data[, j] == 3, j + 1] <- sample(1:3, sum(data[, j]==3), replace = TRUE, 
                        prob = c(r[j]*r[j], 2*r[j]*rbar[j], rbar[j]*rbar[j]))
                }
            }
        }

        geno[[i]] <- list(map=map[[i]],data=data)
        class(geno[[i]]) <- "A"
        class(geno[[i]]$map) <- NULL
    }

# get QTL genotypes from geno[[i]]$data

    if(n.qtl > 0) {
        qtl.chr <- qtl.loc <- NULL
        for (i in 1:n.chr) {
            o <- grep("^QTL[0-9]+", mar.names[[i]])  
            if (length(o) > 0) {
                qtl.chr <- c(qtl.chr, rep(i, length(o)))
                qtl.loc <- c(qtl.loc, o)
            }
        }
        qtl.geno <- matrix(nrow=n.ind,ncol=n.qtl) 
        for(i in 1:n.qtl) qtl.geno[ ,i] <- geno[[qtl.chr[i]]]$data[ ,qtl.loc[i]]
    }
# generate covariates and normal phenotypes

  mu = rep(10,n.pheno)
#  Sigma=sigma %x% diag(n.ind) 
require(MASS)
    pheno.normal <- mvrnorm(n.ind,mu,sigma)   # mu=10, ve=1 
    residual <- pheno.normal - mu

    # Main effects added
if (n.qtl > 0) {
        if(!is.null(qtl.main)) {
           n.main <- nrow(qtl.main)
if(type=="bc" | type=="riself" | type=="risib") {
           if(multiple.trait)
            {
              dimnames(qtl.main) <- list(paste("main", seq(n.main), sep = "."),
                                         c("pheno","qtl", "add"))

              for(m in unique(qtl.main[,"pheno"])) {
#                tmp.main = rbind(qtl.main[qtl.pos[qtl.main[,1],1]==m,])
                tmp.main = rbind(qtl.main[qtl.main[,"pheno"]==m,])
                for(i in 1:nrow(tmp.main)) {
                 xa <- qtl.geno[ ,tmp.main[i,"qtl"]]-1.5
               pheno.normal[,m] <- pheno.normal[,m]+xa*tmp.main[i,3]
                 }
                }
              }
            if(!multiple.trait)
              {  
              dimnames(qtl.main) <- list(paste("main", seq(n.main), sep = "."),
                                         c("qtl", "add"))
                 for(i in 1:n.main) {
                 xa <- qtl.geno[ ,qtl.main[i,1]]-1.5
                 pheno.normal <- pheno.normal+xa*qtl.main[i,2]
                  }
              }
              
           }

           if(type=="f2") {

           if(multiple.trait) 
              {
              dimnames(qtl.main) <- list(paste("main", seq(n.main), sep = "."),
                                         c("pheno","qtl","add","dom"))
              for(m in unique(qtl.main[,"pheno"])) {
                tmp.main = rbind(qtl.main[qtl.main[,"pheno"]==m,])
                for(i in 1:nrow(tmp.main)) {
                  xa <- qtl.geno[ ,tmp.main[i,"qtl"]]-2
                  xd <- (1+xa)*(1-xa)-0.5
           pheno.normal[,m] <- pheno.normal[,m]+xa*tmp.main[i,"add"]+xd*tmp.main[i,"dom"]
                 }
                }
              }
            if(!multiple.trait)
              {  
              dimnames(qtl.main) <- list(paste("main", seq(n.main), sep = "."),
                                         c("qtl", "add","dom"))
               for(i in 1:n.main) {
                  xa <- qtl.geno[ ,qtl.main[i,"qtl"]]-2
                  xd <- (1+xa)*(1-xa)-0.5
                  pheno.normal <- pheno.normal+xa*qtl.main[i,"add"]+xd*qtl.main[i,"dom"]
                   }                 
              }
              
           }
        }
        
        
        
if(!is.null(qtl.epis)) {
           n.epis <- nrow(qtl.epis)   
  if(type=="bc" | type=="riself" | type=="risib"){
    if(multiple.trait)
    {
     dimnames(qtl.epis) <- list(paste("epis", seq(n.epis), sep = "."),
                                   c("pheno","qtl.a", "qtl.b", "aa"))
        for(m in unique(qtl.epis[,"pheno"]))
            {                          
            tmp.epis=NULL               
                 for(k in 1:nrow(qtl.epis)){
                    if(qtl.epis[k,"pheno"]==m) 
                          tmp.epis=rbind(tmp.epis,qtl.epis[k,])
                       }
             for(i in 1:nrow(tmp.epis)) 
                  {
                  xa1 <- qtl.geno[ ,tmp.epis[i,"qtl.a"]]-1.5
                  xa2 <- qtl.geno[ ,tmp.epis[i,"qtl.b"]]-1.5
                  pheno.normal[,m] <- pheno.normal[,m]+xa1*xa2*tmp.epis[i,"aa"] 
                  }  
            }
        }
        if(!multiple.trait)
          {
              dimnames(qtl.epis) <- list(paste("epis", seq(n.epis), sep = "."),
                                         c("qtl.a", "qtl.b", "aa"))
                       for (i in 1:n.epis) {
                  xa1 <- qtl.geno[, qtl.epis[i, "qtl.a"]] - 1.5
                  xa2 <- qtl.geno[, qtl.epis[i, "qtl.b"]] - 1.5
                  pheno.normal <- pheno.normal + xa1 * xa2 * qtl.epis[i, "aa"]
                }
            }
           
         }  
           
           
    if(type=="f2") {
        if(multiple.trait)
          {   
              dimnames(qtl.epis) <- list(paste("epis", seq(n.epis), sep = "."),
                              c("pheno","qtl.a", "qtl.b", "aa", "ad", "da", "dd"))
              for(m in unique(qtl.epis[,"pheno"]))
                  {                          
                    tmp.epis=NULL               
                    for(k in unique(qtl.epis[,"pheno"]))
                         tmp.epis=rbind(tmp.epis,qtl.epis[k,])

              for(i in 1:nrow(tmp.epis)) 
                    {
                  xa1 <- qtl.geno[ ,tmp.epis[i,"qtl.a"]]-2
                  xd1 <- (1+xa1)*(1-xa1)-0.5
                  xa2 <- qtl.geno[ ,tmp.epis[i,"qtl.b"]]-2
                  xd2 <- (1+xa2)*(1-xa2)-0.5
    pheno.normal[,m] <- pheno.normal[,m]+xa1*xa2*tmp.epis[i,"aa"]+xa1*xd2*tmp.epis[i,"ad"]
                                      +xd1*xa2*tmp.epis[i,"da"]+xd1*xd2*tmp.epis[i,"dd"]  
                     }
                 }
         }
         if(!multiple.trait)
         {
              dimnames(qtl.epis) <- list(paste("epis", seq(n.epis), sep = "."),
                              c("qtl.a", "qtl.b", "aa", "ad", "da", "dd"))
              for (i in 1:nrow(qtl.epis)) {
                  xa1 <- qtl.geno[, qtl.epis[i, "qtl.a"]] - 2
                  xd1 <- (1 + xa1) * (1 - xa1) - 0.5
                  xa2 <- qtl.geno[, qtl.epis[i, "qtl.b"]] - 2
                  xd2 <- (1 + xa2) * (1 - xa2) - 0.5
                  pheno.normal <- pheno.normal + xa1*xa2*qtl.epis[i,"aa"]+ 
                  xa1*xd2*qtl.epis[i,"ad"] + xd1*xa2*qtl.epis[i,"da"] + 
                  xd1*xd2*qtl.epis[i,"dd"]
                }
 
          }  
       }  
    }
}
    if(!is.null(covariate)) {
       names(covariate) <- c("fix.cov", "ran.cov")
       fix.cov <- sample(c(0,1),n.ind,replace = TRUE)
       random.cov <- sample(0:4,n.ind,replace = TRUE)
       random.coe <- rnorm(5,0,covariate[2]^0.5)
       random <- rep(0,n.ind)
       for(i in 1:n.ind) {
           for(j in 1:5) 
              if(random.cov[i]==j-1) random[i] <- random.coe[j]
       }
       for(m in 1:ncol(pheno.normal))
       pheno.normal[,m] <- pheno.normal[,m]+fix.cov*covariate[1]+random
       if(!is.null(gbye) & n.qtl>0) {
          n.gbye <- nrow(gbye)
       if(type=="bc" | type=="riself" | type=="risib") {
           if(multiple.trait)
           {
              dimnames(gbye) <- list(paste("GxE", seq(n.gbye), sep = "."),
                                     c("pheno","qtl", "add"))
           for(m in unique(gbye[,"pheno"])) {
              tmp.gbye = rbind(gbye[gbye[,"pheno"]==m,])
              for(i in 1:nrow(tmp.gbye)) {
                xa <- qtl.geno[ ,tmp.gbye[i,"qtl"]]-1.5
                pheno.normal[,m] <- pheno.normal[,m]+xa*fix.cov*tmp.gbye[i,"add"]
               }
              }
          }     
          if(!multiple.trait)
          {  
            dimnames(gbye) <- list(paste("GxE", seq(n.gbye), sep = "."),
                                     c("qtl", "add"))
              for(i in 1:nrow(gbye)) {
                  xa <- qtl.geno[ ,gbye[i,"qtl"]]-1.5
                  pheno.normal <- pheno.normal+xa*fix.cov*gbye[i,"add"] 
              }  
           } 
        }
        if(type=="f2") {
           if(multiple.trait)
           {
             dimnames(gbye) <- list(paste("GxE", seq(n.gbye), sep = "."),
                                     c("pheno","qtl","add","dom"))
           for(m in unique(gbye[,"pheno"])) {
              tmp.gbye = rbind(gbye[gbye[,"pheno"]==m,])
              for(i in 1:nrow(tmp.gbye)) {
                xa <- qtl.geno[ ,tmp.gbye[i,"qtl"]]-2
                xd <- (1+xa)*(1-xa)-0.5
      pheno.normal[,m] <- pheno.normal[,m]+xa*fix.cov*gbye[i,"add"]+xd*fix.cov*gbye[i,"dom"]
               }
              }
          }     
          if(!multiple.trait)
          {  
             dimnames(gbye) <- list(paste("GxE", seq(n.gbye), sep = "."),
                                     c("qtl", "add", "dom"))
              for(i in 1:nrow(gbye)) {
                  xa <- qtl.geno[ ,gbye[i,"qtl"]]-2
                  xd <- (1+xa)*(1-xa)-0.5
             pheno.normal <- pheno.normal+xa*fix.cov*gbye[i,"add"]+xd*fix.cov*gbye[i,"dom"]
                  }   
           } 
       }
    }
    }
    
# generate genotypic values
 gvalue <- pheno.normal - residual

# truncate pheno.nomal to ordinal (binary) phenotypes
 if(!multiple.trait){
    p <- quantile(pheno.normal,probs=cumsum(ordinal))
    pheno.ordinal <- as.numeric(factor(cut(pheno.normal,breaks=c(-Inf,p,Inf)))) - 1
}
  
# save QTL positions, heritabilites of main and epistatic effects into list qtl
# save variances explained by covariates and g-e interactions into qtl
 
    qtl <- list()
    vp <- var(pheno.normal)
#    qtl$var <- sigma
    if(n.qtl > 0) {
       qtl$geno <- qtl.geno
       qtl$pos <- qtl.pos 
  if(!is.null(qtl.main)) {
      if(multiple.trait) v=vp[qtl.main[,"pheno"],qtl.main[,"pheno"]]
      else v=vp
      if(is.matrix(v)) v=diag(v) 
      if(type=="bc" | type=="riself" | type=="risib") 
       qtl.main[ ,"add"] <- qtl.main[ ,"add"]^2/(4*v)
      if(type=="f2") {
                      qtl.main[ ,"add"] <- qtl.main[ ,"add"]^2/(2*v)
                      qtl.main[ ,"dom"] <- qtl.main[ ,"dom"]^2/(4*v)
                       }
        qtl$herit.main <- qtl.main
        }     
       if(!is.null(qtl.epis)) {

          if(multiple.trait) v=vp[qtl.epis[,"pheno"],qtl.epis[,"pheno"]]
          else v=vp
          if(is.matrix(v)) v=diag(v)
          if(type=="bc") qtl.epis[ ,"aa"] <- qtl.epis[ ,"aa"]^2/(16*v)
          if(type=="f2") {
             qtl.epis[ ,"aa"] <- qtl.epis[ ,"aa"]^2/(4*v)
             qtl.epis[ ,"ad"] <- qtl.epis[ ,"ad"]^2/(8*v)
             qtl.epis[ ,"da"] <- qtl.epis[ ,"da"]^2/(8*v)
             qtl.epis[ ,"dd"] <- qtl.epis[ ,"dd"]^2/(16*v)
              }
          qtl$herit.epis <- qtl.epis
       }
    } 

    if(!is.null(covariate)) {
      if(multiple.trait)
      {
        colnames(vp)=paste("Pheno",seq(n.pheno),sep=".")
        rownames(vp)=paste("Pheno",seq(n.pheno),sep=".")
        v=diag(vp)
      } else v=vp  
       qtl$herit.cov <- cbind(0.25*covariate[1]^2/v,covariate[2]/v)
       colnames(qtl$herit.cov)=names(covariate)
       if(!is.null(gbye)) {
          if(multiple.trait)
          {
          v=vp[gbye[,"pheno"],gbye[,"pheno"]]
          if(is.matrix(v)) v=diag(v)
          } else v=vp
          if(type=="bc") gbye[ ,"add"] <- gbye[ ,"add"]^2*0.25/(4*v)
          if(type=="f2") {
             gbye[ ,"add"] <- gbye[ ,"add"]^2*0.25/(2*v)
             gbye[ ,"dom"] <- gbye[ ,"dom"]^2*0.25/(4*v)
          } 
          qtl$herit.gbye <- gbye
       }
    }
    if(multiple.trait){
        colnames(sigma)=paste("Pheno",seq(n.pheno),sep=".")
        rownames(sigma)=paste("Pheno",seq(n.pheno),sep=".")
     qtl$var <- sigma
     }
    class(qtl) <- c("qb.sim", "list")

# remove QTL positions and genotypes from cross$geno

    for (i in 1:n.chr) {
          o <- grep("^QTL[0-9]+", mar.names[[i]])  
         if (length(o) != 0) {
            geno[[i]]$data <- geno[[i]]$data[,-o,drop = FALSE]
            geno[[i]]$map <- geno[[i]]$map[-o]
        }
    }

# randomly give missing marker genotypes

    for (i in 1:n.chr) {
        z <- sample(c(1,NA), n.ind * ncol(geno[[i]]$data), replace = TRUE,
                    prob = c(1-missing.geno, missing.geno))
        dim(z) <- c(n.ind, ncol(geno[[i]]$data))
        geno[[i]]$data <- geno[[i]]$data*z
    }

# randomly give missing phenotypes and covariates
for(m in 1:n.pheno)
    pheno.normal[,m] = pheno.normal[,m]*sample(c(1,NA), n.ind, replace = TRUE,
                  prob = c(1-missing.pheno, missing.pheno))
  
if(!multiple.trait) pheno.ordinal = pheno.ordinal*sample(c(1,NA), n.ind, replace = TRUE,
                      prob = c(1-missing.pheno, missing.pheno))
    if(!is.null(covariate)) {
       fix.cov = fix.cov*sample(c(1,NA), n.ind, replace = TRUE,
         prob = c(1-missing.pheno, missing.pheno))
       random.cov = random.cov*sample(c(1,NA), n.ind, replace = TRUE,
         prob = c(1-missing.pheno, missing.pheno))
    }                       
           
# create a cross dataset
    if(!multiple.trait) pheno <- data.frame(pheno.normal=pheno.normal,pheno.ordinal=pheno.ordinal)
    if(multiple.trait)  pheno <-data.frame(pheno.normal=pheno.normal)
  
    if(!is.null(covariate))
       pheno <- data.frame(fix.cov=fix.cov,random.cov=random.cov,pheno)
  
    cross <- list(geno=geno, pheno=pheno, qtl=qtl, gvalue=gvalue)
    class(cross) <- c(type, "cross")
    for (i in 1:n.chr) storage.mode(cross$geno[[i]]$data) <- "integer"

    cross
}


summary.qb.sim <- function(object, ...) object[-1]








    
