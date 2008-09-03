#####################################################################
##
## $Id: pleiotropy.R,v 1.22 2007/10/3 19:26:43 samban Exp $
##
##     Copyright (C) 2007 Samprit Banerjee
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
##
##############################################################################
qb.inter.pleiotropy <- function(qbObject,n.pheno)
{
 tmp = vector(mode="list",length=n.pheno)
 for(i in 1:n.pheno)
 {
  mainloci <- qb.get(qbObject, "mainloci",pheno=i)
tmp[[i]]=paste(mainloci[, "niter"],paste(mainloci[, "chrom"], mainloci[, "locus"], sep = ":"),sep="->")
 }

 
  tt = match(tmp[[1]],tmp[[2]],nomatch=FALSE)
  if(n.pheno > 2)
   {   for(i in 3:n.pheno)
        tt = match(tmp[[1]][which(tt!=0)],tmp[[i]],nomatch=FALSE)
    }
  H1=tmp[[1]][which(tt!=0)]
index.H1 = as.integer(unlist(lapply(H1,function(x) unlist(strsplit(x,split="->"))[1])))

  H1 = unlist(lapply(H1,function(x) unlist(strsplit(x,split="->"))[2]))

  inter = vector(mode="list",length=2)
  grid = pull.grid(qbObject, offset = TRUE)  
  if(is.null(H1)) H1 <- NA
  inter[[1]] = ordered(H1,paste(grid[, 1], grid[, 2], sep = ":"))

  inter[[2]] = index.H1
  names(inter) <- c("inter","index")
 inter
}

################################################################################
################ function to test for pleiotropic effects ######################
################################################################################
qb.pleiotropy <- function(qbObject,scan = "main",
                                   type = types)
{

  qb.name <- deparse(substitute(qbObject))
  cross <- qb.cross(qbObject) ## YANDELL.
  pheno.name <- names(cross$pheno)[qb.get(qbObject, "pheno.col")]
  geno.names <- names(cross$geno)

  types <- c("posterior","2logBF","BF")
  
  type <- types[pmatch(tolower(type), tolower(types), nomatch = 2)][1]

  is.bc <- (qb.get(qbObject, "cross") == "bc")

  if(scan == "main")
    { vars <- "add"
      element <- paste("var", vars, sep = "")
    }  
  else stop("qb.pleiotropy can handle main effects only",call.=TRUE)

    n.iter <- qb.get(qbObject,"n.iter")
    n.pheno = qb.get(qbObject,"n.pheno")


 main.val=vector(mode="list",length=n.pheno)
    for(i in 1:n.pheno)
    {
      mainloci = qb.get(qbObject,"mainloci",pheno=i)
      main.val[[i]] = mainloci[,element]
     } 

  inter.H1 = qb.inter.pleiotropy(qbObject,n.pheno)
  condition = (main.val[[1]][inter.H1$index] > 0)
    for(i in 2:n.pheno)
      condition = (condition & (main.val[[i]][inter.H1$index] > 0))

 if(!any(is.na(inter.H1$inter)))   temp.H1 <- unlist(tapply(condition, inter.H1$inter, sum, na.rm = TRUE)) else  temp.H1 <- NA
 
  temp.H1[is.na(temp.H1)] <- 0
  stat.H1 = (temp.H1)/n.iter

  x <- matrix(0, length(levels(inter.H1[[1]])), length(scan))
  dimnames(x) <- list(NULL, vars)

 if(type!="posterior")
 { 
  bf.prior <- qb.get(qbObject, "mean.nqtl") / length(unlist(pull.loci(cross)))
  stat <- (stat.H1*(1-bf.prior^n.pheno))/((1-stat.H1)*bf.prior^n.pheno)
     if(type == "2logBF") stat <- 2 * log(pmax(1, stat))
          
 } else stat <- stat.H1                

     x[,"add"] <- stat

  x <- cbind(x, count = temp.H1)
 
  attr(x, "class") <- c("qb.pleiotropy", "matrix")
  attr(x, "type.scan") <- type
  attr(x, "scan") <- scan.save <- vars
  attr(x, "cross.class") <- qb.cross.class(qbObject)
  attr(x, "pheno.name") <- pheno.name
  attr(x, "geno.names") <- geno.names
  attr(x, "qb") <- qb.name

  x
}
################################################################################
################### summary routine for qb.pleiotropy ##########################
################################################################################
summary.qb.pleiotropy <- function(object, threshold = NULL)
{
 grid <- attr(object, "grid")
 if(is.null(threshold))
   threshold <- ifelse(attr(object, "type.scan") == "posterior", 0.001, 1)

 out <- cbind(grid, object)
 out <- out[out[,attr(object, "scan")] > threshold,]
 cat(paste("The", attr(object, "type.scan"), "profile for testing pleiotropic effects"))
 cat("\n\n")
 out
}
################################################################################
print.qb.pleiotropy <- function(x, ...) print(summary(x))
################################################################################

################################################################################
################### plotting routine for qb.pleiotropy #########################
################################################################################
plot.qb.pleiotropy <- function(x,chr=NULL,
                               main = paste(type, "for pleiotropic test for"
                              ,scan.plots),
                            ...)
{
  attr(x, "min.iter") <- min.iter <- 1
  type <- attr(x,"type.scan")
  class(x)[1] <- "qb.scanone"
  scan.plots = attr(x,"scan")
  plot.qb.scanone(x,chr=chr,scan = scan.plots,main = main,
                  weight = x[, "count"], ...)

}                            
