#####################################################################
##
## $Id: scan.R,v 1.11.2.10 2006/12/12 19:23:28 byandell Exp $
##
##     Copyright (C) 2007 Brian S. Yandell
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
## Need to watch out for extra digits in locus when matching up with grid.
join.chr.pos <- function(chrom, locus, digits = 10)
  paste(chrom, signif(locus, digits), sep = ":")
##  paste(chrom, locus, sep = ":")
make.chr.pos <- function(chrom, locus,
                         level.chrom = chrom, level.locus = locus,
                         levels = unique(join.chr.pos(level.chrom,
                           level.locus, ...)),
                         ...)
  ordered(join.chr.pos(chrom, locus, ...), levels)
##############################################################################
qb.inter <- function(qbObject, x = pull.grid(qbObject, offset = TRUE),
                     mainloci = qb.get(qbObject, "mainloci", ...), ...)
{
  ## Create identifier of chrom.locus from mainloci into pseudomarker grid.
  inter <- make.chr.pos(mainloci[, "chrom"], mainloci[, "locus"], x[, 1], x[, 2])
  tmp <- is.na(inter)
  if(any(tmp)) {
    stop(paste("qb.scanone or qb.sliceone mismatch with grid:\n", sum(tmp),
               "missing values generated\n"))
  }
  inter
}
##############################################################################
qb.threshold <- function(out, threshold, pos = 1)
{
  ## Internal routine.
  
  ## Keep only chrs with some value about threshold.
  if(any(threshold != 0)) {
    if(is.null(names(threshold)))
      names(threshold) <- dimnames(out)[[2]][pos + seq(length(threshold))]
    threshold <- threshold[!is.na(match(names(threshold), dimnames(out)[[2]]))]
    
    if(is.null(threshold) | !length(threshold))
      return(out)
    
    use <- threshold >= 0
    if(any(use)) {
      keep <- apply(out[, names(threshold[use]), drop = FALSE], 1,
                    function(x) any(x >= threshold[use]))
    }
    else
      keep <- rep(FALSE, nrow(out))
    use <- threshold <= 0
    if(any(use)) {
      maxout <- apply(out[, names(threshold[use]), drop = FALSE],
                      2, max)
      keep <- keep | apply(out[, names(threshold[use]), drop = FALSE],
                           1,
                           function(x)
                           any(x >= maxout + threshold[use]))
    }
    out[keep, , drop = FALSE]
  }
  else
    out
}
##############################################################################
qb.count <- function(stat, type.scan, n.iter, bf.prior)
{
  if(type.scan == "log10")
    stat <- log10(1 + stat)
  else if(type.scan != "count" & type.scan != "nqtl") {
    ## Else type.scan is posterior or BF.
    stat <- (1 + stat) / (2 + n.iter)
    if(type.scan == "logposterior") {
      stat <- log10(stat)
    }
    else {
      if(match(type.scan, c("2logBF","BF"), nomatch = 0)) {
        stat <- stat * (1 - bf.prior) / ((1 - stat) * bf.prior)
        if(type.scan == "2logBF")
          stat <- 2 * log(pmax(1, stat))
      }
    }
  }
  stat
}
##############################################################################
qb.nind.pheno <- function(qbObject,
                          pheno.name = pheno.names[qb.get(qbObject, "pheno.col")[1]],
                          nfixcov, cross,
                          covar.name = pheno.names[qb.get(qbObject, "covar")])
{
  pheno.names <- qb.pheno.names(qbObject, cross)
  
  not.missing <- apply(cross$pheno[, pheno.name, drop = FALSE], 1,
               function(x) all(!is.na(x) & abs(x) != Inf))
  if(nfixcov) {
    ## Reduce pheno count by missing covariate values.
    not.missing <- not.missing & {
      apply(cross$pheno[, covar.name, drop = FALSE], 1,
            function(x) all(!is.na(x) & abs(x) != Inf))
    }
  }
  sum(not.missing)
}
##############################################################################
qb.npar <- function(var1, var2, nfixcov, nrancov, intcov, iterdiag.nqtl,
                    iterdiag, mainloci, gbye, pairloci)
{
  ## Number of parameters in QTL model averaged over MCMC runs.
  ## Not correct for covariates yet!
  npar <- rep(0, nrow(iterdiag))
  
  ## Main QTL degrees of freedom.
  tmp <- apply(as.matrix(mainloci[, paste("var", var1, sep = "")]), 1,
               function(x) sum(x > 0))
  tmp <- tapply(tmp, mainloci[, "niter"], sum)
  npar[iterdiag.nqtl > 0] <- tmp
  
  ## Covariate degrees of freedom.
  npar <- npar + nfixcov + nrancov

  ## Get GxE degrees of freedom.
  if(sum(intcov)) {
    if(length(gbye)) {
      same <- match(paste(gbye[, "niter"], gbye[, "chrom"],
                          gbye[, "locus"], sep = ":"),
                    paste(mainloci[, "niter"], mainloci[, "chrom"],
                          mainloci[, "locus"], sep = ":"))
      
      tmp <- rep(0, nrow(mainloci))
      tmp[same] <- apply(as.matrix(gbye[, paste("var", var1, sep = "")]), 1,
                         function(x) sum(x > 0))
      tmp <- tapply(tmp, mainloci[, "niter"], sum)
      npar[iterdiag.nqtl > 0] <- npar[iterdiag.nqtl > 0] + tmp
    }
  }

  ## Epistasis degrees of freedom.
  if(length(pairloci)) {
    tmp <- apply(as.matrix(pairloci[, paste("var", var2, sep = "")]), 1, function(x) sum(x > 0))
    tmp <- tapply(tmp, pairloci[, "niter"], sum)
    tmp2 <- match(unique(pairloci[, "niter"]), iterdiag[, "niter"], nomatch = 0)
    npar[tmp2] <- npar[tmp2] + tmp
  }
  npar
}
##############################################################################
qb.reference <- function(qbObject, mainloci, iterdiag, inter, type.scan)
{
  if(any(type.scan == "cellmean")) {
    reference <- unlist(tapply(qb.meancomp(qbObject)[match(mainloci[, "niter"],
                                                     iterdiag[, "niter"]),
                                               "grand.mean"],
                         inter, mean))
    reference[is.na(reference)] <- mean(reference, na.rm = TRUE)
    is.bc <- (qb.cross.class(qbObject) == "bc")
    attr(reference, "genos") <- c("A","H", if (!is.bc) "B")
    reference
  }
  else
    0
}
##############################################################################
qb.scanmain <- function(x, type.scan, is.bc, scans, scan.save, n.iter, bf.prior,
                        reference, covar, covar.means, inter, mainloci,
                        gbye, intcov, nfixcov, sum.scan, qb.coef)
{
  totvar <- 0
  is.count <- type.scan %in% c("count", "log10", "posterior", "logposterior",
                          "2logBF", "BF", "nqtl")

  ## Get non-epistatic components: additive and dominance.
  vars <- c("add", if(!is.bc) "dom")
  if(type.scan != "heritability")
    vars <- vars[match(scans, vars, nomatch = 0)]
  
  if(length(vars)) {
    ## Number of main effect samples per locus.
    if(is.count & any(scan.save == "main")) {
      tmp <- apply(as.matrix(mainloci[, paste("var", vars, sep = "")]), 1,
                   function(x) any(x > 0))
      if(type.scan == "nqtl") {
        nqtl.main <- paste(mainloci[, "niter"], mainloci[, "chrom"], sep = ":")
        tmp <- tapply(tmp, nqtl.main, sum)[nqtl.main]
        tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
      }
      else
        tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
      tmp[is.na(tmp)] <- 0
      x[, "main"] <- qb.count(tmp, type.scan, n.iter, bf.prior)
    }
    else if(type.scan == "cellmean") {
      for(i in attr(reference, "genos"))
        x[, i] <- reference
    }
    for(i in vars) {
      if(type.scan == "estimate" | type.scan == "cellmean")
        ## Parameter estimates of main effects.
        element <- i
      else
        ## Variance components.
        element <- paste("var", i, sep = "")
      
      ## Get samples for this component.
      main.val <- mainloci[, element]

      ## Get GxE samples if any intcov selected.
      if(sum(intcov)) {
        ## Get GxE samples.
        if(any(scan.save == "GxE"))
          tmp2 <- rep(0, nrow(mainloci))
        
        ## Loop over all covariates.
        covars <- seq(nfixcov)[intcov]
        for(j in covars) {
          gbyej <- gbye[gbye[, "covar"] == j, ]
          if(length(gbyej)) {
            same <- match(paste(gbyej[, "niter"], gbyej[, "chrom"],
                                gbyej[, "locus"], sep = ":"),
                          paste(mainloci[, "niter"], mainloci[, "chrom"],
                                mainloci[, "locus"], sep = ":"))
            
            ## Parameter estimates of GxE fixed effects.
            if(match(j, covar, nomatch = 0) | type.scan == "heritability") {
              cname <- paste(paste(i, names(covar.means)[j], sep = "."))
              tmp <- rep(0, length(main.val))
              tmp[same] <- gbyej[[element]]
              ## Include the covariate?
              if(is.count) {
                ## Count times GxE effect is present.
                tmp <- tmp > 0
                if(any(scan.save == "GxE"))
                  tmp2 <- pmax(tmp2, tmp)
                if(any(scan.save == cname)) {
                  if(type.scan == "nqtl") {
                    tmp <- tapply(tmp, nqtl.main, sum)[nqtl.main]
                    tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
                  }
                  else
                    tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
                  tmp[is.na(tmp)] <- 0
                  x[, cname] <- qb.count(tmp, type.scan, n.iter, bf.prior)
                }
              }
              else { ## is.effect
                tmp <- unlist(tapply(tmp, inter, mean))
                tmp[is.na(tmp)] <- 0
                if(type.scan == "heritability")
                  totvar <- totvar + tmp
                if(match(cname, scan.save, nomatch = 0))
                  x[,cname] <- tmp
                if(any(scan.save == "GxE"))
                  x[,"GxE"] <- x[,"GxE"] + tmp
                if(sum.scan != "no")
                  x[,"sum"] <- x[,"sum"] + tmp
              }
            }
            
            if(type.scan == "estimate" | type.scan == "cellmean") {
              ## Offset parameter estimate by covariates.
              if(covar.means[j] != 0)
                main.val[same] <- main.val[same] + covar.means[j] * gbyej[[i]]
            }
          }
        }
        if(is.count & any(scan.save == "GxE")) {
          if(type.scan == "nqtl") {
            tmp2 <- tapply(tmp2, nqtl.main, sum)[nqtl.main]
            tmp2 <- unlist(tapply(tmp2, inter, mean, na.rm = TRUE))
          }
          else
            tmp2 <- unlist(tapply(tmp2, inter, sum, na.rm = TRUE))
          tmp2[is.na(tmp2)] <- 0
          x[, "GxE"] <- qb.count(tmp2, type.scan, n.iter, bf.prior)
        }
      }

      ## Now include the main effect components.
      if(is.count) {
        if(any(scan.save == i)) {
          ## Count times main effect is present.
          if(type.scan == "nqtl") {
            tmp <- tapply(main.val > 0, nqtl.main, sum)[nqtl.main]
            tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
          }
          else
            tmp <- unlist(tapply(main.val > 0, inter, sum, na.rm = TRUE))
          tmp[is.na(tmp)] <- 0
          x[, i] <- qb.count(tmp, type.scan, n.iter, bf.prior)
        }
      }
      else { ## is.effect
        ## Get main effect element and average.
        tmp <- unlist(tapply(main.val, inter, mean, na.rm = TRUE))
        tmp[is.na(tmp)] <- 0

        if(type.scan == "cellmean") {
          ## Add contribution of element to cell mean.
          if(any(names(qb.coef) == i))
            x <- x + outer(tmp, qb.coef[[i]])
        }
        else {
          if(type.scan == "heritability")
            totvar <- totvar + tmp
          if(match(i, scan.save, nomatch = 0))
            x[,i] <- tmp
          if(any(scan.save == "main"))
            x[,"main"] <- x[,"main"] + tmp
          if(sum.scan != "no")
            x[,"sum"] <- x[,"sum"] + tmp
        }
      }
    }
  }
  list(totvar = totvar, x = x)
}
##############################################################################
qb.scanepis <- function(x, type.scan, is.bc, scans, scan.save, n.iter, bf.prior,
                        inter, mainloci, pairloci, sum.scan, qb.coef,
                        half = FALSE, is.slice = FALSE)
{
  totvar <- 0
  is.count <- type.scan %in% c("count", "log10", "posterior", "logposterior",
                          "2logBF", "BF", "nqtl")

  ## Index for epistasis.
  nqtl.pair <- c(paste(pairloci[, "niter"], pairloci[, "chrom1"], sep = ":"),
                 paste(pairloci[, "niter"], pairloci[, "chrom2"], sep = ":"))
  epinter <- make.chr.pos(nqtl.pair,
                          c(pairloci[, "locus1"], pairloci[, "locus2"]))

  ## epii identifies mainloci with epistatic pairs.
  epii <- match(unique(epinter),
                paste(mainloci[, "niter"], inter, sep = ":"),
                nomatch = 0)

  ## Epistatic components.
  vars <- c("aa", if(!is.bc) c("ad","da","dd"))
  if(type.scan != "heritability")
    vars <- vars[match(scans, vars, nomatch = 0)]
  if(length(vars)) {
    tmp <- rep(0, length(inter))
    
    ## Number of epistatic samples per locus.
    if(is.count & any(scan.save == "epistasis")) {
      if(type.scan == "nqtl") {
        tmp2 <- table(nqtl.pair)[nqtl.pair]
        tmp[epii] <- unlist(tapply(tmp2, epinter, mean))[epii > 0]
        tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
      }
      else {
        tmp[epii] <- 1
        tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
      }
      tmp[is.na(tmp)] <- 0
      x[, "epistasis"] <- qb.count(tmp, type.scan, n.iter, bf.prior)
    }
    
    for(i in vars) {
      if(type.scan == "estimate" | (type.scan == "cellmean" & is.slice))
        ## Parameter estimates of epistasis.
        element <- i
      else
        ## Variance components for epistasis.
        element <- paste("var", i, sep = "")
      
      tmp2 <- pairloci[, element]
      if(is.count) {
        if(any(scan.save == i)) {
          ## Count times epistatic element is present.
          if(type.scan == "nqtl") {
            ## NB: This multiply counts locus used in several pairs.
            tmp2 <- unlist(tapply(rep(tmp2 > 0, 2),
                                  nqtl.pair, sum))[nqtl.pair]
            tmp[epii] <- unlist(tapply(tmp2, epinter, mean))[epii > 0]
            tmp2 <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
          }
          else {
            tmp[epii] <- unlist(tapply(rep(tmp2 > 0, 2), epinter,
                                       sum))[epii > 0]
            tmp2 <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
          }
          tmp2[is.na(tmp2)] <- 0
          ## Compute count diagnostic.
          x[, i] <- qb.count(tmp2, type.scan, n.iter, bf.prior)
          
        }
      }
      else { ## is.effect
        ## Get epistatic element and average.
        tmp[epii] <- unlist(tapply(rep(tmp2, 2), epinter, sum))[epii > 0]
        tmp2 <- unlist(tapply(tmp, inter, mean))
        
        if(type.scan == "cellmean" & is.slice) {
          ## Add contribution of element to cell mean.
          if(any(names(qb.coef) == i))
            x <- x + outer(tmp2, qb.coef[[i]])
        }
        else {
          ## Reduce epistatic value by half.
          if(half)
            tmp2 <- tmp2 / 2
          tmp2[is.na(tmp2)] <- 0
          if(type.scan == "heritability")
            totvar <- totvar + tmp2
          if(match(i, scan.save, nomatch = 0))
            x[,i] <- tmp2
          if(any(scan.save == "epistasis"))
            x[,"epistasis"] <- x[,"epistasis"] + tmp2
          if(sum.scan != "no")
            x[,"sum"] <- x[,"sum"] + tmp2
        }
      }
    }
  }
  list(totvar = totvar, x = x)
}
##############################################################################
qb.scanone <- function(qbObject, epistasis = TRUE,
                       scan = c("main", "GxE", "epistasis"),
                       type.scan = type.scans,
                       covar = if(nfixcov) seq(nfixcov) else 0,
                       adjust.covar = NA,
                       chr = NULL,
                       sum.scan = "yes",
                       min.iter = 1,
                       aggregate = TRUE,
                       smooth = 3,
                       weight = c("sqrt","count","none","atten","ratten"),
                       split.chr = qb.get(qbObject, "split.chr"),
                       center.type = c("mode","mean","scan"),
                       half = FALSE,
                       verbose = FALSE,
                       ...)
{
  type.scans <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","logposterior","2logBF","BF","nqtl",
             "npar","rss")
  nfixcov <- qb.get(qbObject, "nfixcov")

  qb.commonone(qbObject, "scanone",, epistasis, scan, type.scan, covar,
               adjust.covar, chr, sum.scan, min.iter, aggregate,
               smooth, weight, split.chr, center.type, half, verbose, ...)
}
###################################################################
check.intcov <- function(intcov, nfixcov)
{
  tmp <- length(intcov)
  if(tmp > nfixcov) {
    if(nfixcov == 0)
      intcov <- NULL
    else
      intcov <- intcov[seq(nfixcov)]
  }
  else if(tmp < nfixcov)
    intcov <- c(intcov, rep(FALSE, nfixcov - tmp))

  ## Check that intcov is of length nfixcov.
  if(length(intcov) != nfixcov)
    stop(paste("mismatch in qb object: intcov length (", sum(intcov),
               ") via qb.model must match fixcov length(",
               nfixcov, ") via qb.data", sep = ""))
  intcov
}
###################################################################
qb.commonone <- function(qbObject,
                         call = "scanone",
                         slice,
                         epistasis = TRUE,
                         scan = c("main", "GxE", "epistasis"),
                         type.scan = type.scans,
                         covar = if(nfixcov) seq(nfixcov) else 0,
                         adjust.covar = NA,
                         chr = NULL,
                         sum.scan = "yes",
                         min.iter = 1,
                         aggregate = TRUE,
                         smooth = 3,
                         weight = c("sqrt","count","none","atten","ratten"),
                         split.chr = qb.get(qbObject, "split.chr"),
                         center.type = c("mode","mean","scan"),
                         half = FALSE,
                         verbose = FALSE,
                         pheno.col = qb.get(qbObject, "pheno.col")[1], ...)
{
  qb.exists(qbObject)
  
  qb.name <- deparse(substitute(qbObject))
  
  is.bc <- (qb.cross.class(qbObject) == "bc")


  ## Force chr and slice["chr"] to index chromosomes used in qbObject.
  chr <- qb.find.chr(qbObject, chr)

  is.slice <- (call == "sliceone")
  if(is.slice) {
    ## Restrict attention to selected chromosomes
    if(missing(slice))
      stop("must specify chromosome to slice upon")
    
    ## Set up slice vector.
    ## slice = c(chr=, upper=TRUE, start=, end=, weight=c(0,1,2))
    ## NB: slice could be list, with first element logical or character.
    if(is.null(names(slice)))
      names(slice) <- c("chr","start","end")[seq(length(slice))]
    slice["chr"] <- qb.find.chr(qbObject, unlist(slice["chr"]))
    slice <- unlist(slice)
    if(is.na(slice["start"]))
      slice["start"] <- 0
    if(is.na(slice["end"]))
      slice["end"] <- max(pull.grid(qbObject)$pos)
    tmp <- names(slice)
    slice <- as.numeric(slice)
    names(slice) <- tmp
  }

  ## Determine variance components.
  var1 <- "add"
  var2 <- ifelse(epistasis, "aa", character(0))
  if(!is.bc) {
    var1 <- c(var1,"dom")
    if(epistasis)
      var2 <- c(var2,"ad","da","dd")
  }

  ## Determine coefficients used in cellmean and ideal heritability.
  if(is.slice) {
    qb.coef <- if(is.bc)
      list(add = rep(c(-0.5,0.5), 2),
           aa = c(0.25,-0.25,-0.25,0.25))
    else
      list(add = rep(c(-1,0,1), 3), dom = rep(c(-0.5,0.5,-0.5), 3),
           aa = c(1,0,-1,0,0,0,-1,0,1),
           ad = c(0.5,0,-0.5,-0.5,0,0.5,0.5,0,-0.5),
           da = c(0.5,-0.5,0.5,0,0,0,-0.5,0.5,-0.5),
           dd = c(0.25,-0.25,0.25,-0.25,0.25,-0.25,0.25,-0.25,0.25))
  }
  else {
    qb.coef <- if(is.bc)
      list(add = c(-0.5,0.5))
    else
      list(add = c(-1,0,1), dom = c(-0.5,0.5,-0.5))
  }
  
  ## Number of fixed covariates
  nfixcov <- qb.get(qbObject, "nfixcov")
  nrancov <- qb.get(qbObject, "nrancov")
  intcov <- as.logical(qb.get(qbObject, "intcov"))
  intcov <- check.intcov(intcov, nfixcov)
  
  ## Determine type of scan.
  type.scans <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","logposterior","2logBF","BF","nqtl",
             "npar","rss")
  type.scan <- type.scans[pmatch(tolower(type.scan), tolower(type.scans), nomatch = 2)][1]

  is.count <- type.scan %in% c("count", "log10", "posterior", "logposterior",
                          "2logBF", "BF", "nqtl")
  is.var <- type.scan %in% type.scans[1:6]
  is.effect <- is.var | type.scan == "estimate"
  is.lod <- type.scan %in% type.scans[c(2:5,16:17)]

  ## Number of individuals for phenotype.
  cross <- qb.cross(qbObject, genoprob = FALSE)
  pheno.name <- qb.pheno.names(qbObject, cross)[pheno.col][1]
  nind.pheno <- qb.nind.pheno(qbObject, pheno.name, nfixcov, cross)

  ## Genotype names.
  geno.names <- qb.geno.names(qbObject, cross)

  ## Following prior used for Bayes factors.
  bf.prior <- qb.get(qbObject, "mean.nqtl") / qb.nloci(qbObject, cross)
  if(is.slice)
    bf.prior <- bf.prior * bf.prior

  rm(cross)
  gc()

  ## Subset on chromosomes.
  qbObject <- subset(qbObject,
                     chr = {
                       if(is.slice)
                         sort(unique(c(chr, slice["chr"])))
                       else
                         chr},
                     restrict.pair = FALSE)

  ## Pull grid of loci.
  grid <- pull.grid(qbObject, offset = TRUE)

  if(is.slice) {
    ## Restrict attention to samples including slice.
    qbObject <- subset(qbObject, region = list(chr=slice["chr"],
                                   start=slice["start"], end=slice["end"]),
                       restrict.pair = FALSE)
  }
  
  ## Get MCMC samples.
  iterdiag <- qb.get(qbObject, "iterdiag", pheno.col = pheno.name)
  mainloci <- qb.get(qbObject, "mainloci", pheno.col = pheno.name)

  ## No QTL at all.
  if(is.null(mainloci))
    return(NULL)
    
  iterdiag.nqtl <- qb.nqtl(qbObject, iterdiag, mainloci)

  ## Need pairloci earlier for slice?
  pairloci <- qb.get(qbObject, "pairloci", pheno.col = pheno.name)
  if(is.null(pairloci))
    epistasis <- FALSE
  else if(is.slice) {
    ## Restrict pairloci to pairs with slice.
    tmp <- pairloci$chrom1 == slice["chr"] | pairloci$chrom2 == slice["chr"]
    if(sum(tmp))
      pairloci <- pairloci[tmp,]
    else
      epistasis <- FALSE
  }
  
  ## Find interaction pattern.
  if(verbose)
    cat("finding loci ...\n")
  inter <- qb.inter(qbObject, grid, mainloci)

  ## Covariate adjustment calculations.
  if(type.scan == "heritability")
    totvar <- rep(0, length(levels(inter)))
  if(nfixcov) {
    ## Covariate means.
    covar.means <- covar.mean(qbObject, adjust.covar,
                              verbose = verbose & (type.scan == "estimate"),
                              pheno.col = pheno.col)
    
    ## Explained covariance for heritability.
    if(type.scan == "heritability") {
      tmp <- apply(qb.varcomp(qbObject, c("fixcov","rancov")), 1, sum)
      tmp <- unlist(tapply(tmp[match(mainloci[, "niter"], iterdiag[, "niter"])],
                           inter, mean))
      tmp[is.na(tmp)] <- 0
      totvar <- tmp
    }
  }

  ## Set up scan names.
  ## Scan can be several variance components (default is all).
  ## scan = original or default scan names
  ##   ("main","epistasis","GxE" or "A","H","B").
  ## scans = all terms needed for analysis (var1, var2, and var1.covar=GxE).
  ## scan.save = scan names to save in returned object.
  ## sum.scan = indicator whether "sum" name is returned ("yes", "no", "only").
  if(type.scan == "cellmean") {
    scan <- c("A","H","B")[seq(3 - is.bc)]
    if(is.slice) {
      scans <- c(var1, var2)
      scan <- c(outer(scan, scan, paste, sep = ""))
    }
    else {
      scans <- var1
    }
    scan.save <- scan
    sum.scan <- "no"
  }
  else {
    if(type.scan == "estimate" & sum.scan == "yes")
      sum.scan <- "no"
    aggregs <- c("main","epistasis","GxE","gbye")
    tmp <- pmatch(tolower(scan), tolower(aggregs), nomatch = 0)
    if(any(tmp))
      scan[tmp > 0] <- aggregs[pmin(tmp, 3)]
    aggregs <- aggregs[1:3]
    if(!(any(pmatch(scan, var1, nomatch = 0)) |
         any(match(scan, var2, nomatch = 0))) &
       any(match(scan, aggregs, nomatch = 0))) {
      if(!epistasis)
        scan <- scan[scan != aggregs[2]]
      if(!sum(intcov))
        scan <- scan[scan != aggregs[3]]
      scans <- scans.save <- NULL
      ## Include main if main or GxE requested.
      if(any(scan == aggregs[1]))
        scans <- scans.save <- var1
      else
        if(any(scan == aggregs[3]) & sum(intcov))
          scans <- var1
      if(any(scan == aggregs[2])) {
        scans <- c(scans, var2)
        scans.save <- c(scans.save, var2)
      }
      if(any(scan == aggregs[3])) {
        if(sum(intcov)) {
          if(length(covar.means)) {
            tmp <- seq(nfixcov)[intcov]
            tmp <- names(covar.means)[covar[match(tmp, covar, nomatch = 0)]]
            n.var1.covar <- length(tmp)
            if(n.var1.covar) {
              scans <- c(scans, outer(var1, tmp, paste, sep = "."))
              scans.save <- c(scans.save, outer(var1, tmp, paste, sep = "."))
            }
          }
        }
        else
          n.var1.covar <- 0
      }
    }
    else {
      scans <- scans.save <- scan
      ## Be sure main present if GxE terms are included in scan.
      if(sum(intcov)) {
        if(length(covar.means)) {
          tmp <- seq(nfixcov)[intcov]
          tmp <- names(covar.means)[covar[match(tmp, covar, nomatch = 0)]]
          if(length(tmp)) {
            for(vari in rev(var1)) {
              if(any(match(scans, outer(vari, tmp, paste, sep = "."),
                           nomatch = 0)))
                if(!any(match(scans, vari, nomatch = 0)))
                  scans <- c(vari, scans)
            }
          }
        }
      }
      ## Include GxE covariates in summary scans even if not saved.
      if(sum.scan != "no" & is.effect) {
        vars <- var1[match(scans, var1, nomatch = 0)]
        if(sum(intcov)) {
          if(length(covar.means) & length(vars)) {
            tmp <- seq(nfixcov)[intcov]
            tmp <- names(covar.means)[covar[match(tmp, covar, nomatch = 0)]]
            if(length(tmp))
              scans <- unique(c(scans, outer(vars, tmp, paste, sep = ".")))
          }
        }
      }
    }
    if(sum.scan %in% c("yes","only") & (length(scans) == 1))
      sum.scan <- "no"
    ## The vector scans contains elements to scan,
    ## either to show directly or to combine in sum.
    ## The vector scan.save indicates which elements to save.
    scan.save <- if(aggregate) scan else scans.save
    scan.save <- switch(sum.scan,
                        no = scan.save,
                        yes = c(scan.save, "sum"),
                        two =, only = "sum")
  }

  if(verbose) {
    cat("\n", type.scan, "of", pheno.name, "for",
        paste(scan.save, collapse = ","), "\n")
    if(min.iter > 1)
      cat("Including only loci pairs with at least", min.iter, "samples.\n")
    cat("\n")
  }

  ## NOW GET SAMPLES AVERAGED AT CHR and POS.
  x <- matrix(0, length(levels(inter)), length(scan.save))
  dimnames(x) <- list(NULL, scan.save)

  ## Extract environmental variance.
  if(type.scan == "heritability" | is.lod |
     (type.scan == "variance" & any(scan.save == "env"))) {
    if(verbose)
      cat("environmental variance ...\n")
    tmp <- unlist(tapply(iterdiag[match(mainloci[, "niter"], iterdiag[, "niter"]),
                                  "envvar"],
                         inter, mean))
    tmp[is.na(tmp)] <- 0
    if(type.scan == "heritability")
      totvar <- totvar + tmp
    else {
      if(is.lod)
        env <- tmp
      else if(any(scan.save == "env"))
        ## Actually want to examine env variance directly.
        x[,"env"] <- tmp
    }
  }

  ## Need n.iter for counts.
  if(is.count)
    n.iter <- qb.niter(qbObject)

  if(verbose)
    cat("non-epistatic components ...\n")
  ## Get non-epistatic components: additive and dominance.
  gbye <- if(sum(intcov))
    qb.get(qbObject, "gbye", pheno.col = pheno.name)
  else
    NULL
  reference <- qb.reference(qbObject, mainloci, iterdiag, inter, type.scan)

  if(is.slice) {
    ## Adjustment for slice: genos considers pairs.
    tmp2 <- attr(reference, "genos")
    attr(reference, "genos") <- c(outer(tmp2, tmp2, paste, sep = ""))
  }
    
  tmp <- qb.scanmain(x, type.scan, is.bc, scans, scan.save, n.iter, bf.prior,
                     reference, covar, covar.means, inter, mainloci,
                     gbye, intcov, nfixcov, sum.scan, qb.coef)
  x <- tmp$x
  if(type.scan == "heritability")
    totvar <- totvar + tmp$totvar
  reference <- mean(reference)

  if(epistasis) {
    if(verbose)
      cat("epistatic components ...\n")
    tmp <- qb.scanepis(x, type.scan, is.bc, scans, scan.save, n.iter, bf.prior,
                       inter, mainloci, pairloci, sum.scan, qb.coef, half, is.slice)
    x <- tmp$x
    if(type.scan == "heritability")
      totvar <- totvar + tmp$totvar
  }

  ## Extract counts averaged over MCMC runs.
  if(is.count) {
    if(verbose)
      cat("counts ...\n")
    if(sum.scan != "no") {
      ## Number of iterations per locus.
      if(type.scan == "nqtl") {
        nqtl.main <- paste(mainloci[, "niter"], mainloci[, "chrom"], sep = ":")
        tmp <- table(nqtl.main)[nqtl.main]
        tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
        tmp[is.na(tmp)] <- 0
      }
      else
        tmp <- unclass(table(inter))
      x[, "sum"] <- qb.count(tmp, type.scan, n.iter, bf.prior)
    }
    if(any(scans == "nqtl")) {
      ## Number of QTL averaged over MCMC runs.
      tmp <- iterdiag.nqtl[match(mainloci$niter, iterdiag$niter)]
      tmp <- unlist(tapply(tmp[match(mainloci$niter, iterdiag$niter)],
                           inter, mean))
      x[, "nqtl"] <- if(type.scan == "count") tmp else log10(tmp)
    }
  }
  else if(type.scan != "cellmean") {
    if(type.scan == "heritability") {
      for(i in scan.save) {
        x[, i] <- 100 * x[,i] / totvar
        x[is.na(x[, i]), i] <- 0
      }
    }
    else if(is.lod) {
      ## Number of model parameters.
      npar <- qb.npar(var1, var2, nfixcov, nrancov, intcov, iterdiag.nqtl,
                      iterdiag, mainloci, gbye, pairloci)

      ## Residual sum of squares (RSS) averaged over MCMC runs.
      tmp <- (nind.pheno - npar - 1) * iterdiag[, "envvar"]
      rss <- unlist(tapply(tmp[match(mainloci[, "niter"], iterdiag[, "niter"])],
                           inter, mean))
      rss[is.na(rss)] <- 0
      
      ## Keep npar for detection probability.
      if(type.scan == "detection" | type.scan == "npar") {
        ## Probability of detection given data.
        ## Number of parameters averages over MCMC runs.
        npar <- unlist(tapply(npar[match(mainloci[, "niter"], iterdiag[, "niter"])],
                              inter, mean))
        npar[is.na(npar)] <- 0
      }
      ## Calculate LPD, LR or deviance.
      ## mostly correct for LPD, but does it get LPD?
      ## also see ideas in Gaffney code
      for(i in scan.save) {
        ## This counts model df.
        ## It only counts epistasis once, even though loci may
        ## interact with multiple other loci.
        nscan <- switch(i, sum = length(scans) - 1,
                        main = length(var1),
                        epistasis = length(var2),
                        GxE = n.var1.covar,
                        1)
        x[, i] <- calc.objective(x[, i], rss, env, nind.pheno,
                                 nscan, npar, type.scan)
      }
    }
  }

  if(is.slice) {
    ## Add column for slice as mean of locus in slice.
    tmp <- (mainloci$chrom == slice["chr"] &
            mainloci[, "locus"] >= slice["start"] &
            mainloci[, "locus"] <= slice["end"])
    tmp <- tapply(mainloci[tmp, "locus"], mainloci[tmp, "niter"],
                  mean, na.rm = TRUE)
    tmp <- rep(tmp, c(table(mainloci[, "niter"])))
    tmp2 <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
    tmp2[is.na(tmp2)] <- mean(tmp2, na.rm = TRUE)
    tmp <- (grid$chr == slice["chr"] &
            grid$pos >= slice["start"] & grid$pos <= slice["end"])
    if(any(tmp))
      tmp2[tmp] <- NA
    x <- cbind(x, slice = tmp2)
    scan.save <- c(scan.save, "slice")
    
    ## Drop slice chromosome if not in chr list.
    if(is.na(match(slice["chr"], chr))) {
      tmp <- dimnames(x)
      x <- as.matrix(x[grid$chr != slice["chr"], ])
      dimnames(x) <- list(NULL, tmp[[2]])
      grid <- grid[grid$chr != slice["chr"], ]
      mainloci <- mainloci[mainloci$chrom != slice["chr"], ]
      pairloci <- pairloci[!(pairloci$chrom1 == slice["chr"] &
                             pairloci$chrom2 == slice["chr"]), ]
    }
  }

  ## Add objects used by plot or summary methods.
  if(sum.scan == "two")
    x
  else {
    vars <- if(is.bc) "varadd" else c("varadd","vardom")
    x <- list(one = x, grid = grid,
              mainloci = mainloci[, c("niter","chrom","locus",vars)],
              pairloci = pairloci[, c("niter","chrom1","locus1","chrom2","locus2")])

    ## Assign attributes passed to plot.qb.scanone.
    attr(x, "class") <- c("qb.scanone", "list")
    attr(x, "type.scan") <- type.scan
    attr(x, "scan") <- scan.save
    attr(x, "cross.class") <- qb.cross.class(qbObject)
    attr(x, "chr") <- chr
    attr(x, "min.iter") <- min.iter
    attr(x, "pheno.name") <- pheno.name
    attr(x, "geno.names") <- geno.names
    attr(x, "reference") <- mean(reference)
    attr(x, "step") <- qb.get(qbObject, "step")
    attr(x, "niter") <- qb.niter(qbObject)
    attr(x, "split.chr") <- split.chr

    ## We are not really using all the above. This can be simplified.
    if(is.slice)
      scan <- dimnames(x$one)[[2]]

    x <- qb.to.scanone(x, chr, smooth, scan.save, weight[1], split.chr,
                       center.type)
    if(is.slice)
      attr(x, "slice") <- slice
    class(x)[1] <- "qb.scanone"
    x
  }
}
###################################################################
calc.objective <- function(x, rss, env, nind.pheno, nscan, npar, type.scan)
{
  ## Empirical approximation to LPD.
  ## Watch for negative values.
  
  tmp <- rss + env * nscan + nind.pheno * x
  tmp[rss == 0] <- 0
  tmp2 <- tmp <= 0
  if(any(!tmp2))
    tmp[!tmp2] <- nind.pheno * log(tmp[!tmp2] / rss[!tmp2])
  if(any(tmp2))
    tmp[tmp2] <- 0
  x <- tmp
  
  if(type.scan == "LPD")
    x <- x / (2 * log(10))
  else if(type.scan == "LR")
    x <- x / 2
  else if(type.scan == "detection") {
    p1 <- exp(x / 2)
    detect.prior = 1 / length(x)
    x <- p1 * detect.prior / (1 + (p1 - 1) * detect.prior)
    x[is.na(x)] <- 0.5
  }
  else {
    if(type.scan == "rss")
      x <- rss
    else if(type.scan == "npar")
      x <- npar
  }
  x[is.na(x)] <- min(x, na.rm = TRUE)
  x
}
###################################################################
summary.qb.scanone <- function(object,
                               chr = NULL,
                               threshold = 0,
                               sort = "no",
                               n.qtl = 0.05,
                               ...)
{
  summary.qb.to.scanone(object, chr, threshold, sort, n.qtl, ...)
}
###################################################################
summary.qb.to.scanone <- function(object,
                               chr = NULL,
                               threshold = 0,
                               sort = "no",
                               n.qtl = 0.05,
                               ...)
{
  scan <- names(object)[-(1:2)]
  chrs <- attr(object, "chr")
#  count <- attr(object, "count")
  n.iter <- attr(object, "niter")
  centers <- attr(object, "centers")
  nqtl <- attr(object, "nqtl")

  geno.names <- levels(object$chr)
  chr.sub <- unclass(chrs)[match(object$chr, geno.names)] %in%
    qb.find.chr(chr = chr, geno.names = levels(chrs))
  object <- object[chr.sub, ]
  tmp <- table(object$chr) > 0
  geno.names <- names(tmp)[tmp]
  chrs <- chrs[tmp]
  nqtl <- nqtl[tmp]

  ## Set up out object; check for epistasis first.
  main.scan <- c("main","add","dom")

  ## qb.slicetwo uses type.scan as label for epistasis.
  type.scans <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","logposterior","2logBF","BF","nqtl",
             "npar","rss")
  ## qb.slicetwo also uses genotype.
  epi.scan <- c("epistasis","aa","ad","da","dd",
                "AA","AH","AB","HA","HH","HB","BA","BH","BB",
                type.scans)

  epistasis <- any(match(scan, epi.scan, nomatch = 0))
  out <- matrix(0, length(geno.names), 2 + ncol(object) + epistasis)
  object.names <- names(object)[-(1:2)]
  dimnames(out) <- list(geno.names,
                        c("chr", "n.qtl", "pos", "m.pos", "e.pos"[epistasis],
                          object.names))
  out <- as.data.frame(out)
  out[, "chr"] <- chrs

  ## Number of QTL per chr or split chr.
  out[, "n.qtl"] <- nqtl[geno.names]

  ## Positions of centers
  out[, "pos"] <- centers[geno.names, "pos"]
  out[, "m.pos"] <- centers[geno.names, "m.pos"]
  if(epistasis)
    out[, "e.pos"] <- centers[geno.names, "e.pos"]

  ## Use positions to find maximium value.
  for(i in seq(length(geno.names))) {
    ii <- (object$chr == geno.names[i] &
           !is.na(object$chr))
    if(any(ii)) {
      for(j in object.names) {
        wh.pos <- "pos"
        if(j %in% epi.scan)
          wh.pos <- "e.pos"
        else if(j %in% main.scan)
          wh.pos <- "m.pos"
        wh <- which.min(abs(object$pos[ii] - out[i, wh.pos]))[1]
        out[i, j] <- object[ii, j][wh]
      }
    }
  }

  ## Keep only chrs with some value about threshold.
  out <- qb.threshold(out, threshold)

  ## Print values rounded to digits places,
  ## ordered by sort column.
  if(is.null(out))
    NULL
  else {
    ## Restrict to pairs with at least n.qtl estimated QTL.
    tmp <- out[, "n.qtl"] >= n.qtl
    if(sum(tmp)) {
      out <- out[tmp,, drop = FALSE]
      if (match(sort, dimnames(out)[[2]], nomatch = 0) & nrow(out) > 1)
        out <- out[order(- out[, sort]), ]
      class(out) <- c("summary.qb.scanone", "data.frame")
      attr(out, "type.scan") <- attr(object, "type.scan")
      attr(out, "pheno.name") <- attr(object, "pheno.name")
      attr(out, "min.iter") <- attr(object, "min.iter")
      attr(out, "scan") <- scan
      attr(out, "threshold") <- threshold
    }
    else
      out <- NULL
  }
  out
}
###################################################################
print.qb.scanone <- function(x, digits = 3, ...)
  print(summary(x, ...), digits = 3)
###################################################################
print.summary.qb.scanone <- function(x, digits = 3, ...)
{
  cat(attr(x, "type.scan"), "of", attr(x, "pheno.name"), "for",
      paste(attr(x, "scan"), collapse = ","), "\n")
  min.iter <- attr(x, "min.iter")
  if(min.iter > 1)
    cat("Including only loci pairs with at least", min.iter, "samples.\n")
  threshold <- attr(x, "threshold")
  if(any(threshold != 0)) {
    cat("Thresholds:",
        paste(names(threshold), c(threshold), collapse = ", ", sep = "="),
        "\n")
  }
  cat("\n")

  print.data.frame(x[, -1], digits = digits)
}
###################################################################
plot.qb.scanone <- function(x, chr = NULL,
                            scan = scan.plots, ylim = ylims,
                            scan.name = scan.pretty, ...)
{
  geno.names <- attr(x, "geno.names")

  ## Figure out how to organize scans.
  scan.names <- names(x)[-(1:2)]
  scan.plots <-
    if(length(scan.names) < 5)
      rev(scan.names)
    else
      c("sum","main","epistasis")
  is.sum <- match("sum", scan.names, nomatch = 0)

  ## Fine subset that matches chr.
  ## Need to be tricky in case of split.chr not NULL.
  chr <- qb.find.chr(chr = chr, geno.names = levels(attr(x,"chr")))
  chr.sub <- unclass(attr(x, "chr"))[match(x$chr, geno.names)] %in% chr

  ## Automate separate plots by main, epistasis, sum.
  split.plots <- any(match(c("main","epistasis"), scan, nomatch = 0)) &
     length(scan.names) > 5

  if(split.plots) {
    scan.main <- scan.names[c(grep("add", scan.names),
                              grep("dom", scan.names))]
    scan.epis <- scan.names[match(c("aa","ad","da","dd"),
                                  scan.names, nomatch = 0)]
    if(!match("sum", scan, nomatch = 0))
      is.sum <- 0
    is.main <- length(scan.main) > 0 & match("main", scan, nomatch = 0)
    is.epis <- length(scan.epis) > 0 & match("epistasis", scan, nomatch = 0)
    tmpar <- par(mfrow = c((is.sum > 0) + is.main + is.epis, 1),
                 mar = c(3.1,4.1,2.1,0.1))
    scans <- scan.names[is.sum]
    if(is.main)
      scans <- c(scans,scan.main)
    if(is.epis)
      scans <- c(scans,scan.epis)
    scans <- scan.names[is.sum]
    if(is.main)
      scans <- c(scans, scan.main)
    if(is.epis)
      scans <- c(scans, scan.epis)

    ## Set limits to be the same for all scans.
    ylims <- range( c(x[chr.sub, scans]), na.rm = TRUE)

    ret <- NULL
    if(is.sum) {
      ret <- plot.qb.to.scanone(x, chr, scan.names[is.sum], ylim, "all effects", ...)
    }
    if(is.main) {
      tmp <- plot.qb.to.scanone(x, chr, scan.main, ylim, "main effects", ...)
      if(is.null(ret))
        ret <- tmp
      else {
        rnames <- c(row.names(ret), row.names(tmp))
        ret <- data.frame(color = c(as.character(ret$color),
                            as.character(tmp$color)),
                          linetype = c(ret$linetype, tmp$linetype))
        row.names(ret) <- rnames
      }
    }
    if(is.epis) {
      tmp <- plot.qb.to.scanone(x, chr, scan.epis, ylim, "epistatic effects", ...)
      if(is.null(ret))
        ret <- tmp
      else {
        rnames <- c(row.names(ret), row.names(tmp))
        ret <- data.frame(color = c(as.character(ret$color),
                            as.character(tmp$color)),
                          linetype = c(ret$linetype, tmp$linetype))
        row.names(ret) <- rnames
      }
    }
    par(tmpar)
    ret$color <- factor(ret$color)
  }
  else {
    ylims <- range( c(x[chr.sub, scan]), na.rm = TRUE)

    ## Work on pretty title.
    if(length(scan) < 4)
      scan.pretty <- paste(scan, collapse="+")
    else
      scan.pretty <- "effects"
    ret <- plot.qb.to.scanone(x, chr, scan, ylim, scan.pretty, ...)
  }

  invisible(ret)
}
###################################################################
plot.qb.to.scanone <- function(x,
                               chr = NULL,
                               scan = names(x)[-(1:2)],
                               ylim = ylims,
                               scan.name = scan.pretty,
                               col = NULL, lty = 1,
                               main = paste(type.scan, "of", pheno.name,
                                 "for", scan.name),
                               sub = subs,
                               verbose = FALSE,
                               add = FALSE,
                               ...)
{
  ## Work on pretty title.
  if(length(scan) < 4)
    scan.pretty <- paste(scan, collapse="+")
  else
    scan.pretty <- "effects"

  ## Print message about plot.
  if(verbose) {
    cat("\n", attr(x, "type.scan"), "of", pheno.name, "for",
        paste(scan, collapse = ","), "\n")
  }

  ## Process the selected scan terms.
  ylims <- range( c(x[, scan]), na.rm = TRUE)

  ## Figure out phenotype name indirectly.
  pheno.name <- attr(x, "pheno.name")
  type.scan <- attr(x, "type.scan")

  ## Find sum, if more than one scan, and color scheme.
  allscan <- length(scan) > 1

  ## Set up color scheme.
  scan.col <- function(x, allscan, supplied.col = NULL) {
    ## Default colors.
    if(!allscan) {
      cols <- "black"
      names(cols) <- x
    }
    else {
      col <- c(sum = "black",
               add = "blue", dom = "red",
               aa = "purple",
               ad = "green", da = "darkgreen",
               dd = "orange",
               main = "blue", epistasis = "purple", GxE = "darkred",
               A = "blue", H = "purple", B = "red",
               AA = "blue", AH = "purple", HA = "green", HH = "red")
      cols <- col[x]
      cols[grep(".add", x)] <- "darkblue"
      cols[grep(".dom", x)] <- "darkred"
      name.col <- names(cols)
      cols[is.na(cols)] <- "black"
      name.col[is.na(name.col)] <- x[is.na(name.col)]
      names(cols) <- name.col
    }
    if(!is.null(supplied.col)) {
      if(is.null(names(supplied.col))) {
        n.col <- length(supplied.col)
        names(supplied.col) <- array(x, n.col)
      }
      tmp <- match(names(supplied.col), names(cols), nomatch = 0)
      if(any(tmp > 0))
        cols[tmp] <- supplied.col[tmp > 0]
    }
    cols
  }

  col <- scan.col(scan, allscan, col)
  
  ## Set any missing colors to "black".
  tmp <- length(scan) - length(col)
  tmp2 <- names(col)
  if(tmp > 0) {
    if(is.null(tmp2)) {
      if(is.character(col))
        col <- c(col, rep("black", tmp))
      else
        col <- c(col, rep(1, tmp))
      names(col) <- scan
    }
    else {
      cols[tmp2] <- col
      col <- cols
    }
  }
  else { ## length(col) <= length(scan)
    if(is.null(tmp2)) {
      col <- col[seq(length(scan))]
      names(col) <- scan
    }
    else {
      tmp <- match(scan, tmp2, nomatch = 0)
      col <- col[tmp]
      tmp <- tmp == 0
      if(any(tmp)) {
        tmp2 <- names(col)
        if(is.character(col))
          col <- c(col, rep("black", sum(tmp)))
        else
          col <- c(col, rep(1, sum(tmp)))
        names(col) <- c(tmp2, scan[tmp])
      }
    }
  }
  subs <- NULL
  if(length(col) > 1 & !all(col == col[1]))
    subs <- paste(names(col), col, sep = "=")

  ## Set up line type.scans.
  if(is.numeric(lty)) {
    lty <- 1 + pmin(6, pmax(0, lty))
    lty <- c("blank", "solid", "dashed", "dotted", "dotdash",
             "longdash", "twodash")[lty]
  }
  if(length(lty) >= length(col))
    lty <- lty[seq(length(col))]
  else
    lty <- c(lty, rep("solid", length(col) - length(lty)))
  names(lty) <- names(col)
  if(length(lty) > 1 & !all(lty == lty[1])) {
    if(is.null(subs))
      subs <- paste(names(lty), lty, sep = "=")
    else
      subs <- paste(subs, lty)
  }
  subs <- paste(subs, collapse = ", ")
  
  ## Plot object after converting to scanone format.
  class(x) <- c("scanone", "data.frame")

  ## Add in breaks for split chromosomes.
  ## Change chr from numeric to character.
  orig.chr <- attr(x, "chr")
  x$chr <- orig.chr[unclass(x$chr)]
  geno.names <- levels(orig.chr)
  chr <- geno.names[chr]
  geno.names <- geno.names[geno.names %in% x$chr]
  chr <- geno.names[geno.names %in% chr]
  x$chr <- ordered(x$chr, geno.names)

  split.chr <- attr(x, "split.chr")
  split.chr <- split.chr[names(split.chr) %in% chr]
  if(length(split.chr)) { ## Some chr to be plotted are split.
    split.x <-
      data.frame(chr = ordered(rep(names(split.chr), sapply(split.chr, length)),
                   geno.names),
                 pos = unlist(split.chr))
    for(i in names(x)[-(1:2)])
      split.x[[i]] <- rep(NA, nrow(split.x))
    ## Make sure row names look like pseudomarkers, not markers.
    row.names(split.x) <- paste("c", as.character(split.x$chr), ".loc0",
                                seq(nrow(split.x)), sep = "")

    x <- rbind(x, split.x)
    x <- x[order(x$chr, x$pos), ]
  }

  for(i in seq(scan)) {
    scani <- scan[i]
    lodcolumn <- match(scani, names(x)) - 2
    if(i == 1)
      dimnames(x)[[2]][lodcolumn + 2] <- type.scan
    ## Call plot.scanone from R/qtl.
    plot(x, lodcolumn = lodcolumn, chr = chr, ..., add = (i > 1) | add,
         ylim = ylim, main = main,
         col = col[scani], lty = lty[scani])
    if(i == 1) {
      if(type.scan == "log10") {
        tmp <- c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000)
        axis(4,log10(tmp),tmp)
      }
      if(type.scan == "estimate")
        abline(h = 0, col = "grey", lty = 3, lwd = 2)
    }
  }
  if(allscan)
    if((length(col) < 5 | !missing(sub)) & sub != "")
      mtext(sub, 1, 2, cex = 0.65)

  ## Annotate axis and add vertical split if one chr and it is split.
  if(length(chr) == 1) {
    if(length(split.chr))
      abline(v = split.chr[[1]], col = "gray", lty = 2)
  }

  invisible(data.frame(color = col, linetype = lty))
}
###################################################################
qb.get.epis <- function(mainloci, pairloci, inter)
{
  if(is.null(pairloci))
    return(NULL)
  
  ## Epistasis counter.
  nqtl.pair <- c(paste(pairloci[, "niter"], pairloci[, "chrom1"], sep = ":"),
                 paste(pairloci[, "niter"], pairloci[, "chrom2"], sep = ":"))
  epinter <- make.chr.pos(nqtl.pair,
                          c(pairloci[, "locus1"], pairloci[, "locus2"]))

  ## epii identifies mainloci with epistatic pairs.
  epii <- match(unique(epinter),
                paste(mainloci[, "niter"], inter, sep = ":"),
                nomatch = 0)
  tmp <- rep(0, length(inter))
  tmp[epii] <- 1
  tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
  tmp[is.na(tmp)] <- 0
  tmp
}    
###################################################################
qb.get.main <- function(mainloci, inter)
{
  ## Main effects counter.
  vars <- c("varadd","vardom")
  vars <- vars[!is.na(match(vars, names(mainloci)))]
  tmp <- apply(as.matrix(mainloci[, vars]), 1,
               function(x) any(x > 0))
  tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
  tmp[is.na(tmp)] <- 0
  tmp
}
###################################################################
qb.centers <- function(object, center.type = c("mode","mean","scan"),
                       mainloci, pairloci, smooth = 3,
                       weight = "sqrt",
                       geno.names = levels(ordered(object$chr)),
                       inter, niter = unclass(table(inter)),
                       nepis = qb.get.epis(mainloci, pairloci, inter),
                       type.scan = attr(object, "type.scan"))
{
  ## Caution: inter must be constructed with object from the un-split chr
  ## so that it matches with mainloci. See call sequence in qb.to.scanone().

  center.type <- match.arg(center.type)
  if(type.scan %in% c("estimate","cellmean","nqtl","rss","npar") &
     center.type == "scan") {
    warning(paste("center.type reset to mode for", type.scan, "scans"))
    center.type <- "mode"
  }
  centers <- list()
  
  main.scan <- c("main","add","dom")
  epi.scan <- c("epistasis","aa","ad","da","dd")
  epistasis <- !is.null(nepis)

  if(center.type == "scan" & !missing(object)) {
    e.wh <- 0
    object.names <- names(object)[-(1:2)]
    epistasis <- epistasis & any(match(object.names, epi.scan, nomatch = 0))

    for(i in seq(length(geno.names))) {
      ii <- (object$chr == geno.names[i] &
             !is.na(object$chr))
      if(any(ii)) {
        tmp <- object.names %in% main.scan
        if(any(tmp)) {
          if("main" %in% object.names)
            m.wh <- which.max(object[ii, "main"])
          else {
            tmp <- object.names[tmp]
            m.wh <- which.max(apply(object[ii, tmp], 1, sum))
          }
        }
        else
          m.wh <- 0
        if(epistasis) {
          if("epistasis" %in% object.names)
            e.wh <- which.max(object[ii, "epistasis"])
          else {
            tmp <- object.names[object.names %in% epi.scan]
            e.wh <- which.max(apply(object[ii, tmp], 1, sum))
          }
          centers$e.pos[geno.names[i]] <- e.wh
        }
        wh <- {
          if("sum" %in% object.names)
            which.max(object[ii, "sum"])
          else ## Otherwise set pos to main or epistasis.
            ifelse(m.wh, m.wh, e.wh)
        }
      }
      centers$pos[geno.names[i]] <- wh
      centers$m.pos[geno.names[i]] <- m.wh
    }
  }
  else {
    nmain <- qb.get.main(mainloci, inter)
    
    if(center.type == "mean") {
      for(i in seq(length(geno.names))) {
        ii <- object$chr == geno.names[i]
        if(any(ii)) {
          pos <- weighted.mean(object$pos[ii], niter[ii])
          centers$pos[geno.names[i]] <- which.min(abs(object$pos[ii] - pos))

          tmp <- nmain[ii] > 0
          if(any(tmp)) {
            m.pos <- weighted.mean(object$pos[ii][tmp], nmain[ii][tmp])
            m.wh <- which.min(abs(object$pos[ii][tmp] - m.pos))
          }
          else
            m.wh <- wh
          centers$m.pos[geno.names[i]] <- m.wh
          
          if(epistasis) {
            tmp <- nepis[ii] > 0
            if(any(tmp)) {
              e.pos <- weighted.mean(object$pos[ii][tmp], nepis[ii][tmp])
              e.wh <- which.min(abs(object$pos[ii][tmp] - e.pos))
            }
            else
              e.wh <- wh
            centers$e.pos[geno.names[i]] <- e.wh
          }
        }
      }
    }
    else { ## default: center.type == "mode" or is.null(object) 
      tmp <- qb.smoothone(niter, object, smooth, niter, weight = weight)
      centers$pos <- unlist(tapply(tmp, object$chr, which.max))
      tmp <- qb.smoothone(nmain, object, smooth, niter, weight = weight)
      centers$m.pos <- unlist(tapply(tmp, object$chr, which.max))
      if(epistasis) {
        tmp <- qb.smoothone(nepis, object, smooth, nepis)
        centers$e.pos <- unlist(tapply(tmp, object$chr, which.max))
      }
    }
  }

  ## Now turn indices into chr positions.
  offset <- cumsum(c(0, unclass(table(object$chr))))
  offset <- offset[-length(offset)]
  names(offset) <- levels(object$chr)
  centers$pos <- object$pos[centers$pos + offset]
  centers$m.pos <- object$pos[centers$m.pos + offset]
  if(epistasis)
    centers$e.pos <- object$pos[centers$e.pos + offset]
  
  centers <- data.frame(centers)
  row.names(centers) <- geno.names
  centers
}
###################################################################
qb.to.scanone <- function(x,
                          chr = NULL,
                          smooth = 3,
                          scan = dimnames(x$one)[[2]],
                          weight = c("sqrt","count","none","atten","ratten"),
                          split.chr = attr(x, "split.chr"),
                          center.type = c("mode","mean","scan"),
                          ...)
{
  if(is.null(x))
    return(NULL)
  
  ## Prepare qb.scanone object as a scanone object.

  weight <- match.arg(weight)
  geno.names <- attr(x, "geno.names")
  reference <- attr(x, "reference")
  n.iter <- attr(x, "niter")
  center.type <- match.arg(center.type)

  ## Set up output grid.
  grid <- x$grid
  mainloci <- x$mainloci

  ## Subset to selected chromosomes.
  chr.sub <- grid$chr %in% chr
  grid <- grid[chr.sub, ]
  mainloci <- mainloci[mainloci$chrom %in% chr, ]
  one <- as.matrix(x$one[chr.sub, scan])
  dimnames(one) <- list(NULL, scan)

  ## Get interaction pattern.
  inter <- qb.inter(, grid, mainloci)

  ## Get number of samples per pos.
  niter <- unclass(table(inter))

  ## Reduce x$one to loci with at least min.iter samples.
  min.iter <- attr(x, "min.iter")
  if(min.iter > 1)
    one[niter < min.iter, ] <- 0

  ## Split chromosomes according to split.chr.
  xout <- qb.chrsplit(grid, mainloci, chr, n.iter, geno.names, split.chr)
  geno.names <- levels(xout$chr)
  chr <- attr(xout, "unsplit")

  ## Values at maximum number of smoothed iterations.
  epi.scan <- c("epistasis","aa","ad","da","dd")
  epistasis <- any(match(scan, epi.scan, nomatch = 0))
  nepis <- qb.get.epis(mainloci, x$pairloci, inter)

  n.qtl <- tapply(niter, xout$chr, sum) / n.iter

  ## Set up xout with scanone object attributes.
  scan <- dimnames(one)[[2]]
  ## Smooth over points and create scanone object scans.
  for(varcomp in scan) {
    tmp <- match(varcomp, epi.scan, nomatch = 0)
    xout[, varcomp] <- qb.smoothone(one[, varcomp], grid,
                                    smooth, if(tmp) nepis else niter,
                                    reference, weight = weight)
  }

  ## Find centers.
  centers <- qb.centers(xout, center.type, mainloci, x$pairloci,
                        smooth, weight, geno.names, inter, niter, nepis,
                        attr(x, "type.scan"))

  class(xout) <- c("qb.to.scanone", "scanone", "data.frame")
  attr(xout, "type.scan") <- attr(x, "type.scan")
  attr(xout, "model") <- "normal"
  attr(xout, "weight") <- weight
  attr(xout, "center.type") <- center.type
  attr(xout, "centers") <- centers
  attr(xout, "chr") <- ordered(attr(x, "geno.names")[chr], attr(x, "geno.names"))
  attr(xout, "min.iter") <- min.iter
  attr(xout, "niter") <- n.iter
  attr(xout, "nqtl") <- n.qtl
  attr(xout, "pheno.name") <- attr(x, "pheno.name")
  attr(xout, "geno.names") <- geno.names
  attr(xout, "split.chr") <- split.chr
  attr(xout, "cross.class") <- attr(x, "cross.class")

  xout
}
############################################################################## 
qb.smoothchr <- function(x, smooth, niter, reference = 0, weight = "sqrt")
{
  ## Weighted average of x.
  switch(weight,
         count = {w <- niter},
         none = {w <- rep(1, length(x))},
         sqrt =, {w <- sqrt(niter)})
  
  nmap <- length(x)
  re.na <- is.na(x)
  x[re.na] <- reference
  if(nmap > 3) {
    o <- (x != reference & !is.na(x))
    for(i in seq(smooth)) {
      wtlod <- w * x
      x <- wtlod[c(1, seq(nmap - 1))] + wtlod[c(seq(2, nmap), nmap)]
      
      ## Double weight at observation if not zero.
      if(any(o))
        x[o] <- x[o] + 2 * wtlod[o]
      wtlod <- w[c(1, seq(nmap - 1))] + w[c(seq(2, nmap), nmap)]
      if(any(o))
        wtlod[o] <- wtlod[o] + 2 * w[o]
      x <- x / wtlod
    }
    x[is.na(x)] <- reference
  }
  x[re.na] <- NA
  x
}
##############################################################################
make.atten <- function(pos, smooth, weight)
{
  wt <- exp(-as.matrix(dist(pos)) / smooth)
  if(weight == "ratten") {
    ## Use weight matrix as sqrt of distances.
    wt <- svd(wt)
    wt <- wt$u %*% diag(sqrt(wt$d)) %*% t(wt$v)
  }
  t(apply(wt, 1, function(x,y) x / y, apply(wt, 2, sum)))
}
qb.smoothone <- function(x, grid, smooth, niter, reference = 0,
                         weight = "sqrt")
{
  if(smooth) {
    ## Attenuation smoothes using distance (assumes 3 as default).
    ## Rescales to 100 * cM.
    if(weight %in% c("atten","ratten"))
      smooth <- smooth * 2 / 3
    for(chr in unique(grid$chr)) {
      rows <- chr == grid$chr
      if(sum(rows)) {
        if(weight %in% c("atten","ratten")) {
          ## This approach ignores niter, which might be important.
          wt <- make.atten(grid[rows, 2], smooth, weight)
          x[rows] <- matrix(x[rows], 1) %*% wt
        }
        else
          x[rows] <- qb.smoothchr(x[rows], smooth, niter[rows], reference,
                                  weight)
      }
    }
  }
  x
}
##############################################################################
qb.indextwo <- function(iterdiag, mainloci, nqtl = qb.nqtl(, iterdiag, mainloci))
{
  ## index of pairs of loci for each iteration
  ## need iterdiag.nqtl as well as mainloci!
  unlist(apply(as.matrix(seq(nrow(mainloci))[!duplicated(mainloci[, "niter"])]),
               1,
               function(x,y) {
                 n <- y[x]
                 if(n > 1) {
                   values <- seq(x, len = n)
                   ## get all possible pairs of values
                   j <- matrix(values, n, n)
                   j <- j[row(j) > col(j)]
                   rbind(rep(values[-n], (n-1):1), j)
                 }
                 else
                   matrix(0, 2, 0)
               },
               nqtl[match(mainloci[, "niter"], iterdiag[, "niter"])]))
}
##############################################################################
qb.intertwo <- function(min.iter = 1, mainloci, iterdiag, pairloci)
{
  nqtl <- qb.nqtl(, iterdiag, mainloci)
  index <- qb.indextwo(iterdiag, mainloci, nqtl)
  gridtwo <- matrix(t(mainloci[index, c("chrom","locus")]), 4)
  inter <- paste(gridtwo[1, ], gridtwo[2, ],
                 gridtwo[3, ], gridtwo[4, ], sep = ":")
  inter <- ordered(inter, inter[!duplicated(inter)])
  
  ## Set up epistatic count, which requires several other things.
  nqtl <- nqtl * (nqtl - 1) / 2
  epi <- rep(0, length(inter))
  epi[match(paste(pairloci[, "niter"],
                  pairloci[, "chrom1"], pairloci[, "locus1"],
                  pairloci[, "chrom2"], pairloci[, "locus2"],
                  sep = ":"), 
            paste(rep(iterdiag[, "niter"], nqtl),
                  as.character(inter),
                  sep = ":"))] <- 1

  ## Return columns as chrom1, locus1, chrom2, locus2, niter, nepis.
  gridtwo <- rbind(gridtwo[, !duplicated(inter)],
                   niter = unclass(table(inter)),
                   nepis = unlist(tapply(epi, inter, sum)))
  gridtwo[, gridtwo["niter", ] >= min.iter]
}
###################################################################
qb.scantwo <- function(qbObject, epistasis = TRUE,
                       scan = list(upper = upper.scan, lower = lower.scan),
                       type.scan = c(
                         upper = "heritability",
                         lower = "heritability"),
                       upper.scan = "epistasis",
                       lower.scan = "full",
                       covar = {
                         if(nfixcov) seq(qb.get(qbObject, "nfixcov"))
                         else 0},
                       adjust.covar = NA,
                       chr = NULL,
                       min.iter = 1,
                       verbose = FALSE, ...)
{
  qb.exists(qbObject)
  
  ## Need to add aggregate facilities and redo counts as in qb.scanone.

  ## Following prior used for Bayes factors.
  ## Need to do this before subsetting on chr.
  bf.prior <- qb.get(qbObject, "mean.nqtl") / qb.nloci(qbObject)
  bf.prior <- bf.prior * bf.prior
  
  qb.name <- deparse(substitute(qbObject))
  chr <- qb.find.chr(qbObject, chr)
  qbObject <- subset(qbObject, chr = chr, restrict.pair = FALSE)
  
  nfixcov <- qb.get(qbObject, "nfixcov")
  nrancov <- qb.get(qbObject, "nrancov")
  intcov <- as.logical(qb.get(qbObject, "intcov"))
  intcov <- check.intcov(intcov, nfixcov)

  pairloci <- qb.get(qbObject, "pairloci", ...)
  if(is.null(pairloci))
    epistasis <- FALSE

  ## Determine type of qb.scan.
  type.scans <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","logposterior","BF","2logBF","nqtl")
  tmp <- names(type.scan)
  type.scan <- type.scans[pmatch(tolower(type.scan), tolower(type.scans), nomatch = 2,
                       duplicates.ok = TRUE)]
  type.scan <- array(type.scan, 2)
  if(is.null(tmp))
    tmp <- c("upper","lower")
  names(type.scan) <- tmp

  if(any(type.scan == "cellmean"))
    stop("cellmean type not implemented: use qb.slice")
  is.count <- match(type.scan,
                    c("count", "log10", "posterior", "logposterior",
                      "BF", "2logBF", "nqtl"),
                    nomatch = 0)
  names(is.count) <- names(type.scan)

  is.var <- match(type.scan, type.scans[1:6], nomatch = 0)
  is.effect <- is.var | type.scan == "estimate"
  is.lod <- match(type.scan, type.scans[2:5], nomatch = 0)
  names(is.var) <- names(is.effect) <- names(is.lod) <- names(type.scan)
  
  ## Number of individuals for phenotype.
  cross <- qb.cross(qbObject, genoprob = FALSE)
  tmp <- list(...)
  if("pheno.col" %in% names(tmp))
    pheno.col <- tmp$pheno.col
  else
    pheno.col <- qb.get(qbObject, "pheno.col")
  pheno.name <- qb.pheno.names(qbObject, cross)[pheno.col[1]]
  nind.pheno <- qb.nind.pheno(qbObject, pheno.name, nfixcov, cross)

  ## Genotype names.
  map <- pull.map(cross)
  geno.names <- qb.geno.names(qbObject, cross)
  rm(cross)
  gc()
  
  ## Get MCMC samples.
  iterdiag <- qb.get(qbObject, "iterdiag", ...)
  mainloci <- qb.get(qbObject, "mainloci", ...)
  if(is.null(mainloci))
    return(NULL)
  
  iterdiag.nqtl <- qb.nqtl(qbObject, iterdiag, mainloci)
  npair <- iterdiag.nqtl
  npair <- npair * (npair - 1) / 2

  ## Determine variance components.
  is.bc <- (qb.cross.class(qbObject) == "bc")
  var1 <- "add"
  if(epistasis) {
    var2 <- "aa"
    if(!is.bc) {
      var1 <- c(var1,"dom")
      var2 <- c(var2,"ad","da","dd")
    }
  }
  else
    var2 <- character()

  ## Set up index into mainloci for pairs.
  index <- qb.indextwo(iterdiag, mainloci, iterdiag.nqtl)
  nindex <- length(index) / 2

  ## Find all pairs of loci.
  if(verbose)
    cat("finding all pairs of loci ...\n")
  tmp <- matrix(t(mainloci[index, c("chrom","locus")]), 4)
  inter <- paste(tmp[1, ], tmp[2, ],
                 tmp[3, ], tmp[4, ], sep = ":")
  inter <- ordered(inter, inter[!duplicated(inter)])

  ## Set up index for number of linked qtl.
  if(any(type.scan == "nqtl")) {
                                        #    nqtl.main <- paste(mainloci[index[seq(by = 2, length = nindex)], "niter"],
                                        #                       tmp[1, ], tmp[3, ], sep = ":")
    nqtl.main <- paste(mainloci[, "niter"], mainloci[, "chrom"], sep = ":")
  }

  ## Covariate adjustment calculations.
  if(any(type.scan == "heritability"))
    totvar <- rep(0, length(inter))
  if(nfixcov) {
    ## Covariate means.
    covar.means <- covar.mean(qbObject, adjust.covar,
                              verbose = verbose & (any(type.scan == "estimate")),
                              ...)

    ## Explained covariance for heritability.
    if(any(type.scan == "heritability"))
      totvar <- rep(apply(qb.varcomp(qbObject, c("fixcov","rancov")),
                          1, sum),
                    npair)
  }

  ## Set up lists (elements upper, lower) of scan names.
  ## Scan can be several variance components (default is all).
  ## upper, lower = original or default scan names
  ##   ("epistasis","full" or "main","GxE").
  ## scan.names = names retained for returned object.
  ## scan = all terms needed for analysis (var1, var2, and var1.covar=GxE).
  if(missing(scan))
    scan <- list(upper = upper.scan, lower = lower.scan)
  if(!is.list(scan))
    scan <- as.list(scan)
  if(length(scan) != 2)
    stop("scan must be list of length two with names upper and lower")
  if(any(is.na(match(names(scan), c("lower","upper")))))
    stop("scan list names must be upper and lower")
  
  ## Match key words if present (but convert joint to full).
  for(tri in c("lower","upper")) {
    keys <- c("main","epistasis","full","joint")
    keysfull <- c("main","epistasis","full","full")
    tmp <- pmatch(scan[[tri]], keys, nomatch = 0, duplicates.ok = TRUE)
    if(any(tmp))
      scan[[tri]][tmp > 0] <- keysfull[tmp]
  }
  ## If no epistasis and missing upper, then set upper scan to main.
  if(!epistasis) {
    if(missing(scan) & missing(upper.scan))
      scan[["upper"]] <- "main"
  }
  ## Save scan names for later plot.
  if(any(type.scan == "estimate")) {
    if(missing(scan) & missing(upper.scan) & type.scan["upper"] == "estimate")
      scan$upper <- "aa"
    if(missing(scan) & missing(lower.scan) & type.scan["lower"] == "estimate")
      scan$lower <- "add"
  }
  scan.names <- scan
  ## Now convert scan to the terms needed for analysis.
  for(tri in c("lower","upper")) {
    scan[[tri]] <- switch(scan[[tri]],
                          main = {
                            var1
                          },
                          epistasis = {
                            if(epistasis)
                              var2
                            else
                              var1
                          },
                          full = {
                            c(var1, var2)
                          },
                          scan[[tri]])
    vars <- var1[match(scan[[tri]], var1, nomatch = 0)]
    if(sum(intcov)) {
      if(length(covar.means) & length(vars)) {
        tmp <- seq(nfixcov)[intcov]
        tmp <- names(covar.means)[covar[match(tmp, covar, nomatch = 0)]]
        if(length(tmp))
          scan[[tri]] <- unique(c(scan[[tri]],
                                  outer(vars, tmp, paste, sep = ".")))
      }
    }
  }

  ## Weight for summary.qb.scantwo
  weight <- character()
  for(tri in c("lower","upper")) {
    weight[tri] <- if(any(match(scan[[tri]],
                                c("main","sum","full","add","dom"),
                                nomatch = 0)))
      "main"
    else
      "epistasis"
  }

  if(verbose) {
    cat(paste(c("\nupper:","lower:"), type.scan[c("upper","lower")], "of",
              pheno.name, "for",
              sapply(scan, paste, collapse = "+")[c("upper","lower")],
              collapse = "\n"),
        "\n")
    if(min.iter > 1)
      cat("Including only loci pairs with at least", min.iter, "samples.\n")
    cat("\n")
  }

  ## Extract environmental variance.
  if(any(type.scan == "heritability") | any(is.lod)) {
    if(verbose)
      cat("environmental variance ...\n")
    tmp <- rep(iterdiag[, "envvar"], npair)
    if(any(type.scan == "heritability"))
      totvar <- totvar + tmp
    else if(any(is.lod))
      env <- unlist(tapply(tmp, inter, mean))
  }

  var.elem <- function(type.scan, vari) {
    ifelse(type.scan == "estimate", vari, paste("var", vari, sep = ""))
  }

  accum <- matrix(0, nindex, 2)
  dimnames(accum) <- list(NULL, c("lower","upper"))
  
  ## Number of main effect samples per locus.
  is.full <- (is.count &
              sapply(scan, function(x) match("full", x, nomatch = 0)))
  names(is.full) <- names(is.count)
  if(any(is.full)) {
    for(tri in names(is.full)[is.full])
      accum[, tri] <- 1
  }
  
  ## Get non-epistatic components: additive and dominance.
  if(any(type.scan == "heritability"))
    vars <- var1
  else
    vars <- var1[match(unique(unlist(scan)), var1, nomatch = 0)]
  if(verbose & length(vars))
    cat("non-epistatic components ...\n")
  for(vari in vars) {
    if(any(type.scan == "estimate"))
      main.val <- mainloci[, vari]
    
    ## Loop over all interacting covariates.
    if(sum(intcov)) {
      ## Get GxE samples.
      gbye <- qb.get(qbObject, "gbye", ...)
      
      cov.val <- rep(0, nrow(mainloci))
      covars <- seq(nfixcov)[intcov]
      for(covj in covars) {
        gbyej <- gbye[gbye[, "covar"] == covj, ]
        if(length(gbyej)) {
          same <- match(paste(gbyej[, "niter"], gbyej[, "chrom"],
                              gbyej[, "locus"], sep = ":"),
                        paste(mainloci[, "niter"], mainloci[, "chrom"],
                              mainloci[, "locus"], sep = ":"))
          
          ## Parameter estimates of GxE fixed effects.
          done.her <- done.est <- FALSE
          for(tri in c("lower","upper")) {
            if(match(covj, covar, nomatch = 0) |
               type.scan[tri] == "heritability") {
              if(tri == "lower" | type.scan[tri] != type.scan["lower"]) {
                cname <- paste(vari, paste(names(covar.means)[covj],
                                           sep = "."))
                cov.val[same] <- gbyej[[var.elem(type.scan[tri], vari)]]
                tmp <- apply(array(cov.val[index], c(2,nindex)), 2, sum)
              }
              if(type.scan[tri] == "heritability" & !done.her) {
                totvar <- totvar + tmp
                done.her <- TRUE
              }
              if(match(cname, scan[[tri]], nomatch = 0))
                accum[, tri] <- accum[, tri] + tmp
            }
            
            ## Offset parameter estimate by covariates.
            if((type.scan[tri] == "estimate" | type.scan[tri] == "cellmean") &
               covar.means[covj] != 0 & !done.est) {
              main.val[same] <- main.val[same] + covar.means[covj] * gbyej[[vari]]
              done.est <- TRUE
            }
          }
        }
      }
    }

    ## Now include the component.
    done.her <- FALSE
    for(tri in c("lower","upper")) {
      if(tri == "lower" | type.scan[tri] != type.scan["lower"]) {
        if(type.scan[tri] == "estimate")
          tmp <- main.val
        else ## variance components and counts
          tmp <- mainloci[, paste("var", vari, sep = "")]
        tmp <- apply(array(tmp[index], c(2,nindex)), 2, sum)
      }
      if(type.scan[tri] == "heritability" & !done.her) {
        totvar <- totvar + tmp
        done.her <- TRUE
      }
      if(match(vari, scan[[tri]], nomatch = 0))
        accum[, tri] <- accum[, tri] + tmp
    }
  }

  ## Epistatic components.
  if(epistasis) {
    if(verbose)
      cat("epistatic components ...\n")
    if(any(type.scan == "heritability"))
      vars <- var2
    else
      vars <- var2[match(unique(unlist(scan)), var2, nomatch = 0)]

    if(length(vars)) {
      epi.match <- match(paste(pairloci[, "niter"],
                               pairloci[, "chrom1"], pairloci[, "locus1"],
                               pairloci[, "chrom2"], pairloci[, "locus2"],
                               sep = ":"), 
                         paste(rep(iterdiag[, "niter"], npair),
                               as.character(inter),
                               sep = ":"),
                         nomatch = 0)
      if(any(type.scan == "nqtl")) {
        nqtl.pair <- c(paste(pairloci[, "niter"], pairloci[, "chrom1"],
                             pairloci[, "chrom2"], sep = ":"))
      }

      ## Number of epistatic samples per locus.
      is.epis <- (is.count &
                  sapply(scan, function(x) match("epistasis", x, nomatch = 0)))
      names(is.epis) <- names(is.count)
      if(any(is.epis))
        accum[epi.match, names(is.epis)[is.epis]] <- 1

      ## Get element and average.
      for(vari in vars) {
        done.her <- FALSE
        for(tri in c("lower","upper")) {
          element <- var.elem(type.scan[tri], vari)
          if(type.scan[tri] == "heritability" & !done.her) {
            totvar[epi.match] <-
              totvar[epi.match] + pairloci[epi.match > 0, element]
            done.her <- TRUE
          }
          if(match(vari, scan[[tri]], nomatch = 0))
            accum[epi.match, tri] <-
              accum[epi.match, tri] + pairloci[epi.match > 0, element]
        }
      }
    }
  }

  ## Need n.iter for counts.
  if(any(is.count))
    n.iter <- nrow(iterdiag)

  ## Sum counts and/or average effects for each loci pair.
  n.inter <- length(levels(inter))
  for(tri in c("lower","upper")) {
    if(is.count[tri]) {
      if(type.scan[tri] == "nqtl") {
        if(weight[tri] == "main") {
          tmp <- unlist(table(nqtl.main))[nqtl.main]
          tmp2 <- accum[, tri] > 0
          accum[, tri] <- apply(array(tmp[index], c(2,nindex)), 2, sum)
          tmp <- apply(array(nqtl.main[index], c(2,nindex)), 2,
                       function(x) identical(x[1], x[2]))
          accum[tmp, tri] <- accum[tmp, tri] / 2
          accum[!tmp2, tri] <- 0
          rm(tmp2)
          gc()
                                        #         accum[, tri] <- unlist(tapply(accum[, tri] > 0, nqtl.main, sum))[nqtl.main]
        }
        else { ## epistasis
          accum[epi.match, tri] <-
            unlist(tapply(accum[epi.match, tri] > 0,
                          nqtl.pair, sum))[nqtl.pair][epi.match > 0]
        }
        tmp <- unlist(tapply(accum[, tri], inter, mean, na.rm = TRUE))
      }
      else ## other counts
        tmp <- unlist(tapply(accum[, tri] > 0, inter, sum))
      tmp[is.na(tmp)] <- 0
      accum[seq(n.inter), tri] <-
        qb.count(tmp, type.scan[tri], n.iter, bf.prior)
    }
    else { ## is.effect
      tmp <- unlist(tapply(accum[, tri], inter, mean))
      tmp[is.na(tmp)] <- 0
      accum[seq(n.inter), tri] <- tmp
    }
  }
  accum <- accum[seq(n.inter), ]
  
  if(any(type.scan == "heritability")) {
    totvar <- unlist(tapply(totvar, inter, mean))
    totvar[is.na(totvar)] <- 0
  }

  
  ## Compute heritability.
  if(any(type.scan == "heritability") & verbose)
    cat("heritability ...\n")
  for(tri in c("lower","upper")) {
    if(type.scan[tri] == "heritability") {
      accum[, tri] <- 100 * accum[, tri] / totvar
      accum[, tri][is.na(accum[, tri])] <- 0
    }
  }
  if(any(is.lod)) {
    if(verbose)
      cat("LPD or other diagnostics ...\n")

    ## Number of model parameters.
    npar <- qb.npar(var1, var2, nfixcov, nrancov, intcov, iterdiag.nqtl,
                    iterdiag, mainloci, gbye, pairloci)

    ## Residual sum of squares.
    rss <- (nind.pheno - npar - 1) * iterdiag[, "envvar"]
    ## there is probably a better way to connect mainloci[, "niter"] to inter
    tmp <- c(array(mainloci[index, "niter"], c(2,nindex))[1, ])
    tmp <- match(tmp, iterdiag[, "niter"])
    rss <- unlist(tapply(rss[tmp], inter, mean))
    
    if(any(type.scan == "detection")) {
      npar <- unlist(tapply(npar[tmp], inter, mean))
      npar[is.na(npar)] <- 0
      for(tri in c("lower","upper")) {
      }
    }
    
    ## calculate LPD or other diagnostic.
    for(tri in c("lower","upper")) {
      ## Count df for main effects twice, for both loci.
      ## Only count interacting covariates once (could be a mistake).
      ## Would be more involved to count these properly.
      tmp <- sum(!is.na(match(scan[[tri]], var1)))
      nscan <- tmp + length(scan[[tri]])
      if(is.lod[tri])
        accum[, tri] <- calc.objective(accum[, tri], rss, env, nind.pheno,
                                       nscan, npar, type.scan[tri])
    }
  }

  ## Keep only samples with at least min.iter iterations.
  accum <- accum[unclass(table(inter)) >= min.iter, ]

  ## Diagonal from qb.scanone.
  if(is.count["lower"])
    scan.one <- scan$lower[1]
  else
    scan.one <- scan$lower[unlist(apply(as.matrix(var1),1,grep,scan$lower))]
  qb.scan <- list(two = accum)
  grid <- pull.grid(qbObject, offset = TRUE, spacing = TRUE)
  if(length(scan.one)) {
    if(verbose)
      cat("qb.scanone on diagonal with", paste(scan.one, collapse = ","),
          "...\n")
    qb.scan$one <- qb.scanone(qbObject, epistasis,
                              scan.one, type.scan["lower"], sum.scan = "two",
                              covar = covar, min.iter = min.iter,
                              verbose = FALSE)
  }
  else {
    if(verbose)
      cat("diagonal set to zero\n")
    qb.scan$one <- rep(0, nrow(grid))
  }

  qb.scan$grid <- grid
  qb.scan$iterdiag <- iterdiag
  qb.scan$mainloci <- mainloci
  qb.scan$pairloci <- pairloci

  ## Assign attributes.
  attr(qb.scan, "class") <- c("qb.scantwo", "list")
  attr(qb.scan, "type.scan") <- type.scan
  attr(qb.scan, "scan") <- scan.names
  attr(qb.scan, "min.iter") <- min.iter
  attr(qb.scan, "cross.class") <- qb.cross.class(qbObject)
  attr(qb.scan, "chr") <- chr
  attr(qb.scan, "weight") <- weight
  attr(qb.scan, "pheno.name") <- pheno.name
  attr(qb.scan, "geno.names") <- geno.names
  attr(qb.scan, "map") <- map
#  attr(qb.scan, "qb") <- qb.name
  attr(qb.scan, "niter") <- qb.niter(qbObject)
  attr(qb.scan, "reference") <- mean(qb.reference(qbObject,
                                                  qb.scan$mainloci,
                                                  qb.scan$iterdiag,
                                                  inter, type.scan))
  attr(qb.scan, "split.chr") <- qb.get(qbObject, "split.chr")
  qb.scan
}
###################################################################
summary.qb.scantwo <- function(object,
                               chr = NULL, ## Must be integer for now.
                               threshold = 0,
                               sort = "no",
                               which.pos = "upper",
                               min.iter = attr(object, "min.iter"),
                               refine = FALSE, width = 10, smooth = 3,
                               n.qtl = 0.05,
                               weight = c("sqrt","count","none","atten","ratten"),
                               ...)
{
  ## new intertwo needs to be checked out
  ## need pos1 and pos2 for lower and upper separately
  ## chr not working?
  weight <- match.arg(weight)
  
  pheno.name <- attr(object, "pheno.name")

  ## Get position pairs.
  gridtwo <- qb.intertwo(min.iter, object$mainloci, object$iterdiag, object$pairloci)

  geno.names <- attr(object, "geno.names")
  inter <- paste(geno.names[gridtwo[1, ]], geno.names[gridtwo[3, ]], sep = ":")

  ## Get unique pairs of chromosomes.
  ## This assumes chromosome names are unique!
  tmp <- order(gridtwo[1, ], gridtwo[3, ])
  chr.pair <- unique(inter[tmp])
  chrs <- as.matrix(gridtwo[c(1,3), tmp[!duplicated(inter[tmp])]])
  chr <- qb.find.chr(chr = chr, geno.names = geno.names)

  tmp <- !is.na(match(chr.pair,
                      c(outer(geno.names[chr], geno.names[chr], paste, sep = ":"))))
  chr.pair <- chr.pair[tmp]
  chrs <- as.matrix(chrs[, tmp])
  keep <- !is.na(match(gridtwo[1, ], chr)) & !is.na(match(gridtwo[3, ], chr))

  out <- matrix(0, length(chr.pair), 9)
  dimnames(out) <- list(chr.pair, c("chr1", "chr2",
                                    "n.qtl",
                                    "l.pos1", "l.pos2", "lower",
                                    "u.pos1", "u.pos2", "upper"))

  n.iter <- attr(object, "niter")

  out[, c("chr1", "chr2")] <- t(chrs)

  tmp <- tapply(gridtwo["niter", keep], inter[keep], sum) / n.iter
  out[, "n.qtl"] <- out[,"n.qtl"] <- tmp[chr.pair]
  rm(keep)

  ## Restrict to pairs with at least n.qtl estimated QTL.
  tmp <- out[, "n.qtl"] >= n.qtl
  if(sum(tmp)) {
    out <- out[tmp,, drop = FALSE]
    chr.pair <- chr.pair[tmp]
  }
  else
    out <- NULL
  
  if(!is.null(out)) {
    x2 <- qb.scantwo.smooth(object, chr, smooth,
                            qb.intertwo(min.iter, object$mainloci, object$iterdiag,
                                        object$pairloci),
                            weight = weight, ...)

    type.scan <- attr(x2, "type.scan")

    ## Center as mean for variance, estimate, cellmean.
    ## Center as mode for all other types of scans.
    center <- character()
    for(tri in names(type.scan)) {
      center[tri] <- "mean"
      if(is.na(match(type.scan[tri], c("variance","estimate","cellmean"))))
        center[tri] <- "mode"
    }
    
    ## Weighted means by chr.
    tmp <- upper.tri(x2$lod)
    tmpx <- x2
    tmpx$lod[tmp] <- t(x2$lod)[tmp] - x2$lod[tmp]
    tmpx$map$chr <- ordered(geno.names[tmpx$map$chr], geno.names)
    tmp <- summary(tmpx)
    tmp2 <- paste(tmp$chr1, tmp$chr2, sep = ":")
    tmp2 <- match(chr.pair, tmp2)
    ## Fix any NA, due to reversal of chr1 and chr2.
    if(any(is.na(tmp2))) {
      tmp3 <- paste(geno.names[tmp$chr2], geno.names[tmp$chr1], sep = ":")
      tmp2[is.na(tmp2)] <- match(chr.pair[is.na(tmp2)], tmp3)
    }
    
    ## R/qtl column names changing with 1.04-48.
    if(compareVersion(qtlversion(), "1.04-48") < 0)
      stop("old version of R/qtl: please update now")
    
    ## The following uses R/qtl's summary.scanone to get mode
    ## for each pair of chromosomes in upper and/or lower triangle.
    ## Want to have option to get mean (weighted by gridtwo[,"nepis"]).
    
    ## Lower triangle (full).
    for(i in c("pos1","pos2"))
      out[, paste("l", i, sep = ".")] <- tmp[tmp2, paste(i, "f", sep = "")]
    out[, "lower"] <- tmp[tmp2, "lod.full"]

    ## Upper triangle (int).
    tmp <- summary(tmpx, what = "int")
    for(i in c("pos1","pos2"))
      out[, paste("u", i, sep = ".")] <- tmp[tmp2, i]
    out[, "upper"] <- tmp[tmp2, "lod.int"]
    
    rm(tmp,tmp2,tmpx)
    gc()
    
    ## Drop loci pairs with all zeroes.
    if(nrow(out) > 1) {
      keep <- apply(out[, 5:6], 1, function(x) !all(x == 0))
      if(sum(keep) > 1)
        out <- out[keep, ]
      else {
        dim.out <- dimnames(out)
        dim.out[[1]] <- dim.out[[1]][keep]
        matrix(out[keep,], 1, dimnames = dim.out)
      }
    }
    
    ## Keep only chrs with some value about threshold.
    out <- qb.threshold(out, threshold, 2)
  }
  
  ## Refine estimates for LPD.
  if(!is.null(out) & refine) if(any(center == "mode")) {
    map <- attr(object, "map")
    n.sum <- nrow(out)

    if(n.sum) for(i in seq(n.sum)) {
      chr <- as.vector(out[i,1:2])
      pos <- as.vector(out[i,3:4])
      for(j in 1:2) {
        ## Make sure you use 2*width at ends.
        rng <- range(map[[chr[j]]])
        pos[j] <- max(rng[1] + width - 0.1,
                      min(rng[2] - width + 0.1, pos[j]))
      }

      for(j in 1:2) {
        ## Refine both triangular parts.
        for(tri in c("upper","lower")) {
          grid <- qb.scantwo.slice(x2, chr[j],
                                   slice=c(chr=chr[3-j],
                                     start = pos[3-j] - width,
                                     end = pos[3-j] + width,
                                     upper = (tri == "upper")),
                                   type.scan, smooth, weight)
          tmp <- max(grid[, 3])
          if(tmp > out[i, tri]) {
            ## Return position for the chosen triangular part.
            out[i, c(paste("chr", j, sep = ""),
                     paste(substring(tri, 1, 1), ".pos", j, sep = ""))] <-
                       unlist(grid[which.max(grid[, 3]), c("chr","pos")])
            ## Update LPD or other mode.
            out[i, tri] <- tmp
          }
        }
      }
    }
  }
  
  ## Order by sort column.
  if(!is.null(out)) {
    if(nrow(out) > 1 & match(sort, dimnames(out)[[2]], nomatch = 0))
      out <- out[order(- out[, sort]), ]
  }
  
  ## Print values rounded to digits places,
  ## ordered by sort column.
  if(!is.null(out)) {
    out <- as.data.frame(out)
    class(out) <- c("summary.qb.scantwo", "data.frame")
    attr(out, "type.scan") <- type.scan
    attr(out, "pheno.name") <- pheno.name
    attr(out, "scan") <- attr(x2, "scan")
    attr(out, "min.iter") <- min.iter
    attr(out, "threshold") <- threshold
  }
  out
}
###################################################################
print.qb.scantwo <- function(x, ...) print(summary(x, ...), ...)
###################################################################
print.summary.qb.scantwo <- function(x, digits = 3, ...)
{
  z <- as.character(unlist(x[, 1]))
  if (max(nchar(z)) == 1) 
    rownames(x) <- apply(x[, 1:2], 1, function(a) {
      paste("c", a, collapse = ":", sep = "")
    })
  else rownames(x) <- apply(x[, 1:2], 1, function(a) {
    paste(sprintf("c%-2s", a), collapse = ":")
  })
  
  cat(paste(c("upper:","lower:"), attr(x, "type.scan")[c("upper","lower")],
            "of", attr(x, "pheno.name"),
            "for", attr(x, "scan")[c("upper","lower")],
            collapse = "\n"),
      "\n")
  min.iter <- attr(x, "min.iter")
  if(min.iter > 1)
    cat("Including only loci pairs with at least", min.iter, "samples.\n")
  threshold <- attr(x, "threshold")
  if(any(threshold != 0)) {
    cat("Thresholds:",
        paste(names(threshold), c(threshold), collapse = ", ", sep = "="),
        "\n")
  }
  cat("\n")

  print.data.frame(x[, -(1:2)], digits = digits)
}
###################################################################
qb.scantwo.slice <- function(x2, chr, slice, type.scan, smooth, weight = "sqrt")
{
  ## Get grid.
  grid <- x2$map[, 1:2]
  names(grid) <- c("chr", "pos")
  
  ## slice = c(chr=, upper=TRUE, start=, end=, weight=c(0,1,2))
  if(is.null(names(slice)))
    names(slice) <- c("chr","upper","start","end","weight")[seq(length(slice))]
  if(is.na(slice["upper"]))
    slice["upper"] <- 1
  if(is.na(slice["start"]))
    slice["start"] <- 0
  if(is.na(slice["end"]))
    slice["end"] <- max(grid$pos)
  if(is.na(slice["weight"]))
    slice["weight"] <- 2
  if(slice["upper"])
    x2$lod <- t(x2$lod)
  ## Make symmetric.
  lod2 <- t(x2$lod)
  lod2[row(lod2) > col(lod2)] <- x2$lod[row(lod2) > col(lod2)]
  
  ## Set diagonal to average of off-diagonal.
  ## Could use 1-D scan for lower?
  diaglod <- lod2[row(lod2) == 1 + col(lod2)]
  nlod <- length(diaglod)
  diag(lod2) <- (diaglod[c(1, seq(nlod))] + diaglod[c(seq(nlod), nlod)]) / 2

  ## Now get desired row(s)
  is.slice <- grid$chr == slice["chr"] & grid$pos >= slice["start"] &
  grid$pos <= slice["end"]
  if(!sum(is.slice)) {
    stop(paste("slice is invalid:",
               paste(names(slice), "=", slice, collapse = ", ")))
  }

  ## And get only desired chromosomes.
  is.chr <- !is.na(match(grid$chr, chr))

  if(sum(is.slice) == 1)
    lod2 <- lod2[is.slice, is.chr]
  else {
    ## Weighted average depending on choice of weights.
    ## For now use nitertwo, most interesting.
    if(slice["weight"] == 2) {
      lod2 <- apply(lod2[is.slice, is.chr] * x2$nitertwo[is.slice, is.chr],
                    2, sum, na.rm = TRUE) /
                      apply(x2$nitertwo[is.slice, is.chr],
                            2, sum, na.rm = TRUE)
    }
    else if(slice["weight"] == 1) {
      lod2 <- apply(lod2[is.slice, is.chr], 2, weighted.mean,
                    x2$niterone[is.slice], na.rm = TRUE)
    }
    else {
      lod2 <- apply(lod2[is.slice, is.chr], 2, mean, na.rm = TRUE)
    }
  }
  
  type.slice <- type.scan[c("lower", "upper")[(1 + slice["upper"])]]
  lod2[is.na(lod2)] <- 0
  grid[[type.slice]] <- rep(NA, nrow(grid))
  grid[is.chr, type.slice] <- qb.smoothone(lod2, grid[is.chr, ], smooth,
                                    x2$niterone[is.chr], weight = weight)
  
  ## Add smooth estimate of locus on slice chromosome.
  chr.name <- paste("chr", slice["chr"], sep = ".")
  if(sum(is.slice) == 1)
    grid[[chr.name]] <- rep(grid$pos[is.slice], nrow(grid))
  else {
    tmp <- matrix(grid$pos[is.slice], sum(is.slice), sum(is.chr))
    tmp <- apply(tmp * x2$nitertwo[is.slice, is.chr], 2, sum, na.rm = TRUE) /
      apply(x2$nitertwo[is.slice, is.chr], 2, sum, na.rm = TRUE)
    tmp[is.na(tmp)] <- mean(tmp, na.rm = TRUE)
    if(any(is.na(tmp))) ## no samples here
      grid[[chr.name]] <- rep(0, nrow(grid))
    else {
      grid[[chr.name]] <- rep(NA, nrow(grid))
      grid[is.chr, chr.name] <- qb.smoothone(tmp, grid[is.chr, ], smooth,
                                             x2$niterone[is.chr],
                                             weight = weight)
    }
  }

  ## Reduce down to desired chr.
  grid <- grid[is.chr, ]
  
  ## Make grid a scanone object.
  class(grid) <- c("scanone", "data.frame")
  attr(grid, "type.scan") <- type.slice
  attr(grid, "model") <- "normal"
  grid
}
###################################################################
qb.scantwo.smooth <- function(x, chr = NULL, smooth = 3, gridtwo, ...)
{
  type.scan <- attr(x, "type.scan")
  scan <- attr(x, "scan")
  weight <- attr(x, "weight")

  ## Force getting of 2-D sampling grid as well.
  i.lower <- gridtwo[3:4,]

  gridone <- x$grid
  mainloci <- x$mainloci

  ## Subset index for selected chromosomes.
  ## Note careful handshaking below to match up chr.sub.
  if(!is.null(chr)) {
    if(!is.numeric(chr))
      stop("chr must be numeric index to chromosomes")
    tmp <- ordered(gridone$chr)
    chr.names <- levels(tmp)[match(chr, levels(tmp))]
    chr.sub <- match(gridone$chr, chr.names)
    chr.sub <- !is.na(chr.sub)

    if(!sum(chr.sub))
      stop(paste("no samples for chromosomes",
                 chr[sort(unique(chr.sub))],
                 collapse = ","))

    gridone <- gridone[gridone$chr %in% chr.names, ]
    mainloci <- mainloci[mainloci$chrom %in% chr.names, ]
  }
  else {
#    chr <- sort(unique(gridone$chr))
    chr.sub <- rep(TRUE, nrow(gridone))
  }

  nmap <- sum(chr.sub)

  ## Get indices into lod matrix.
  tmp <- unclass(make.chr.pos(gridtwo[1,], gridtwo[2, ],
                            gridone$chr, gridone$map))
  i.lower <- unclass(make.chr.pos(i.lower[1,], i.lower[2,],
                                gridone$chr, gridone$map))
  i.upper <- tmp + (i.lower - 1) * nmap
  i.lower <- i.lower + (tmp - 1) * nmap
  chr.sub2 <- !is.na(i.lower)

  ## lod matrix has upper triangle as 2-D 
  ##                lower triangle as 2-D 
  ##                diagonal       as 1-D 
  lod <- matrix(0, nmap, nmap)
  lod[i.lower[chr.sub2]] <- x$two[chr.sub2, "lower"]
  if(all(x$two[chr.sub2, "upper"] == 0)) {
    lod[i.upper[chr.sub2]] <- x$two[chr.sub2, "lower"]
    type.scan["upper"] <- type.scan["lower"]
    scan[["upper"]] <- scan[["lower"]]
  }
  else
    lod[i.upper[chr.sub2]] <- x$two[chr.sub2, "upper"]
  diag(lod) <- x$one[chr.sub]

  ## Now get number of iterations.
  nitertwo <- matrix(0, nmap, nmap)
  tmp <- c("niter","nepis")[1 + (weight == "epistasis")]
  names(tmp) <- names(weight)
  nitertwo[i.lower[chr.sub2]] <- gridtwo[tmp["lower"], chr.sub2]
  nitertwo[i.upper[chr.sub2]] <- gridtwo[tmp["upper"], chr.sub2]
  niterone <- unclass(table(qb.inter(, gridone, mainloci)))

  ## Smooth lod matrix by chromosome.
  lod <- qb.smoothtwo(gridone, nitertwo, niterone, lod, smooth, ...)

  ## Make a scantwo object.
  lst <- list(lod = lod, map = gridone, scanoneX = NULL,
              niterone = niterone, nitertwo = nitertwo)
  attr(lst, "class") <- "scantwo"
  attr(lst, "type.scan") <- type.scan
  attr(lst, "scan") <- scan
  invisible(lst)
}
###################################################################
plot.qb.scantwo <- function(x,
                            chr = NULL,
                            smooth = 3,
                            main = mains,
                            offset = offsets,
                            nodiag = all(diag(x2$lod) == 0),
                            slice = NULL,
                            show.locus = TRUE,
                            weight = c("sqrt","count","none","atten","ratten"),
                            verbose = FALSE,
                            split.chr = attr(x, "split.chr"),
                            ...)
{
  weight <- match.arg(weight)
  geno.names <- attr(x, "geno.names")

  ## Find numerical indices for chr and slice.
  chrs <- chr <- qb.find.chr(chr = chr, geno.names = geno.names)
  if(!is.null(slice)) {
    slice <- qb.find.chr(chr = slice[1], geno.names = geno.names)
    chrs <- c(chrs, slice)
  }
  
  min.iter <- attr(x, "min.iter")
  x2 <- qb.scantwo.smooth(x, chrs, smooth,
                          qb.intertwo(min.iter, x$mainloci, x$iterdiag, x$pairloci),
                          weight = weight, ...)

  pheno.name <- attr(x, "pheno.name")
  type.scan <- attr(x2, "type.scan")
  scan <- attr(x2, "scan")
  min.iter <- attr(x, "min.iter")
  
  if(verbose) {
    cat(paste(c("\nupper:","lower:"), type.scan[c("upper","lower")], "of",
              pheno.name, "for", scan[c("upper","lower")],
              collapse = "\n"),
        "\n")
    if(min.iter > 1)
      cat("Including only loci pairs with at least", min.iter, "samples.\n")
  }
  mains <- paste(type.scan[c("upper","lower")], "of",
                 scan[c("upper","lower")], collapse = " / ")
  
  if(is.null(slice)) {
    ## Rescale values if type is estimate.
    tmpfn <- function(x) {
      max(x, -x, na.rm = TRUE)
    }
    offsets <- c(lower = 0, upper = 0)
    if(type.scan["upper"] == "estimate")
      offsets["upper"] <- tmpfn(lod[row(lod) < col(lod)])
    if(type.scan["lower"] == "estimate") {
      offsets["lower"] <- tmpfn(lod[row(lod) >= col(lod)])
    }
    if(is.null(names(offset)))
      names(offset) <- names(offsets)
    if(type.scan["upper"] == "estimate") {
      lod[row(lod) < col(lod)] <- 
        1 + (lod[row(lod) < col(lod)] / offset["upper"])
    }
    if(type.scan["lower"] == "estimate") {
      lod[row(lod) >= col(lod)] <- 
        1 + (lod[row(lod) >= col(lod)] / offset["upper"])
    }
    if(any(offset > 0)) {
      if(verbose) {
        cat("NOTE: estimate rescaled to 0 = -max, 1 = 0, 2 = max:\n",
            "max for",
            paste(names(offset), c(signif(offset,3)), collapse = ", ",
                  sep = " = "),
            "\n")
      }
      mains <- paste(mains, "\nestimate rescaled by",
                     paste(names(offset), c(signif(offset,3)),
                           collapse = ", ", sep = " = "))
    }

    ## Make sure chr is ordered with geno.names.
    grid <- data.frame(chr = x2$map$chr, pos = x2$map$map)
    
    ## Add extra marker at split points.
    split.chr <- split.chr[names(split.chr) %in% geno.names[chr]]

    if(length(split.chr)) {
      xout <- qb.chrsplit(grid, x$mainloci, chr, attr(x, "niter"), geno.names, split.chr)
      geno.names <- levels(xout$chr)
      chr <- attr(xout, "unsplit")
      x2$map$chr <- ordered(xout$chr, geno.names)
    }
    else
      x2$map$chr <- ordered(geno.names[x2$map$chr], geno.names)

    ## Plot scantwo object.
    if(compareVersion(qtlversion(), "1.04-48") < 0)
      stop("old version of R/qtl: please update now")
    
    ## plot.scantwo 1.04 assumes upper triangle is add.
    ## plot.scantwo 1.03 assumes upper triangel is epis.
    tmp <- upper.tri(x2$lod)
    x2$lod[tmp] <- t(x2$lod)[tmp] - x2$lod[tmp]

    plot(x2, nodiag = nodiag, main = main,
         incl.markers = TRUE, ...)
    if(verbose)
      cat("\n")
    invisible(x2)
  }
  else { ## 1-D slice through 2-D surface
    grid <- qb.scantwo.slice(x2, chr, slice, type.scan, smooth, weight)

    ## Plot slice.
    if(var(grid[[4]]) > 0 & show.locus) {
      tmpar <- par(mfrow=c(2,1), mar=c(2.1,4.1,0.1,0.1))
      on.exit(par(tmpar))
    }
    grid$chr <- ordered(geno.names[grid$chr], geno.names)
    plot(grid, ylim = range(grid[[3]], na.rm = TRUE), ...)
    abline(h = 0, lty = 3, lwd = 2, col = "red")
    if(var(grid[[4]]) > 0 & show.locus) {
      plot(grid, lodcolumn = 2,
           ylim = range(grid[[4]], na.rm = TRUE), ...)
      rug(attr(x, "map")[[slice[1]]], 0.02, 2, quiet = TRUE)
    }
    invisible(grid)
  }
}
###################################################################
qb.smoothtwo <- function(grid, nitertwo, niterone, x, smooth,
                         offdiag = 0.5, weight = "sqrt", ...)
{
  ## Weighted average of x.
  ## Use either local weighting or attenuated weighting.
  smooth1 <- smooth
  is.atten <- weight %in% c("atten","ratten")
  smoothpair <- function(x, smooth, w, w2 = w, offdiag)
    qb.smoothcM(x, w, w2)
  if(is.atten) {
    ## Rescale smooth (default 3 changed to 10) for sqrt attenuation.
    ## Smoothing should probably depend on niter.
    ## Should smoothing be weighted by niterone?
    if(weight == "ratten")
      smooth <- smooth * 4
  }
  else {
    switch(weight,
           count = {w <- nitertwo},
           none = {w <- array(1, dim(nitertwo))},
           sqrt =, {w <- sqrt(nitertwo)})
    
    smoothpair <- function(x, smooth, w, w2 = w, offdiag)
      qb.smoothpair(x, smooth, w, offdiag)
  }
  
  if(smooth) {
    chrs <- as.character(unique(grid$chr))
    n.chr <- length(chrs)
    offdiag <- min(1, max(0, offdiag))
    if(is.atten)
      wts <- list()
    if(n.chr == 1) {
      if(is.atten)
        wts[[chrs[1]]] <- make.atten(grid[, 2], smooth, weight)
      x <- qb.smoothsame(x, smooth,
                         if(is.atten) wts[[chrs[1]]]
                         else w,
                         offdiag, smoothpair)
    }
    else {
      for(i in seq(n.chr)) {
        rows <- chrs[i] == grid$chr
        if(sum(rows)) {
          ## Set up attenuation weights if used.
          if(is.atten) {
            if(is.null(wts[[chrs[i]]]))
              wts[[chrs[i]]] <- make.atten(grid[rows, 2], smooth, weight)
          }
          ## Diagonal matrices.
          x[rows,rows] <- qb.smoothsame(x[rows,rows], smooth,
                                        if(is.atten) wts[[chrs[i]]]
                                        else w[rows,rows],
                                        offdiag, smoothpair)
          
          ##*** This is in transistion from grid to wts.
          ## Now off diagonal matrices.
          if(i < n.chr) for(j in seq(i + 1, n.chr)) {
            cols <- chrs[j] == grid$chr
            if(sum(cols)) {
              if(is.atten) {
                if(is.null(wts[[chrs[j]]]))
                  wts[[chrs[j]]] <- make.atten(grid[cols, 2], smooth, weight)
              }
              ## Lower triangle matrix.
              x[rows,cols] <- smoothpair(x[rows,cols], smooth,
                                         if(is.atten) wts[[chrs[i]]]
                                         else w[rows,cols],
                                         if(is.atten) wts[[chrs[j]]]
                                         else NULL,
                                         offdiag)
              ## Upper triangle matrix.
              x[cols,rows] <- smoothpair(x[cols,rows], smooth,
                                         if(is.atten) wts[[chrs[j]]]
                                         else w[cols,rows],
                                         if(is.atten) wts[[chrs[i]]]
                                         else NULL,
                                         offdiag)
            }
          }
        }
      }
    }
    diag(x) <- qb.smoothone(diag(x), grid, smooth1, niterone, weight = weight)
  }
  x
}
qb.smoothcM <- function(x, wt1, wt2)
{
  ## Smooth by attenuating with cM distance between loci.
  t(wt1) %*% x %*% wt2
}
qb.smoothpair <- function(x, smooth, w, offdiag = 0.5)
{
  n.map <- dim(x)
  if(min(n.map) > 3) {
    nr <- n.map[1]
    nc <- n.map[2]
    o <- (x != 0)
    for(i in seq(smooth)) {
      ## Set up numerator = weighted sum of xs.
      wt <- w * x
      x <- (wt[, c(1, seq(nc - 1))] + wt[, c(seq(2, nc), nc)] +
            wt[c(1, seq(nr - 1)), ] + wt[c(seq(2, nr), nr), ])
      if(offdiag > 0)
        x <- x + offdiag *
          (wt[c(1, seq(nr - 1)), c(1, seq(nc - 1))] +
           wt[c(1, seq(nr - 1)), c(seq(2, nc), nc)] +
           wt[c(seq(2, nr), nr), c(seq(2, nc), nc)] +
           wt[c(seq(2, nr), nr), c(1, seq(nc - 1))])
      
      if(any(o))
        x[o] <- (x + 4 * (1 + offdiag) * wt)[o]
      
      ## Now get denominator = sum of weights.
      wt <- (w[, c(1, seq(nc - 1))] + w[, c(seq(2, nc), nc)] +
             w[c(1, seq(nr - 1)), ] + w[c(seq(2, nr), nr), ])
      if(offdiag > 0)
        wt <- wt + offdiag *
          (w[c(1, seq(nr - 1)), c(1, seq(nc - 1))] +
           w[c(1, seq(nr - 1)), c(seq(2, nc), nc)] +
           w[c(seq(2, nr), nr), c(seq(2, nc), nc)] +
           w[c(seq(2, nr), nr), c(1, seq(nc - 1))])
      ## Off-diagonal elements.
      if(any(o))
        wt[o] <- (wt + 4 * (1 + offdiag) * w)[o]
      x <- x / wt
      x[is.na(x)] <- 0
    }
  }
  x
}
qb.smoothsame <- function(x, smooth, w, offdiag = 0.5,
                          pairfun)
{
  ## Smooth upper and lower half of x, leaving diagonal unchanged.
  is.upper <- row(x) > col(x)
  is.lower <- row(x) < col(x)
  
  ## Make mat symmetric using mirror of lower triangle.
  tmpfn <- function(x, smooth, w, is.upper) {
    tmpfn2 <- function(x) {
      mat <- x
      mat[is.upper] <- t(x)[is.upper]
      ##*** Diagonal is not working properly, at least for atten.
      diagmat <- mat[row(mat) == 1 + col(mat)]
      ndiag <- length(diagmat)
      diag(mat) <- (diagmat[c(1, seq(ndiag))] +
                    diagmat[c(seq(ndiag), ndiag)]) / 2
      mat
    }
    pairfun(tmpfn2(x), smooth, tmpfn2(w),, offdiag)
  }
  x[is.lower] <- tmpfn(x, smooth, w, is.upper)[is.lower]
  x[is.upper] <- t(tmpfn(t(x), smooth, t(w), is.upper))[is.upper]
  x
}
