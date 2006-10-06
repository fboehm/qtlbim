#####################################################################
##
## $Id: slice.R,v 1.12.2.8 2006/10/06 15:17:25 byandell Exp $
##
##     Copyright (C) 2006 Brian S. Yandell
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
qb.sliceone <- function(qbObject, slice, epistasis = TRUE,
                        scan = c("main", "GxE", "epistasis"),
                        type = types,
                        covar = if(nfixcov) seq(nfixcov) else 0,
                        adjust.covar = NA,
                        chr = NULL,
                        sum.scan = "yes",
                        min.iter = 1,
                        aggregate = TRUE,
                        verbose = FALSE)
{
  ## 1-D slice through 2-D surface.

  qb.name <- deparse(substitute(qbObject))

  ## Restrict attention to selected chromosomes
  if(missing(slice))
    stop("must specify chromosome to slice upon")

  if(is.null(names(slice)))
    names(slice) <- c("chr","start","end")[seq(length(slice))]
  if(!is.null(chr))
    qbObject <- subset(qbObject, chr = c(chr, slice["chr"]),
                       restrict.pair = FALSE)
  else
    chr <- seq(length(qb.cross(qbObject)$geno))

  ## Set up slice vector.
  ## slice = c(chr=, upper=TRUE, start=, end=, weight=c(0,1,2))
  if(is.na(slice["start"]))
    slice["start"] <- 0
  if(is.na(slice["end"]))
    slice["end"] <- max(pull.grid(qbObject)$pos)

  ## Restrict attention to samples including slice.
  qbObject <- subset(qbObject, region = list(chr=slice["chr"],
                       start=slice["start"], end=slice["end"]))
  
  is.bc <- (qb.get(qbObject, "cross") == "bc")

  ## Determine variance components.
  pairloci <- qb.get(qbObject, "pairloci")
  if(is.null(pairloci))
    epistasis <- FALSE
  else { ## Restrict pairloci to pairs with slice.
    tmp <- pairloci$chrom1 == slice["chr"] | pairloci$chrom2 == slice["chr"]
    if(sum(tmp))
      pairloci <- pairloci[tmp,]
    else
      epistasis <- FALSE
  }

  var1 <- "add"
  if(epistasis)
    var2 <- "aa"
  else
    var2 <- character(0)
  if(!is.bc) {
    var1 <- c(var1,"dom")
    if(epistasis)
      var2 <- c(var2,"ad","da","dd")
  }

  ## Determine coefficients used in cellmean and ideal heritability.
  qb.coef <- if(is.bc)
    list(add = rep(c(-0.5,0.5), 2),
         aa = c(0.25,-0.25,-0.25,0.25))
  else
    list(add = rep(c(-1,0,1), 3), dom = rep(c(-0.5,0.5,-0.5), 3),
         aa = c(1,0,-1,0,0,0,-1,0,1),
         ad = c(0.5,0,-0.5,-0.5,0,0.5,0.5,0,-0.5),
         da = c(0.5,-0.5,0.5,0,0,0,-0.5,0.5,-0.5),
         dd = c(0.25,-0.25,0.25,-0.25,0.25,-0.25,0.25,-0.25,0.25))

  ## Number of fixed covariates
  nfixcov <- qb.get(qbObject, "nfixcov")
  intcov <- qb.get(qbObject, "intcov")
  
  ## Determine type of scan.
  types <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","2logBF","BF")
  if(missing(type) & length(grep("^n[pq].*", scan)) > 0)
    type <- "count"
  else {
    if (missing(type) & length(grep("^n.*", scan)) > 0) 
      type <- "log10"
    else
      type <- type[1]
  }
  type <- types[pmatch(tolower(type), tolower(types), nomatch = 1)]
  is.count <- any(match(type, c("count","log10","posterior","2logBF","BF"),
                        nomatch = 0)) 
  is.var <- match(type, types[1:6], nomatch = 0)
  is.effect <- is.var | type == "estimate"
  is.lod <- match(type, types[2:5], nomatch = 0)

  ## Number of individuals for phenotype.
  cross <- qb.cross(qbObject)
  pheno.name <- names(cross$pheno)[qb.get(qbObject, "pheno.col")]
  nind.pheno <- qb.nind.pheno(qbObject, pheno.name, nfixcov, cross)

  ## Genotype names.
  geno.names <- names(cross$geno)

  ## Following prior used for Bayes factors.
  bf.prior <- qb.get(qbObject, "mean.nqtl") /
    length(unlist(pull.loci(cross)))
  bf.prior <- bf.prior * bf.prior

  rm(cross)
  gc()
  
  ## Get MCMC samples.
  iterdiag <- qb.get(qbObject, "iterdiag")
  mainloci <- qb.get(qbObject, "mainloci")
  iterdiag.nqtl <- qb.nqtl(qbObject, iterdiag, mainloci)
  
  ## Find interaction pattern.
  if(verbose)
    cat("finding loci ...\n")
  inter <- qb.inter(qbObject)
  n.x <- length(levels(inter))

  ## Covariate adjustment calculations.
  if(type == "heritability")
    totvar <- rep(0, length(levels(inter)))
  if(nfixcov) {
    ## Covariate means.
    covar.means <- covar.mean(qbObject, adjust.covar,
                              verbose = verbose & (type == "estimate"))
    ## Explained covariance for heritability.
    if(type == "heritability") {
      tmp <- apply(qb.varcomp(qbObject, c("fixcov","rancov")), 1, sum)
      tmp <- unlist(tapply(tmp[match(mainloci$niter, iterdiag$niter)],
                           inter, mean))
      tmp[is.na(tmp)] <- 0
      totvar <- tmp
    }
  }

  ## Scan can be several variance components (default is all).
  ## Set up scan names for summaries?
  ## scans = all terms for plot
  if(type == "cellmean") {
    scans <- c(var1,var2)
    tmp <- c("A","H","B")[seq(3 - is.bc)]
    scan <- scan.save <- c(outer(tmp, tmp, paste, sep = ""))
    sum.scan <- "no"
  }
  else {
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
      scans <- NULL
      if(any(scan == aggregs[1]))
        scans <- var1
      if(any(scan == aggregs[2]))
        scans <- c(scans, var2)
      if(any(scan == aggregs[3]) & sum(intcov)) {
        if(length(covar.means))
          scans <- c(scans,
                     outer(var1, names(covar.means)[covar[intcov]],
                           paste, sep = "."))
      }
    }
    else {
      scans <- scan
      if(sum.scan != "no" & is.effect) {
        vars <- var1[match(scan, var1, nomatch = 0)]
        if(sum(intcov)) {
          if(length(covar.means) & length(vars)) {
            scans <- unique(c(scans,
                              outer(vars,
                                    names(covar.means)[covar[intcov]],
                                    paste, sep = ".")))
          }
        }
      }
    }
    if(sum.scan != "no" & (length(scans) == 1))
      sum.scan <- "no"
    ## The vector scans contains elements to scan,
    ## either to show directly or to combine in sum.
    ## The vector scan.save indicates which elements to save.
    scan.save <- if(aggregate) scan else scans
    scan.save <- switch(sum.scan,
                        no = scan.save,
                        yes = c(scan.save, "sum"),
                        only = "sum")
  }

  if(verbose) {
    cat("\n", type, "of", pheno.name, "for",
        paste(scan.save, collapse = ","), "\n")
    if(min.iter > 1)
      cat("Including only loci pairs with at least", min.iter, "samples.\n")
    cat("\n")
  }

  ## NOW GET SAMPLES AVERAGED AT CHR and POS.
  x <- matrix(0, n.x, length(scan.save))
  dimnames(x) <- list(NULL, scan.save)

  ## Extract environmental variance.
  if(type == "heritability" | is.lod |
     (type == "variance" & any(scan.save == "env"))) {
    if(verbose)
      cat("environmental variance ...\n")
    tmp <- unlist(tapply(iterdiag$envvar[match(mainloci$niter,
                                                 iterdiag$niter)],
                         inter, mean))
    tmp[is.na(tmp)] <- 0
    if(type == "heritability")
      totvar <- totvar + tmp
    else {
      if(is.lod)
        env <- tmp
      else if(any(scan.save == "env"))
        ## Actually want to examine env variance directly.
        x[,"env"] <- tmp
    }
  }

  ## Internal routine to compute count diagnostics.
  if(is.count)
    n.iter <- nrow(iterdiag)

  ## Get non-epistatic components: additive and dominance.
  if(type == "heritability")
    vars <- var1
  else
    vars <- var1[match(scans, var1, nomatch = 0)]
  if(length(vars)) {
    if(verbose)
      cat("non-epistatic components ...\n")
    ## Number of main effect samples per locus.
    if(is.count & any(scan.save == "main")) {
      tmp <- apply(as.matrix(mainloci[, paste("var", vars, sep = "")]), 1,
                   function(x) any(x > 0))
      tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
      tmp[is.na(tmp)] <- 0
      x[, "main"] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
    }
    else if(type == "cellmean") {
      tmp <- unlist(tapply(qb.meancomp(qbObject)[match(mainloci$niter,
                                                   iterdiag$niter),
                                             "grand.mean"],
                           inter, mean))
      tmp[is.na(tmp)] <- mean(tmp, na.rm = TRUE)
      tmp2 <- c("A","H","B")[seq(3 - is.bc)]
      tmp2 <- c(outer(tmp2, tmp2, paste, sep = ""))
      for(i in tmp2)
        x[, i] <- tmp
    }
    for(i in vars) {
      if(type == "estimate" | type == "cellmean")
        ## Parameter estimates of main effects.
        element <- i
      else
        ## Variance components.
        element <- paste("var", i, sep = "")
      
      ## Get samples for this component.
      main.val <- mainloci[[element]]

      ## Get GxE samples if any intcov selected.
      if(sum(intcov)) {
        if(any(scan.save == "GxE"))
          tmp2 <- rep(0, nrow(mainloci))

        if(sum(intcov)) {
          ## Get GxE samples.
          gbye <- qb.get(qbObject, "gbye")

          ## Loop over all covariates.
          for(j in intcov) {
            gbyej <- gbye[gbye$covar == j, ]
            if(length(gbyej)) {
              same <- match(paste(gbyej$niter, gbyej$chrom, gbyej$locus,
                                  sep = ":"),
                            paste(mainloci$niter, mainloci$chrom, mainloci$locus,
                                  sep = ":"))
              
              ## Parameter estimates of GxE fixed effects.
              if(match(j, covar, nomatch = 0) | type == "heritability") {
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
                    tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
                    tmp[is.na(tmp)] <- 0
                    x[, cname] <- qb.count(qbObject, tmp, type, n.iter,
                                           bf.prior)
                  }
                }
                else { ## is.effect
                  tmp <- unlist(tapply(tmp, inter, mean))
                  tmp[is.na(tmp)] <- 0
                  if(type == "heritability")
                    totvar <- totvar + tmp
                  if(match(cname, scan.save, nomatch = 0))
                    x[,cname] <- tmp
                  if(any(scan.save == "GxE"))
                    x[,"GxE"] <- x[,"GxE"] + tmp
                  if(sum.scan != "no")
                    x[,"sum"] <- x[,"sum"] + tmp
                }
              }
              
              if(type == "estimate" | type == "cellmean") {
                ## Offset parameter estimate by covariates.
                if(covar.means[j] != 0)
                  main.val[same] <- main.val[same] + covar.means[j] * gbyej[[i]]
              }
            }
            if(is.count & any(scan.save == "GxE")) {
              tmp2 <- unlist(tapply(tmp2, inter, sum, na.rm = TRUE))
              tmp2[is.na(tmp2)] <- 0
              x[, "GxE"] <- qb.count(qbObject, tmp2, type, n.iter, bf.prior)
            }
          }
        }
      }
      
      ## Now include the main effect components.
      if(is.count) {
        if(any(scan.save == i)) {
          ## Count times main effect is present.
          tmp <- unlist(tapply(main.val > 0, inter, sum, na.rm = TRUE))
          tmp[is.na(tmp)] <- 0
          x[, i] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
        }
      }
      else { ## is.effect
        ## Get main effect element and average.
        tmp <- unlist(tapply(main.val, inter, mean, na.rm = TRUE))
        tmp[is.na(tmp)] <- 0

        if(type == "cellmean") {
          ## Add contribution of element to cell mean.
          if(any(names(qb.coef) == i))
            x <- x + outer(tmp, qb.coef[[i]])
        }
        else {
          if(type == "heritability")
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
  
  ## Index for epistasis.
  if(epistasis) {
    epinter <- c(paste(pairloci$niter,
                       pairloci$chrom1, pairloci$locus1, sep = ":"),
                 paste(pairloci$niter,
                       pairloci$chrom2, pairloci$locus2, sep = ":"))
    tmp <- !duplicated(epinter)
    epinter <- ordered(epinter, epinter[tmp])
    ## epii identifies mainloci with epistatic pairs.
    epii <- match(epinter[tmp],
                  paste(mainloci$niter, inter, sep = ":"),
                  nomatch = 0)
  }

  ## Epistatic components.
  if(epistasis) {
    if(type == "heritability")
      vars <- var2
    else
      vars <- var2[match(scans, var2, nomatch = 0)]
    if(length(vars)) {
      if(verbose)
        cat("epistatic components ...\n")
      tmp <- rep(0, length(inter))

      ## Number of epistatic samples per locus.
      if(is.count & any(scan.save == "epistasis")) {
        tmp[epii] <- 1
        tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
        tmp[is.na(tmp)] <- 0
        x[, "epistasis"] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
      }

      for(i in vars) {
        if(type == "estimate" | type == "cellmean")
          ## Parameter estimates of epistasis.
          element <- i
        else
          ## Variance components for epistasis.
          element <- paste("var", i, sep = "")

        tmp2 <- pairloci[, element]
        if(is.count) {
          if(any(scan.save == i)) {
            ## Count times epistatic element is present.
            tmp[epii] <- unlist(tapply(rep(tmp2 > 0, 2), epinter, sum))[epii > 0]
            tmp2 <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
            tmp2[is.na(tmp2)] <- 0
            ## Compute count diagnostic.
            x[, i] <- qb.count(qbObject, tmp2, type, n.iter, bf.prior)
          }
        }
        else { ## is.effect
          ## Get epistatic element and average.
          tmp[epii] <- unlist(tapply(rep(tmp2, 2), epinter, sum))[epii > 0]
          tmp2 <- unlist(tapply(tmp, inter, mean))
          tmp2[is.na(tmp2)] <- 0
          
          if(type == "cellmean") {
            ## Add contribution of element to cell mean.
            if(any(names(qb.coef) == i))
              x <- x + outer(tmp2, qb.coef[[i]])
          }
          else {
            if(type == "heritability")
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
  }

  ## Extract counts averaged over MCMC runs.
  if(is.count) {
    if(verbose)
      cat("counts ...\n")
    if(sum.scan != "no") {
      ## Number of iterations per locus.
      tmp <- unclass(table(inter))
      x[, "sum"] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
    }
    
    if(any(scan == "nqtl")) {
      ## Number of QTL averaged over MCMC runs.
      tmp <- iterdiag.nqtl[match(mainloci$niter, iterdiag$niter)]
      tmp <- unlist(tapply(tmp[match(mainloci$niter, iterdiag$niter)],
                           inter, mean))
      x[, "nqtl"] <- if(type == "count") tmp else log10(tmp)
    }
  }
  else if(type != "cellmean") {
    if(type == "heritability") {
      for(i in scan.save) {
        x[, i] <- 100 * x[,i] / totvar
        x[is.na(x[, i]), i] <- 0
      }
    }
    else if(is.lod) {
      ## Number of parameters in QTL model averaged over MCMC runs.
      ## Not correct for covariates yet!
      npar <- rep(0, nrow(iterdiag))
      tmp <- apply(as.matrix(mainloci[ ,var1]), 1, function(x) sum(x > 0))
      tmp <- tapply(tmp, mainloci$niter, sum)
      npar[iterdiag.nqtl > 0] <- tmp
      if(epistasis) {
        tmp <- apply(as.matrix(pairloci[ ,var2]), 1, function(x) sum(x > 0))
        tmp <- tapply(tmp, pairloci$niter, sum)
        tmp2 <- match(unique(pairloci$niter), iterdiag$niter, nomatch = 0)
        npar[tmp2] <- npar[tmp2] + tmp
      }
      
      ## Residual sum of squares (RSS) averaged over MCMC runs.
      tmp <- (nind.pheno - npar - 1) * iterdiag$envvar
      rss <- unlist(tapply(tmp[match(mainloci$niter, 
                                     iterdiag$niter)], inter, mean))
      
      ## Keep npar for detection probability.
      if(type == "detection") {
        ## Probability of detection given data.
        ## Number of parameters averages over MCMC runs.
        npar <- unlist(tapply(npar[match(mainloci$niter, iterdiag$niter)],
                              inter, mean))
        npar[is.na(npar)] <- 0
      }
      ## Calculate LPD, LR or deviance.
      ## mostly correct for LPD, but does it get LPD?
      ## also see ideas in Gaffney code
      nscan <- 1
      for(i in scan.save) {
        if(i == "sum")
          nscan <- length(scans) - 1
        x[, i] <- nind.pheno * log((rss + env * nscan +
                                   nind.pheno * x[,i]) / rss)
        if(type == "LPD")
          x[, i] <- x[, i] / (2 * log(10))
        else if(type == "LR")
          x[, i] <- x[, i] / 2
        else if(type == "detection") {
          p1 <- exp(x[, i] / 2)
          detect.prior = 1 / nrow(x)
          x[, i] <- p1 * detect.prior / (1 + (p1 - 1) * detect.prior)
          x[is.na(x[, i]), i] <- 0.5
        }
        x[is.na(x[, i]), i] <- min(x[, i], na.rm = TRUE)
      }
    }
  }

  ## Add column for slice.
  pos <- rep(NA, nrow(x))
  tmp2 <- pairloci[, "locus2"]
  tmp <- pairloci[, "chrom1"] == slice["chr"]
  tmp2[tmp] <- pairloci[tmp, "locus1"]
  tmp <- rep(0, length(inter))
  tmp[epii] <- unlist(tapply(rep(tmp2, 2), epinter, mean,
                             na.rm = TRUE))[epii > 0]
  tmp2 <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
  tmp2[is.na(tmp2)] <- mean(tmp2, na.rm = TRUE)
  x <- cbind(x, slice = tmp2)

  ## Drop slice chromosome if not in chr list.
  if(is.na(match(slice["chr"], chr))) {
    tmp <- dimnames(x)
    x <- as.matrix(x[pull.grid(qbObject)$chr != slice["chr"], ])
    dimnames(x) <- list(NULL, tmp[[2]])
  }

  ## Assign attributes passed to generic plot and summary.
  attr(x, "class") <- c("qb.sliceone", "qb.scanone", "matrix")
  attr(x, "method") <- type
  attr(x, "scan") <- scan.save
  attr(x, "type") <- qbObject$cross
  attr(x, "chr") <- chr
  attr(x, "slice") <- slice
  attr(x, "min.iter") <- min.iter
  attr(x, "pheno.name") <- pheno.name
  attr(x, "geno.names") <- geno.names
  attr(x, "qb") <- qb.name
  x
}
###################################################################
summary.qb.sliceone <- function(object, chr = attr(object, "chr"), ...)
{
  attr(object, "class") <- c("qb.scanone", "matrix")
  summary(object, chr = chr, ...)
}
###################################################################
print.qb.sliceone <- function(x, ...) print(summary(x, ...))
###################################################################
plot.qb.sliceone <- function(x, ..., scan, auto.par = TRUE)
{
  slice <- attr(x, "slice")

  if(attr(x, "method") != "cellmean") {
    if(missing(scan))
      scan <- dimnames(x)[[2]][-ncol(x)]
    return(invisible(plot.qb.scanone(x, ..., scan = scan)))
  }
  if(missing(scan)) {
    qb <- get(attr(x, "qb"))
    is.bc <- (qb.get(qb, "cross") == "bc")
    if(is.bc) {
      if(auto.par) {
        tmpar <- par(mfcol=c(2,1), mar = c(4.1,4.1,3.1,0.1))
        on.exit(par(tmpar))
      }
      plot.qb.scanone(x, ..., scan = c("AA","HA"),
                    col = c("blue","purple"))
      plot.qb.scanone(x, ..., scan = c("AH","HH"),
                    col = c("blue","purple"))
    }
    else {
      if(auto.par) {
        tmpar <- par(mfcol=c(3,1), mar = c(4.1,4.1,3.1,0.1))
        on.exit(par(tmpar))
      }
      plot.qb.scanone(x, ..., scan = c("AA","HA","BA"),
                    col = c("blue","purple","red"))
      plot.qb.scanone(x, ..., scan = c("AH","HH","BH"),
                    col = c("blue","purple","red"))
      plot.qb.scanone(x, ..., scan = c("AB","HB","BB"),
                    col = c("blue","purple","red"))
    }
  }
  else {
    if(length(scan) == 1 & scan[1] == "slice")
      attr(x, "method") <- paste("chr", slice["chr"], "(cM)")
    plot.qb.scanone(x, scan = scan, ...)
  }
  invisible()
}
######################################################
qb.slicetwo <- function(qbObject, chr, pos, type = "2logBF", width = 10)
{
  qb.name <- deparse(substitute(qbObject))
  chr <- unlist(chr)
  pos <- unlist(pos)
  names(chr) <- names(pos) <- NULL

  slice <- list()
  for(i in 1:2) {
    ## Slice of objective function.
    ii <- paste("obj", i, sep = "")
    slice[[ii]] <-
      qb.sliceone(qbObject, chr=chr[i],
               slice=c(chr=chr[3-i],
                 start = pos[3-i] - width,
                 end = pos[3-i] + width),
               sum.scan = "no",
               type = type, scan = "epistasis")
    attr(slice[[ii]], "qb") <- qb.name

    ## Slice of Cockerham parameter estimates.
    ii <- paste("est", i, sep = "")
    slice[[ii]] <-
      qb.sliceone(qbObject, chr=chr[i],
               slice=c(chr=chr[3-i],
                 start = pos[3-i] - width,
                 end = pos[3-i] + width),
               type = "estimate", scan = "epistasis")
    attr(slice[[ii]], "qb") <- qb.name

    ## Slice of Cell means.
    ii <- paste("mean", i, sep = "")
    slice[[ii]] <-
      qb.sliceone(qbObject, chr=chr[i],
               slice=c(chr=chr[3-i],
                 start = pos[3-i] - width,
                 end = pos[3-i] + width),
               type = "cellmean")
    attr(slice[[ii]], "qb") <- qb.name
  }
  class(slice) <- c("qb.slicetwo", "list")
  attr(slice, "cross") <- qb.get(qbObject, "cross.name")
  attr(slice, "pheno.col") <- qb.get(qbObject, "pheno.col")
  attr(slice, "type") <- type
  attr(slice, "chr") <- chr
  attr(slice, "pos") <- pos
  slice
}
######################################################
summary.qb.slicetwo <- function(object, ...)
{
  out <- list(rbind(summary(object$obj1),
                    summary(object$obj2)))
  names(out) <- attr(object, "type")
  out$cellmean <- list(rbind(summary(object$mean1),
                             summary(object$mean2)))
  out$estimate <- list(rbind(summary(object$est1),
                             summary(object$est2)))
  out
}
######################################################
print.qb.slicetwo <- function(x, ...) print(summary(x, ...))
######################################################
plot.qb.slicetwo <- function(x, byrow = TRUE,
                             figs = fig.options,
                             auto.par = TRUE,
                             ...)
{
  chr <- attr(x, "chr")
  pos <- attr(x, "pos")

  cross <- get(attr(x, "cross"))
  markers <- find.marker(cross, names(cross$geno)[chr], pos)

  is.bc <- class(cross)[1] == "bc"
  if(is.bc) {
    scans <- c("AA","AH","HA","HH")
    cols <- c("blue","purple","green","red")
  }
  else {
    scans <- list(A = c("AA","HA","BA"),
                  H = c("AH","HH","BH"),
                  B = c("AB","HB","BB"))
    cols <- c("blue","purple","red")
  }
  fig.options <- c("profile", "effects", "cellmean", "effectplot")
  is.profile <- any(pmatch(tolower(figs), "profile", nomatch = 0))
  is.effects <- any(pmatch(tolower(figs), "effects", nomatch = 0))
  is.cellmean <- any(pmatch(tolower(figs), "cellmean", nomatch = 0))
  is.effectplot <- any(pmatch(tolower(figs), "effectplot", nomatch = 0))
  num.plots <-
    is.profile + is.effects + (3 - 2 * is.bc) * is.cellmean + is.effectplot

  if(auto.par) {
    if(byrow)
      tmpar <- par(mfrow=c(2, num.plots), mar=c(4.1,4.1,3.1,0.1))
    else
      tmpar <- par(mfcol=c(num.plots, 2), mar=c(4.1,4.1,3.1,0.1))
    on.exit(par(tmpar))
  }
  
  for(i in 1:2) {
    if(is.profile) {
      ## Slice of objective function.
      ii <- paste("obj", i, sep = "")
      plot(x[[ii]], scan = "epistasis",
           main = paste(attr(x[[ii]], "method"),
             "\n\nchr", chr[i], "by", chr[3-i]))
      abline(v = pos[i], lty = 2, col = "red")
    }
    if(is.effects) {
      ## Slice of Cockerham parameter estimates.
      ii <- paste("est", i, sep = "")
      plot(x[[ii]], scan = "epistasis",
           main = paste(attr(x[[ii]], "method"),
             "\n\nchr", chr[i], "by", chr[3-i]))
      abline(v = pos[i], lty = 2, col = "red")
    }
    if(is.cellmean) {
      ## Slice of Cell means.
      ii <- paste("mean", i, sep = "")
      if(is.bc) {
        plot(x[[ii]], scan = scans[c(1,i+1,4-i,4)],
             col = cols[c(1,i+1,4-i,4)],
             main = paste(attr(x[[ii]], "method"),
               "\n\nchr", chr[i], "by", chr[3-i]))
        abline(v = pos[i], lty = 2, col = "red")
      }
      else {
        for(j in c("A","H","B")) {
          plot(x[[ii]],
               scan = scans[[j]], col = cols,
               main = paste(attr(x[[ii]], "method"),
                 "\n\nchr", chr[i], "by", chr[3-i]))
          abline(v = pos[i], lty = 2, col = "red")
        }
      }
    }
    if(is.effectplot) {
      effectplot(cross, attr(x, "pheno.col"),
                 mname1 = markers[i], mname2 = markers[3-i],
                 main = paste("interaction\n\nchr",
                   paste(chr[c(i,3-i)], collapse = ", ")))
    }
  }
  invisible()
}
