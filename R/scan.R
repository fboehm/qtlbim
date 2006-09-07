#####################################################################
##
## $Id: scan.R,v 1.11.2.3 2006/09/06 01:50:22 byandell Exp $
##
##     Copyright (C) 2005 Brian S. Yandell
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
pull.loci <- function(cross, step = attr(cross$geno[[1]]$prob, "step"))
{
  if(is.null(step)) {
    warning("First running qb.genoprob with default step size",
            call. = FALSE, immediate. = TRUE)
    cross <- qb.genoprob(cross)
    step <- attr(cross$geno[[1]]$prob, "step")
  }

  ## Assume no off.end and stepwidth is variable if no prob provided.
  if(is.null(cross$geno[[1]]$prob)) {
    off.end <- 0
    stepwidth <- "variable"
  }
  else {
    off.end <- attr(cross$geno[[1]]$prob, "off.end")
    if(is.null(off.end))
      off.end <- 0
    stepwidth <- attr(cross$geno[[1]]$prob, "stepwidth")
    if(is.null(stepwidth))
      stepwidth <- "variable"
  }
  
  tmpfn <- function(x, step, off.end, stepwidth) {
    create.map(x$map, step, off.end, stepwidth)
  }
  loci <- lapply(cross$geno, tmpfn, step, off.end, stepwidth)
  class(loci) <- "map"
  loci
}
##############################################################################
pull.grid <- function (qbObject, offset = FALSE, spacing = FALSE,
                       mask.region = TRUE) 
{
  cross <- qb.cross(qbObject)
  step <- qb.get(qbObject, "step")
  if(is.null(step))
    step <- attr(cross$geno[[1]]$prob, "step")
  grid.map <- pull.loci(cross, step)

  ## Subset by region.
  cross.map <- pull.map(cross)
  if(mask.region) {
    region <- qb.get(qbObject, "subset")$region
    for(i in region$chr) {
      grid.map[[i]] <-
        grid.map[[i]][grid.map[[i]] >= region$start[i] - 0.1 &
                      grid.map[[i]] <= region$end[i] + 0.1]
    }
  }
  ## Construct map position with optional offset from 0 start.
  pos <- unlist(grid.map)
  len <- sapply(grid.map, length)
  if(!offset) {
    m <- sapply(cross.map, function(x) x[1])
    pos <- pos - rep(m, len)
  }

  ## Construct grid object with chr as first column.
  grid <- data.frame(chr = rep(seq(grid.map), len),
    row.names = paste("c", names(pos), sep = ""))

  if(spacing) {
    ## If spacing, add columns for map (=pos), eq.spacing, xchr.
    ## This is used only to create scantwo object in qb.scantwo, plot.qb.scantwo.
    grid$map <- pos
    nmap <- length(pos)
    grid$eq.spacing <- unlist(lapply(grid.map, function(x) {
      lx <- length(x)
      if(lx) {
        ## kludge to determine equal spacing
        d <- diff(x)
        tbl <- table(d)
        maxtbl <- max(tbl)
        dmode <- as.numeric(names(tbl)[tbl == maxtbl])
        if(maxtbl * 2 > lx)
          c(1, d == dmode)
        else
          rep(0, lx)
      }
      else
        integer()
    }))
    xclass <- sapply(cross.map, attr, "class")
    grid$xchr <- rep(xclass == "X", len)
  }
  else {
    ## Otherwise second column is pos.
    grid$pos <- pos
  }
  grid
}
##############################################################################
qb.nqtl <- function(qbObject,
                     iterdiag = qb.get(qbObject, "iterdiag"),
                     mainloci = qb.get(qbObject, "mainloci"))
{
  ## Fix nqtl in samples:
  ## iterdiag[, "nqtl"] may be wrong due to subsetting earlier.
  iterdiag.nqtl <- rep(0,nrow(iterdiag))
  tmp <- table(mainloci[, "niter"])
  iterdiag.nqtl[match(names(tmp), iterdiag[, "niter"])] <- tmp
  iterdiag.nqtl
}
##############################################################################
qb.inter <- function(qbObject, x = pull.grid(qbObject, offset = TRUE))
{
  ## Create identifier of chrom.locus from mainloci into pseudomarker grid.
  mainloci <- qb.get(qbObject, "mainloci")
  inter <- ordered(paste(mainloci[, "chrom"], mainloci[, "locus"], sep = ":"),
                   paste(x[, 1], x[, 2], sep = ":"))
  tmp <- is.na(inter)
  if(any(tmp)) {
    stop(paste("qb.scanone mismatch with grid:\n", sum(tmp),
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
#    if(any(is.na(match(names(threshold), dimnames(out)[[2]]))))
#      stop(paste("threshold labels should be from",
#                 dimnames(out)[[2]][-seq(pos)], collapse = ", "))
    
    use <- threshold >= 0
    if(any(use)) {
      if(length(threshold[use]) == 1)
        keep <- out[, names(threshold[use])] >= threshold[use]
      else
        keep <- apply(out[, names(threshold[use])], 1,
                      function(x) any(x >= threshold[use]))
    }
    else
      keep <- rep(FALSE, nrow(out))
    use <- threshold <= 0
    if(any(use)) {
      if(length(threshold[use]) == 1) {
        maxout <- max(out[, names(threshold[use])])
        keep <- keep | out[, names(threshold[use])] >=
          maxout + threshold[use]
      }
      else {
        maxout <- apply(out[, names(threshold[use])], 2, max)
        keep <- keep | apply(out[, names(threshold[use])],
                             1,
                             function(x)
                             any(x >= maxout + threshold[use]))
      }
    }
    if(sum(keep) > 1)
      out[keep,]
    else if(sum(keep) == 1) {
      dim.out <- dimnames(out)
      dim.out[[1]] <- dim.out[[1]][keep]
      matrix(out[keep,], 1, dimnames = dim.out)
    }
    else
      NULL
  }
  else
    out
}
##############################################################################
qb.count <- function(qbObject, stat, type, n.iter, bf.prior)
{
  if(type == "log10")
    stat <- log10(1 + stat)
  else if(type != "count" & type != "nqtl") {
    ## Else type is posterior or BF.
    stat <- (1 + stat) / (2 + n.iter)
    if(type == "logposterior") {
      stat <- log10(stat)
    }
    else {
      if(match(type, c("2logBF","BF"), nomatch = 0)) {
        stat <- stat * (1 - bf.prior) / ((1 - stat) * bf.prior)
        if(type == "2logBF")
          stat <- 2 * log(pmax(1, stat))
      }
    }
  }
  stat
}
##############################################################################
qb.scanone <- function(qbObject, epistasis = TRUE,
                     scan = c("main", "GxE", "epistasis"),
                     type = types,
                     covar = if(nfixcov) seq(nfixcov) else 0,
                     chr = NULL,
                     sum.scan = "yes",
                     min.iter = 1,
                     aggregate = TRUE,
                     half = FALSE,
                     verbose = FALSE)
{
  ## WARNING: Check covariates for npar and rss computations.
  ## Rethink qb.scan for nqtl. [Could use count already in place to do this
  ## readily, but not worth it for this freeze.]

  qb.name <- deparse(substitute(qbObject))
  
  is.bc <- (qb.get(qbObject, "cross") == "bc")

  ## Determine variance components.
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
    list(add = c(-0.5,0.5))
  else
    list(add = c(-1,0,1), dom = c(-0.5,0.5,-0.5))

  ## Number of fixed covariates
  nfixcov <- qb.get(qbObject, "nfixcov")
  intcov <- qb.get(qbObject, "intcov")
  
  ## Determine type of scan.
  types <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","logposterior","2logBF","BF","nqtl")
  type <- types[pmatch(tolower(type), tolower(types), nomatch = 2)][1]

  is.count <- any(match(type,
                        c("count", "log10", "posterior", "logposterior",
                          "2logBF", "BF", "nqtl"),
                        nomatch = 0)) 
  is.var <- match(type, types[1:6], nomatch = 0)
  is.effect <- is.var | type == "estimate"
  is.lod <- match(type, types[2:5], nomatch = 0)

  ## Number of individuals for phenotype.
  cross <- qb.cross(qbObject)
  pheno.name <- names(cross$pheno)[qb.get(qbObject, "pheno.col")]
  nind.pheno <- sum(!is.na(cross$pheno[[qb.get(qbObject, "pheno.col")]]))
  geno.names <- names(cross$geno)

  ## Following prior used for Bayes factors.
  bf.prior <- qb.get(qbObject, "mean.nqtl") / length(unlist(pull.loci(cross)))

  rm(cross)
  gc()
  
  if(!is.null(chr))
    qbObject <- subset(qbObject, chr = chr, restrict.pair = FALSE)
  
  ## Get MCMC samples.
  iterdiag <- qb.get(qbObject, "iterdiag")
  mainloci <- qb.get(qbObject, "mainloci")
  iterdiag.nqtl <- qb.nqtl(qbObject, iterdiag, mainloci)
  pairloci <- qb.get(qbObject, "pairloci")
  if(is.null(pairloci))
    epistasis <- FALSE
  
  ## Find interaction pattern.
  if(verbose)
    cat("finding loci ...\n")
  inter <- qb.inter(qbObject)

  ## Covariate adjustment calculations.
  if(type == "heritability")
    totvar <- rep(0, length(levels(inter)))
  if(nfixcov) {
    ## Covariate means.
    covar.means <- covar.mean(qbObject,
                              verbose = verbose & (type == "estimate"))
    ## Explained covariance for heritability.
    if(type == "heritability") {
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
  if(type == "cellmean") {
    scans <- var1
    scan <- scan.save <- c("A","H","B")[seq(3 - is.bc)]
    sum.scan <- "no"
  }
  else {
    if(type == "estimate" & sum.scan == "yes")
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
      if(any(scan == aggregs[3]) & sum(intcov)) {
        if(length(covar.means)) {
          tmp <- seq(nfixcov)[intcov]
          tmp <- names(covar.means)[covar[match(tmp, covar, nomatch = 0)]]
          if(length(tmp)) {
             scans <- c(scans, outer(var1, tmp, paste, sep = "."))
             scans.save <- c(scans.save, outer(var1, tmp, paste, sep = "."))
           }
        }
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
    if(sum.scan != "no" & (length(scans) == 1))
      sum.scan <- "no"
    ## The vector scans contains elements to scan,
    ## either to show directly or to combine in sum.
    ## The vector scan.save indicates which elements to save.
    scan.save <- if(aggregate) scan else scans.save
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
  x <- matrix(0, length(levels(inter)), length(scan.save))
  dimnames(x) <- list(NULL, scan.save)

  ## Extract environmental variance.
  if(type == "heritability" | is.lod |
     (type == "variance" & any(scan.save == "env"))) {
    if(verbose)
      cat("environmental variance ...\n")
    tmp <- unlist(tapply(iterdiag[match(mainloci[, "niter"], iterdiag[, "niter"]),
                                  "envvar"],
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

  ## Need n.iter for counts.
  if(is.count)
    n.iter <- nrow(iterdiag)

  if(type == "nqtl")
    nqtl.main <- paste(mainloci[, "niter"], mainloci[, "chrom"], sep = ":")

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
      if(type == "nqtl") {
        tmp <- tapply(tmp, nqtl.main, sum)[nqtl.main]
        tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
      }
      else
        tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
      tmp[is.na(tmp)] <- 0
      x[, "main"] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
    }
    else if(type == "cellmean") {
      tmp <- unlist(tapply(qb.meancomp(qbObject)[match(mainloci[, "niter"],
                                                   iterdiag[, "niter"]),
                                             "grand.mean"],
                           inter, mean))
      tmp[is.na(tmp)] <- mean(tmp, na.rm = TRUE)
      x[, "A"] <- x[, "H"] <- tmp
      if(!is.bc)
        x[, "B"] <- tmp
    }
    for(i in vars) {
      if(type == "estimate" | type == "cellmean")
        ## Parameter estimates of main effects.
        element <- i
      else
        ## Variance components.
        element <- paste("var", i, sep = "")
      
      ## Get samples for this component.
      main.val <- mainloci[, element]

      ## Get GxE samples if any intcov selected.
      if(sum(intcov)) {
        if(any(scan.save == "GxE"))
          tmp2 <- rep(0, nrow(mainloci))

        if(sum(intcov)) {
          ## Get GxE samples.
          gbye <- qb.get(qbObject, "gbye")

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
                    if(type == "nqtl") {
                      tmp <- tapply(tmp, nqtl.main, sum)[nqtl.main]
                      tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
                    }
                    else
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
              
              if(type == "estimate") {
                ## Offset parameter estimate by covariates.
                if(covar.means[j] != 0)
                main.val[same] <- main.val[same] + covar.means[j] * gbyej[[i]]
              }
            }
            if(is.count & any(scan.save == "GxE")) {
              if(type == "nqtl") {
                tmp2 <- tapply(tmp2, nqtl.main, sum)[nqtl.main]
                tmp2 <- unlist(tapply(tmp2, inter, mean, na.rm = TRUE))
              }
              else
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
          if(type == "nqtl") {
            tmp <- tapply(main.val > 0, nqtl.main, sum)[nqtl.main]
            tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
          }
          else
            tmp <- unlist(tapply(main.val > 0, inter, sum, na.rm = TRUE))
          tmp[is.na(tmp)] <- 0
          x[, i] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
        }
      }
      else { ## is.effect
        tmp <- unlist(tapply(main.val, inter, mean, na.rm = TRUE))
        tmp[is.na(tmp)] <- 0
        if(type == "cellmean") {
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
    epinter <- c(paste(pairloci[, "niter"],
                       pairloci[, "chrom1"], pairloci[, "locus1"], sep = ":"),
                 paste(pairloci[, "niter"],
                       pairloci[, "chrom2"], pairloci[, "locus2"], sep = ":"))
    if(type == "nqtl") {
      nqtl.pair <- c(paste(pairloci[, "niter"], pairloci[, "chrom1"], sep = ":"),
                     paste(pairloci[, "niter"], pairloci[, "chrom2"], sep = ":"))
    }
    tmp <- !duplicated(epinter)
    epinter <- ordered(epinter, epinter[tmp])
    ## epii identifies mainloci with epistatic pairs.
    epii <- match(epinter[tmp],
                  paste(mainloci[, "niter"], inter, sep = ":"),
                  nomatch = 0)
  }

  ## Epistatic components.
  if(epistasis & type != "cellmean") {
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
        if(type == "nqtl") {
          tmp2 <- table(nqtl.pair)[nqtl.pair]
          tmp[epii] <- unlist(tapply(tmp2, epinter, mean))[epii > 0]
          tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
        }
        else {
          tmp[epii] <- 1
          tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
        }
        tmp[is.na(tmp)] <- 0
        x[, "epistasis"] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
      }

      for(i in vars) {
        if(type == "estimate")
          ## Parameter estimates of epistasis.
          element <- i
        else
          ## Variance components for epistasis.
          element <- paste("var", i, sep = "")

        tmp2 <- pairloci[, element]
        if(is.count) {
          if(any(scan.save == i)) {
            ## Count times epistatic element is present.
            if(type == "nqtl") {
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
            x[, i] <- qb.count(qbObject, tmp2, type, n.iter, bf.prior)
            
          }
        }
        else { ## is.effect
          ## Get epistatic element and average.
          tmp[epii] <- unlist(tapply(rep(tmp2, 2), epinter, sum))[epii > 0]
          tmp2 <- unlist(tapply(tmp, inter, mean))
          
          ## Reduce epistatic value by half.
          if(half)
            tmp2 <- tmp2 / 2
          tmp2[is.na(tmp2)] <- 0
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

  ## Extract counts averaged over MCMC runs.
  if(is.count) {
    if(verbose)
      cat("counts ...\n")
    if(sum.scan != "no") {
      ## Number of iterations per locus.
      if(type == "nqtl") {
        tmp <- table(nqtl.main)[nqtl.main]
        tmp <- unlist(tapply(tmp, inter, mean, na.rm = TRUE))
        tmp[is.na(tmp)] <- 0
      }
      else
        tmp <- unclass(table(inter))
      x[, "sum"] <- qb.count(qbObject, tmp, type, n.iter, bf.prior)
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
      tmp <- apply(as.matrix(mainloci[, var1]), 1, function(x) sum(x > 0))
      tmp <- tapply(tmp, mainloci[, "niter"], sum)
      npar[iterdiag.nqtl > 0] <- tmp
      if(epistasis) {
        tmp <- apply(as.matrix(pairloci[, var2]), 1, function(x) sum(x > 0))
        tmp <- tapply(tmp, pairloci[, "niter"], sum)
        tmp2 <- match(unique(pairloci[, "niter"]), iterdiag[, "niter"], nomatch = 0)
        npar[tmp2] <- npar[tmp2] + tmp
      }
      
      ## Residual sum of squares (RSS) averaged over MCMC runs.
      tmp <- (nind.pheno - npar - 1) * iterdiag[, "envvar"]
      rss <- unlist(tapply(tmp[match(mainloci[, "niter"], iterdiag[, "niter"])],
                           inter, mean))
      
      ## Keep npar for detection probability.
      if(type == "detection") {
        ## Probability of detection given data.
        ## Number of parameters averages over MCMC runs.
        npar <- unlist(tapply(npar[match(mainloci[, "niter"], iterdiag[, "niter"])],
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

  ## Assign attributes passed to plot.qb.scanone.
  attr(x, "class") <- c("qb.scanone", "matrix")
  attr(x, "method") <- type
  attr(x, "scan") <- scan.save
  attr(x, "type") <- qbObject$cross
  attr(x, "chr") <- chr
  attr(x, "min.iter") <- min.iter
  attr(x, "pheno.name") <- pheno.name
  attr(x, "geno.names") <- geno.names
  attr(x, "qb") <- qb.name
  x
}
###################################################################
summary.qb.scanone <- function(object,
                             chr, ## Must be integer for now.
                             threshold = 0,
                             sort = "no",
                             digits = 3,
                             smooth = 3,
                             n.qtl = 0.05,
                             min.iter = attr(object, "min.iter"),
                             verbose = FALSE,
                             ...)
{
  qbObject <- get(attr(object, "qb"))
  chr.qb <- attr(object, "chr")
  if(!is.null(chr.qb))
    qbObject <- subset(qbObject, chr = chr.qb)
  type <- attr(object, "method")
  
  pheno.name <- attr(object, "pheno.name")

  scan <- attr(object, "scan")
  if(verbose) {
    cat(type, "of", pheno.name, "for",
        paste(scan, collapse = ","), "\n")
    if(min.iter > 1)
      cat("Including only loci pairs with at least", min.iter, "samples.\n")
    if(any(threshold != 0)) {
      cat("Thresholds:",
          paste(names(threshold), c(threshold), collapse = ", ", sep = "="),
          "\n")
    }
    cat("\n")
  }

  ## Get chr and pos.
  x <- pull.grid(qbObject, offset = TRUE)

  ## Get interaction pattern.
  inter <- qb.inter(qbObject, x)

  ## Get number of samples per pos.
  niter <- unclass(table(inter))
  n.iter <- nrow(qb.get(qbObject, "iterdiag"))

  mainloci <- qb.get(qbObject, "mainloci")

  main.scan <- c("main","add","dom")
  ## Epistasis counter.
  epi.scan <- c("epistasis","aa","ad","da","dd")
  epistasis <- any(match(scan, epi.scan, nomatch = 0))
  if(epistasis) {
    qb.get.epis <- function(qbObject, inter, mainloci)
    {
      pairloci <- qb.get(qbObject, "pairloci")
      epinter <- c(paste(pairloci[, "niter"], pairloci[, "chrom1"],
                         pairloci[, "locus1"], sep = ":"),
                   paste(pairloci[, "niter"], pairloci[, "chrom2"],
                         pairloci[, "locus2"], sep = ":"))
      rm(pairloci)
      gc()
      tmp <- !duplicated(epinter)
      epinter <- ordered(epinter, epinter[tmp])
      ## epii identifies mainloci with epistatic pairs.
      epii <- match(epinter[tmp],
                    paste(mainloci[, "niter"], inter, sep = ":"),
                    nomatch = 0)
      tmp <- rep(0, length(inter))
      tmp[epii] <- 1
      tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
      tmp[is.na(tmp)] <- 0
      tmp
    }    
    nepis <- qb.get.epis(qbObject, inter, mainloci)
  }
  qb.get.main <- function(qbObject, inter, mainloci)
  {
    vars <- c("varadd","vardom")
    vars <- vars[!is.na(match(vars, names(mainloci)))]
    tmp <- apply(as.matrix(mainloci[, vars]), 1,
                   function(x) any(x > 0))
    tmp <- unlist(tapply(tmp, inter, sum, na.rm = TRUE))
    tmp[is.na(tmp)] <- 0
    tmp
  }
  x.main <- qb.get.main(qbObject, inter, mainloci)
  rm(mainloci)
  gc()
  
  ## Get chromosome names and set up matrix.
  chrs = unique(x$chr)
  if(!missing(chr))
    chrs <- chrs[match(chr, chrs, nomatch = 0)]
  
  out <- matrix(0, length(chrs), 4 + ncol(object) + epistasis)
  object.names <- dimnames(object)[[2]]
  dimnames(out) <- list(attr(object, "geno.names")[chrs],
                        c("chr", "n.qtl", "pos", "m.pos", "e.pos"[epistasis],
                          object.names))

  out[, "chr"] <- chrs

  
  ## Center as mean for variance, estimate, cellmean.
  ## Center as mode for all other types.
  center <- "mean"
  if(is.na(match(type, c("variance","estimate","cellmean"))))
    center <- "mode"

  ## Values at maximum number of smoothed iterations.
  x.iter <- qb.smoothone(niter, x, smooth, niter)
  x.main <- qb.smoothone(x.main, x, smooth, niter)
  if(epistasis)
    x.epis <- qb.smoothone(nepis, x, smooth, nepis)

  for(i in seq(ncol(object))) {
    if(match(object.names[i], epi.scan, nomatch = 0))
      object[,i] <- qb.smoothone(object[,i], x, smooth, nepis)
    else
      object[,i] <- qb.smoothone(object[,i], x, smooth, niter)
  }
  for(i in seq(length(chrs))) {
    ii <- (x$chr == chrs[i] & niter >= min.iter)
    if(any(ii)) {
      wh <- which.max(x.iter[ii])
      out[i, "n.qtl"] <- sum(niter[ii]) / n.iter
      out[i, "pos"] <- x$pos[ii][wh]
      m.wh <- which.max(x.main[ii])
      out[i, "m.pos"] <- x$pos[ii][m.wh]
      if(epistasis) {
        e.wh <- which.max(x.epis[ii])
        out[i, "e.pos"] <- x$pos[ii][e.wh]
      }
      for(j in object.names) {
        out[i, j] <-
          if(match(j, epi.scan, nomatch = 0))
            object[ii, j][e.wh]
          else {
            if(match(j, main.scan, nomatch = 0))
              object[ii, j][m.wh]
            else
              object[ii, j][wh]
          }               
      }
#        if(center == "mean")
        ## Weighted means by chr.
#          c(apply(cbind(x$pos[ii], object[ii, ]), 2, weighted.mean, 
#                  niter[ii]), sum(niter[ii]))
#        else { ## center == "mode"
#          c(x$pos[ii][which.max(object[ii, which.pos])],
#            apply(object[ii, ], 2, max), sum(niter[ii]))
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
      if(sum(tmp) == 1) {
        outn <- dimnames(out)
        out <- matrix(out[tmp,], 1)
        dimnames(out) <- list(outn[[1]][tmp], outn[[2]])
      }
      else
        out <- out[tmp, ]
      round({if (match(sort, dimnames(out)[[2]], nomatch = 0) & nrow(out) > 1)
               out[order(- out[, sort]), ] else out
           }, digits)
    }
    else
      NULL
  }
}
###################################################################
print.qb.scanone <- function(x, ...) print(summary(x, ...))
###################################################################
plot.qb.scanone <- function(x,
                          chr = NULL,
                          smooth = 0,
                          scan = scan.plots,
                          ylim = ylims,
                          scan.name = scan.pretty,
                          col = cols,
                          main = paste(type, "of", pheno.name,
                            "for", scan.name),
                          verbose = FALSE,
                          ...)
{
  ## Now need this to pick up sum
  ## and also to get colors right
  ## and to get ylim right when type is estimate

  type <- attr(x, "method")

  ## Set up output grid.
  qbObject <- get(attr(x, "qb"))
  chr.qb <- attr(x, "chr")
  if(!is.null(chr.qb))
    qbObject <- subset(qbObject, chr = chr.qb)
  
  xout <- pull.grid(qbObject, offset = TRUE)

  ## Subset index for selected chromosomes.
  if(!is.null(chr)) {
    if(!is.numeric(chr))
      stop("chr must be numeric index to chromosomes")
    chr.sub <- match(xout$chr, chr)
    chr <- chr[sort(unique(chr.sub))]
    chr.sub <- !is.na(chr.sub)
    if(!sum(chr.sub))
      stop(paste("no samples for chromosomes", chr, collapse = ","))
    qbObject <- subset(qbObject, chr = chr)
    xout <- pull.grid(qbObject, offset = TRUE)
  }
  else
    chr.sub <- rep(TRUE, nrow(xout))
  
  ## Figure out how to organize scans.
  scan.names <- dimnames(x)[[2]]
  scan.plots <-
    if(length(scan.names) < 5)
      scan.names
    else
      c("sum","main","epistasis")
  is.sum <- match("sum", scan.names, nomatch = 0)
  ## Colors for plots (does not allow for covariates yet.
  scan.col <- function(x, allscan) {
    if(!allscan)
      return("black")
    col <- c(sum = "black",
             add = "blue", dom = "red",
             aa = "purple",
             ad = "green", da = "darkgreen",
             dd = "orange",
             main = "blue", epistasis = "purple", GxE = "darkred",
             A = "blue", H = "purple", B = "red")
    cols <- col[x]
    cols[grep(".add", x)] <- "darkblue"
    cols[grep(".dom", x)] <- "darkred"
    name.col <- names(cols)
    cols[is.na(cols)] <- "black"
    name.col[is.na(name.col)] <- x[is.na(name.col)]
    names(cols) <- name.col
    cols
  }

  ## Automate separate plots by main, epistasis, sum.
  if(any(match(c("main","epistasis"), scan, nomatch = 0)) &
     length(scan.names) > 5) {
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
    cols <- scan.col(scans, length(scans) > 1)
    scans <- scan.names[is.sum]
    if(is.main)
      scans <- c(scans, scan.main)
    if(is.epis)
      scans <- c(scans, scan.epis)
    ylims <- range( c(x[chr.sub, scans]), na.rm = TRUE)
    if(is.sum)
      plot.qb.scanone(x, chr, smooth, scan.names[is.sum], ylim, "all effects",
                    col = col, ...)
    if(is.main)
      plot.qb.scanone(x, chr, smooth, scan.main, ylim, "main effects",
                    col = col, ...)
    if(is.epis)
      plot.qb.scanone(x, chr, smooth, scan.epis, ylim, "epistatic effects",
                    col = col, ...)
    par(tmpar)
    return(invisible(col))
  }

  min.iter <- attr(x, "min.iter")
  if(min.iter > 1) {
    ## Get number of samples per pos.
    niter <- unclass(table(qb.inter(qbObject, xout)))

    ## Reduce x to loci with at least min.iter samples.
    x[chr.sub & niter < min.iter, scan] <- 0
  }

  ## Process the selected scan terms.
  ylims <- range( c(x[chr.sub, scan]), na.rm = TRUE)
      
  ## Work on pretty title.
  if(length(scan) < 4)
    scan.pretty <- paste(scan, collapse="+")
  else
    scan.pretty <- "effects"

  ## Figure out phenotype name indirectly.
  pheno.name <- attr(x, "pheno.name")

  ## Print message about plot.
  if(verbose) {
    cat("\n", attr(x, "method"), "of", pheno.name, "for",
        paste(scan, collapse = ","), "\n")
  }

  ## Set up scanone object attributes.
  class(xout) <- c("scanone", "data.frame")
  attr(xout, "method") <- type
  attr(xout, "type") <- attr(x, "type")
  attr(xout, "model") <- "normal"

  ## Find sum, if more than one scan, and color scheme.
  is.sum <- match("sum", scan, nomatch = 1)
  allscan <- length(scan) > 1

  ## Set up color scheme.
  cols <- scan.col(scan, allscan)
  names(cols) <- scan
  ## Set any missing colors to "black".
  tmp <- length(scan) - length(col)
  tmp2 <- names(col)
  if(tmp > 0) {
    if(is.null(tmp2)) {
      col <- c(col, rep("black", tmp))
      names(col) <- scan
    }
    else {
      cols[tmp2] <- col
      col <- cols
    }
  }
  else {
    if(is.null(tmp2)) {
      col <- col[seq(length(scan))]
      names(col) <- scan
    }
    else
      col <- col[match(scan, tmp2, nomatch = 0)]
  }
    
  ## Smooth over points?
  ## This should be done by chr, possibly using smooth.spline.
  ## Need to start plotting with first scan object.
  niter <- unclass(table(qb.inter(qbObject, xout)))
  xout[ ,type] <- qb.smoothone(x[chr.sub, scan[is.sum]], xout, smooth, niter)
  plot(xout, ..., ylim = ylim, main = main, col = col[is.sum])
  if(type == "log10") {
    tmp <- c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000)
    axis(4,log10(tmp),tmp)
  }
  if(type == "estimate")
    abline(h = 0, col = "grey", lty = 3, lwd = 2)

  ## All other scan objects added to plot.
  if(allscan) {
    for(varcomp in rev(scan[-is.sum])) {
      xout[, type] <- qb.smoothone(x[chr.sub, varcomp], xout, smooth, niter)
      plot(xout, ..., add = TRUE, col = col[varcomp])
    }
    if(length(col) < 5)
      mtext(paste(names(col), col, sep = "=", collapse = ", "), 1, 2,
            cex = 0.65)
  }
  invisible(col)
}
############################################################################## 
qb.smoothone <- function(x, xout, smooth, niter)
{
  if(smooth) {
    qb.smoothchr <- function(x, smooth, niter) {
      nmap <- length(x)
      if(nmap > 3) {
        o <- (x != 0)
        for(i in seq(smooth)) {
          wtlod <- niter * x
          x <- wtlod[c(1, seq(nmap - 1))] + wtlod[c(seq(2, nmap), nmap)]
          if(any(o))
            x[o] <- x[o] + 2 * wtlod[o]
          wtlod <- niter[c(1, seq(nmap - 1))] + niter[c(seq(2, nmap), nmap)]
          if(any(o))
            wtlod[o] <- wtlod[o] + 2 * niter[o]
          x <- x / wtlod
          x[is.na(x)] <- 0
        }
      }
      x
    }
    for(chr in unique(xout$chr)) {
      rows <- chr == xout$chr
      if(sum(rows))
        x[rows] <- qb.smoothchr(x[rows], smooth, niter[rows])
    }
  }
  x
}
##############################################################################
##############################################################################
## need iterdiag.nqtl as well as mainloci!
qb.indextwo <- function(qbObject,
                         iterdiag = qb.get(qbObject, "iterdiag"),
                         mainloci = qb.get(qbObject, "mainloci"),
                         nqtl = qb.nqtl(qbObject, iterdiag, mainloci))
{
  ## index of pairs of loci for each iteration
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
qb.intertwo <- function(qbObject,
                         min.iter = 1,
                         mainloci = qb.get(qbObject, "mainloci"),
                         index = qb.indextwo(qbObject, iterdiag, mainloci,
                           nqtl),
                         iterdiag = qb.get(qbObject, "iterdiag"),
                         pairloci = qb.get(qbObject, "pairloci"),
                         nqtl = qb.nqtl(qbObject, iterdiag, mainloci))
{
  gridtwo <- matrix(t(mainloci[index, c("chrom","locus")]), 4)
  inter <- paste(gridtwo[1, ], gridtwo[2, ],
                 gridtwo[3, ], gridtwo[4, ], sep = ":")
  inter <- ordered(inter, inter[!duplicated(inter)])
  
  ## Set up epistatic count, which requires several other things.
  npair <- qb.nqtl(qbObject, iterdiag, mainloci)
  npair <- npair * (npair - 1) / 2
  epi <- rep(0, length(inter))
  epi[match(paste(pairloci[, "niter"],
                  pairloci[, "chrom1"], pairloci[, "locus1"],
                  pairloci[, "chrom2"], pairloci[, "locus2"],
                  sep = ":"), 
            paste(rep(iterdiag[, "niter"], npair),
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
                     type = c(
                       upper = "heritability",
                       lower = "heritability"),
                     upper.scan = "epistasis",
                     lower.scan = "joint",
                     covar = if(nfixcov) seq(qbObject$nfixcov) else 0,
                     chr = NULL,
                     min.iter = 1,
                     verbose = FALSE)
{
  ## Need to add aggregate facilities and redo counts as in qb.scanone.

  ## Following prior used for Bayes factors.
  ## Need to do this before subsetting on chr.
  bf.prior <- qb.get(qbObject, "mean.nqtl") /
    length(unlist(pull.loci(qb.cross(qbObject))))
  
  qb.name <- deparse(substitute(qbObject))
  if(!is.null(chr))
    qbObject <- subset(qbObject, chr = chr, restrict.pair = FALSE)
  
  nfixcov <- qb.get(qbObject, "nfixcov")
  intcov <- qb.get(qbObject, "intcov")

  pairloci <- qb.get(qbObject, "pairloci")
  if(is.null(pairloci))
    epistasis <- FALSE

  ## Determine type of qb.scan.
  types <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","logposterior","BF","2logBF","nqtl")
  tmp <- names(type)
  type <- types[pmatch(tolower(type), tolower(types), nomatch = 2,
                       duplicates.ok = TRUE)]
  type <- array(type, 2)
  if(is.null(tmp))
    tmp <- c("upper","lower")
  names(type) <- tmp

  if(any(type == "cellmean"))
    stop("cellmean type not implemented: use qb.slice")
  is.count <- match(type,
                    c("count", "log10", "posterior", "logposterior",
                      "BF", "2logBF", "nqtl"),
                    nomatch = 0)
  names(is.count) <- names(type)

  is.var <- match(type, types[1:6], nomatch = 0)
  is.effect <- is.var | type == "estimate"
  is.lod <- match(type, types[2:5], nomatch = 0)
  names(is.var) <- names(is.effect) <- names(is.lod) <- names(type)
  
  ## Number of individuals for phenotype.
  cross <- qb.cross(qbObject)
  pheno.name <- names(cross$pheno)[qb.get(qbObject, "pheno.col")]
  nind.pheno <- sum(!is.na(cross$pheno[[qb.get(qbObject, "pheno.col")]]))
  map <- pull.map(cross)
  rm(cross)
  gc()
  
  ## Get MCMC samples.
  iterdiag <- qb.get(qbObject, "iterdiag")
  mainloci <- qb.get(qbObject, "mainloci")
  iterdiag.nqtl <- qb.nqtl(qbObject, iterdiag, mainloci)
  npair <- iterdiag.nqtl
  npair <- npair * (npair - 1) / 2

  ## Determine variance components.
  is.bc <- (qb.get(qbObject, "cross") == "bc")
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
  index <- qb.indextwo(qbObject, iterdiag, mainloci, iterdiag.nqtl)
  nindex <- length(index) / 2

  ## Find all pairs of loci.
  if(verbose)
    cat("finding all pairs of loci ...\n")
  tmp <- matrix(t(mainloci[index, c("chrom","locus")]), 4)
  inter <- paste(tmp[1, ], tmp[2, ],
                 tmp[3, ], tmp[4, ], sep = ":")
  inter <- ordered(inter, inter[!duplicated(inter)])

  ## Set up index for number of linked qtl.
  if(any(type == "nqtl")) {
#    nqtl.main <- paste(mainloci[index[seq(by = 2, length = nindex)], "niter"],
#                       tmp[1, ], tmp[3, ], sep = ":")
    nqtl.main <- paste(mainloci[, "niter"], mainloci[, "chrom"], sep = ":")
  }

  ## Covariate adjustment calculations.
  if(any(type == "heritability"))
    totvar <- rep(0, length(inter))
  if(nfixcov) {
    ## Covariate means.
    covar.means <- covar.mean(qbObject,
                              verbose = verbose & (any(type == "estimate")))

    ## Explained covariance for heritability.
    if(any(type == "heritability"))
      totvar <- rep(apply(qb.varcomp(qbObject, c("fixcov","rancov")),
                          1, sum),
                    npair)
  }

  ## Set up lists (elements upper, lower) of scan names.
  ## Scan can be several variance components (default is all).
  ## upper, lower = original or default scan names
  ##   ("epistasis","joint" or "main","GxE").
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
  
  ## Match key words if present.
  for(tri in c("lower","upper")) {
    keys <- c("main","epistasis","joint")
    tmp <- pmatch(scan[[tri]], keys, nomatch = 0, duplicates.ok = TRUE)
    if(any(tmp))
      scan[[tri]][tmp > 0] <- keys[tmp]
  }
  ## If no epistasis and missing upper, then set upper scan to main.
  if(!epistasis) {
    if(missing(scan) & missing(upper.scan))
      scan[["upper"]] <- "main"
  }
  ## Save scan names for later plot.
  if(any(type == "estimate")) {
    if(missing(scan) & missing(upper.scan) & type["upper"] == "estimate")
      scan$upper <- "aa"
    if(missing(scan) & missing(lower.scan) & type["lower"] == "estimate")
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
                        joint = {
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
                                c("main","sum","joint","add","dom"),
                                nomatch = 0)))
      "main"
    else
      "epistasis"
  }

  if(verbose) {
    cat(paste(c("\nupper:","lower:"), type[c("upper","lower")], "of",
              pheno.name, "for",
              sapply(scan, paste, collapse = "+")[c("upper","lower")],
              collapse = "\n"),
        "\n")
    if(min.iter > 1)
      cat("Including only loci pairs with at least", min.iter, "samples.\n")
    cat("\n")
  }

  ## Extract environmental variance.
  if(any(type == "heritability") | any(is.lod)) {
    if(verbose)
      cat("environmental variance ...\n")
    tmp <- rep(iterdiag[, "envvar"], npair)
    if(any(type == "heritability"))
      totvar <- totvar + tmp
    else if(any(is.lod))
      env <- unlist(tapply(tmp, inter, mean))
  }

  var.elem <- function(type, vari) {
    ifelse(type == "estimate", vari, paste("var", vari, sep = ""))
  }

  accum <- matrix(0, nindex, 2)
  dimnames(accum) <- list(NULL, c("lower","upper"))
  
  ## Number of main effect samples per locus.
  is.joint <- (is.count &
               sapply(scan, function(x) match("joint", x, nomatch = 0)))
  names(is.joint) <- names(is.count)
  if(any(is.joint)) {
    for(tri in names(is.joint)[is.joint])
      accum[, tri] <- 1
  }
  
  ## Get non-epistatic components: additive and dominance.
  if(any(type == "heritability"))
    vars <- var1
  else
    vars <- var1[match(unique(unlist(scan)), var1, nomatch = 0)]
  if(verbose & length(vars))
    cat("non-epistatic components ...\n")
  for(vari in vars) {
    if(any(type == "estimate"))
      main.val <- mainloci[, vari]
    
    ## Loop over all interacting covariates.
    if(sum(intcov)) {
      ## Get GxE samples.
      gbye <- qb.get(qbObject, "gbye")
  
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
               type[tri] == "heritability") {
              if(tri == "lower" | type[tri] != type["lower"]) {
                cname <- paste(vari, paste(names(covar.means)[covj],
                                           sep = "."))
                cov.val[same] <- gbyej[[var.elem(type[tri], vari)]]
                tmp <- apply(array(cov.val[index], c(2,nindex)), 2, sum)
              }
              if(type[tri] == "heritability" & !done.her) {
                totvar <- totvar + tmp
                done.her <- TRUE
              }
                if(match(cname, scan[[tri]], nomatch = 0))
                  accum[, tri] <- accum[, tri] + tmp
            }
            
            ## Offset parameter estimate by covariates.
            if(type[tri] == "estimate" & covar.means[covj] != 0 &
               !done.est) {
              main.val[same] <- main.val[same] +
                covar.means[covj] * gbyej[[vari]]
              done.est <- TRUE
            }
          }
          }
      }
    }

    ## Now include the component.
    done.her <- FALSE
    for(tri in c("lower","upper")) {
      if(tri == "lower" | type[tri] != type["lower"]) {
        if(type[tri] == "estimate")
          tmp <- main.val
        else ## variance components and counts
          tmp <- mainloci[, paste("var", vari, sep = "")]
        tmp <- apply(array(tmp[index], c(2,nindex)), 2, sum)
      }
      if(type[tri] == "heritability" & !done.her) {
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
    if(any(type == "heritability"))
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
      if(any(type == "nqtl")) {
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
          element <- var.elem(type[tri], vari)
          if(type[tri] == "heritability" & !done.her) {
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
      if(type[tri] == "nqtl") {
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
        qb.count(qbObject, tmp, type[tri], n.iter, bf.prior)
    }
    else { ## is.effect
      tmp <- unlist(tapply(accum[, tri], inter, mean))
      tmp[is.na(tmp)] <- 0
      accum[seq(n.inter), tri] <- tmp
    }
  }
  accum <- accum[seq(n.inter), ]
  
  if(any(type == "heritability")) {
    totvar <- unlist(tapply(totvar, inter, mean))
    totvar[is.na(totvar)] <- 0
  }

  
  ## Compute heritability.
  if(any(type == "heritability") & verbose)
    cat("heritability ...\n")
  for(tri in c("lower","upper")) {
    if(type[tri] == "heritability") {
      accum[, tri] <- 100 * accum[, tri] / totvar
      accum[, tri][is.na(accum[, tri])] <- 0
    }
  }
  if(any(is.lod)) {
    if(verbose)
      cat("LPD or other diagnostics ...\n")
    ## Number of parameters in QTL model averaged over MCMC runs.
    tmp <- apply(as.matrix(mainloci[ ,var1]), 1, function(x) sum(x > 0))
    tmp <- tapply(tmp, mainloci[, "niter"], sum)
    npar <- rep(0, nrow(iterdiag))
    npar[iterdiag.nqtl > 0] <- tmp
    if(epistasis) {
      tmp <- apply(as.matrix(pairloci[, var2]), 1, function(x) sum(x > 0))
      tmp <- tapply(tmp, pairloci[, "niter"], sum)
      tmp2 <- match(unique(pairloci[, "niter"]), iterdiag[, "niter"], nomatch = 0)
      npar[tmp2] <- npar[tmp2] + tmp
    }
    ## Residual sum of squares.
    rss <- (nind.pheno - npar - 1) * iterdiag[, "envvar"]
    ## there is probably a better way to connect mainloci[, "niter"] to inter
    tmp <- c(array(mainloci[index, "niter"], c(2,nindex))[1, ])
    tmp <- match(tmp, iterdiag[, "niter"])
    rss <- unlist(tapply(rss[tmp], inter, mean))
    
    if(any(type == "detection")) {
      npar <- unlist(tapply(npar[tmp], inter, mean))
      npar[is.na(npar)] <- 0
      for(tri in c("lower","upper")) {
      }
    }
    
    ## calculate LPD or other diagnostic.
    for(tri in c("lower","upper")) {
      nscan <- length(scan[tri])
      if(is.lod[tri]) {
        accum[, tri] <- nind.pheno *
          log((rss + env * nscan + nind.pheno * accum[, tri]) / rss)
        if(type[tri] == "LPD")
          accum[, tri] <- accum[, tri] / (2 * log(10))
        else if(type[tri] == "LR")
          accum[, tri] <- accum[, tri] / 2
        else if(type[tri] == "detection") {
          ## This is restricts detection prior to sampled loci pairs.
          p1 <- exp(accum[, tri] / 2)
          detect.prior = 1 / nrow(accum)
          accum[, tri] <- p1 * detect.prior / (1 + (p1 - 1) * detect.prior)
          accum[is.na(accum[, tri]), tri] <- 0.5
        }
      }
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
  if(length(scan.one)) {
    if(verbose)
      cat("qb.scanone on diagonal with", paste(scan.one, collapse = ","),
          "...\n")
    qb.scan$one <- qb.scanone(qbObject, epistasis,
                          scan.one, type["lower"], sum.scan = "only",
                          covar = covar, min.iter = min.iter,
                          verbose = FALSE)
  }
  else {
    if(verbose)
      cat("diagonal set to zero\n")
    qb.scan$one <- rep(0, nrow(pull.grid(qbObject)))
  }

  ## Assign attributes.
  attr(qb.scan, "class") <- c("qb.scantwo", "list")
  attr(qb.scan, "method") <- type
  attr(qb.scan, "scan") <- scan.names
  attr(qb.scan, "min.iter") <- min.iter
  attr(qb.scan, "type") <- qb.get(qbObject, "cross")
  attr(qb.scan, "chr") <- chr
  attr(qb.scan, "weight") <- weight
  attr(qb.scan, "pheno.name") <- pheno.name
  attr(qb.scan, "map") <- map
  attr(qb.scan, "qb") <- qb.name
  qb.scan
}
###################################################################
summary.qb.scantwo <- function(object,
                             chr = NULL, ## Must be integer for now.
                             threshold = 0,
                             sort = "no",
                             digits = 3,
                             which.pos = "upper",
                             min.iter = attr(object, "min.iter"),
                             refine = FALSE, width = 10, smooth = 3,
                             n.qtl = 0.05,
                             verbose = FALSE,
                             ...)
{
  ## new intertwo needs to be checked out
  ## need pos1 and pos2 for lower and upper separately
  ## chr not working?
  
  qbObject <- get(attr(object, "qb"))
  chr.qb <- attr(object, "chr")
  if(!is.null(chr.qb))
    qbObject <- subset(qbObject, chr = chr.qb)
  
  pheno.name <- attr(object, "pheno.name")

  ## Get position pairs.
  gridtwo <- qb.intertwo(qbObject, min.iter)
  inter <- paste(gridtwo[1, ], gridtwo[3, ], sep = ".")

  ## Get unique pairs of chromosomes.
  tmp <- order(gridtwo[1, ], gridtwo[3, ])
  uinter <- unique(inter[tmp])
  chrs <- as.matrix(gridtwo[c(1,3), tmp[!duplicated(inter[tmp])]])
  if(!missing(chr)) {
    tmp <- !is.na(match(uinter, c(outer(chr, chr, paste, sep = "."))))
    uinter <- uinter[tmp]
    chrs <- as.matrix(chrs[, tmp])
    keep <- !is.na(match(gridtwo[1, ], chr)) & !is.na(match(gridtwo[3, ], chr))
  }
  else
    keep <- rep(TRUE, ncol(gridtwo))
  
  out <- matrix(0, length(uinter), 9)
  dimnames(out) <- list(uinter, c("chr1", "chr2",
                                  "n.qtl",
                                  "l.pos1", "l.pos2", "lower",
                                  "u.pos1", "u.pos2", "upper"))

  n.iter <- nrow(qb.get(qbObject, "iterdiag"))

  out[, c("chr1", "chr2")] <- t(chrs)
  out[, "n.qtl"] <- tapply(gridtwo["niter", keep], inter[keep], sum) / n.iter
  rm(keep)

  ## Restrict to pairs with at least n.qtl estimated QTL.
  tmp <- out[, "n.qtl"] >= n.qtl
  if(sum(tmp)) {
    if(sum(tmp) == 1) {
      outn <- dimnames(out)
      out <- matrix(out[tmp,], 1)
      dimnames(out) <- list(outn[[1]][tmp], outn[[2]])
    }
    else
      out <- out[tmp, ]
  }
  else
    out <- NULL
  
#  if(any(center == "mean")) for(i in uinter) {
#    ii <- (i == inter)
#    if(any(ii)) {
#      out[i, 3:6] <- if(sum(ii) == 1)
#        c(gridtwo[c(2,4), ii], object$two[ii, ])
#      else {
#        apply(cbind(t(gridtwo[c(2,4), ii]), object$two[ii, ]), 2,
#              weighted.mean, gridtwo["niter", ii])
#      }
#    }
#  }
#
#  ## Use mode for LPD, heritability, posterior, ...
#  if(any(center == "mode")) {
#    ## Convert object to scantwo after smoothing.
#    x2 <- qb.scantwo.smooth(object, chr, smooth, ...)
#    for(tri in c("lower","upper")) if(center[tri] == "mode") {
#      if(tri == "upper")
#        x2$lod <- t(x2$lod)
#      tmp <- summary(x2)
#      tmp2 <- paste(tmp$chr1,tmp$chr2,sep=".")
#      tmp2 <- match(dimnames(out)[[1]],tmp2)
#      if(tri == "upper") {
#        out[, tri] <- tmp[tmp2, "lod.int"]
#        x2$lod <- t(x2$lod)
#      }
#      else
#        out[, tri] <- tmp[tmp2, "lod.joint"]
#      if (which.pos == tri) 
#        for (i in c("pos1", "pos2"))
#          out[, i] <- tmp[tmp2, i]
#    }
#  }

  if(!is.null(out)) {
    x2 <- qb.scantwo.smooth(object, chr, smooth, gridtwo, ...)

    type <- attr(x2, "method")

    ## Center as mean for variance, estimate, cellmean.
    ## Center as mode for all other types.
    center <- character()
    for(tri in names(type)) {
      center[tri] <- "mean"
      if(is.na(match(type[tri], c("variance","estimate","cellmean"))))
        center[tri] <- "mode"
    }
    
    ## Weighted means by chr.
    ## Could do better row.names using names(cross$geno).
    if(verbose) {
      cat(paste(c("upper:","lower:"), type[c("upper","lower")],
                "of", pheno.name,
                "for", attr(x2, "scan")[c("upper","lower")],
                collapse = "\n"),
          "\n")
      if(min.iter > 1)
        cat("Including only loci pairs with at least", min.iter, "samples.\n")
      if(any(threshold != 0)) {
        cat("Thresholds:",
            paste(names(threshold), c(threshold), collapse = ", ", sep = "="),
            "\n")
      }
      cat("\n")
    }

    x.iter <- qb.smoothtwo(x2$map, x2$nitertwo, x2$niterone, x2$nitertwo,
                            smooth)
    ## somehow x.iter does not agree with posterior
    ## and I am not sure I got the transposes right below
    lod <- x2$lod
    x2$lod <- t(lod)
    x2$lod[row(x.iter) > col(x.iter)] <- x.iter[row(x.iter) > col(x.iter)]
    tmp <- summary(x2) ## lod.int has lower, lod.joint has max main posterior
    tmp2 <- paste(tmp$chr1,tmp$chr2,sep=".")
    tmp2 <- match(dimnames(out)[[1]],tmp2)
    for(i in c("pos1","pos2"))
      out[, paste("l", i, sep = ".")] <- tmp[tmp2, i]
    out[, "lower"] <- tmp[tmp2, "lod.int"]
    
    x2$lod <- lod
    x2$lod[row(x.iter) > col(x.iter)] <- t(x.iter)[row(x.iter) > col(x.iter)]
    tmp <- summary(x2) ## lod.int has upper, lod.joint has max epis posterior
    tmp2 <- paste(tmp$chr1,tmp$chr2,sep=".")
    tmp2 <- match(dimnames(out)[[1]],tmp2)
    for(i in c("pos1","pos2"))
      out[, paste("u", i, sep = ".")] <- tmp[tmp2, i]
    out[, "upper"] <- tmp[tmp2, "lod.int"]
    
    x2$lod <- lod
    
    rm(x.iter, lod, tmp, tmp2)
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
  if(!is.null(out) & refine & any(center == "mode")) {
    if(verbose)
      cat("Refining estimates of mode.\n")

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
                                 type, smooth)
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
  if(is.null(out))
    out
  else
    round(out, digits)
}
###################################################################
print.qb.scantwo <- function(x, ...) print(summary(x, ...))
###################################################################
qb.scantwo.slice <- function(x2, chr, slice, type, smooth)
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
  
  typ <- type[c("lower", "upper")[(1 + slice["upper"])]]
  lod2[is.na(lod2)] <- 0
  grid[[typ]] <- rep(NA, nrow(grid))
  grid[is.chr, typ] <- qb.smoothone(lod2, grid[is.chr, ], smooth,
                                     x2$niterone[is.chr])
  
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
                                              x2$niterone[is.chr])
    }
  }

  ## Reduce down to desired chr.
  grid <- grid[is.chr, ]
  
  ## Make grid a scanone object.
  class(grid) <- c("scanone", "data.frame")
  attr(grid, "method") <- typ
  attr(grid, "type") <- type
  attr(grid, "model") <- "normal"
  grid
}
###################################################################
qb.scantwo.smooth <- function(x,
                            chr = NULL,
                            smooth = 3,
                            gridtwo = qb.intertwo(qbObject, min.iter),
                            ...)
{
  qbObject <- get(attr(x, "qb"))
  chr.qb <- attr(x, "chr")
  if(!is.null(chr.qb))
    qbObject <- subset(qbObject, chr = chr.qb)

  type <- attr(x, "method")
  scan <- attr(x, "scan")
  min.iter <- attr(x, "min.iter")
  weight <- attr(x, "weight")

  ## Get sampling grid.
  gridone <- pull.grid(qbObject, offset = TRUE, spacing = TRUE)
  ## Force getting of 2-D sampling grid as well.
  i.lower <- paste(gridtwo[3,], gridtwo[4, ], sep = ":")
  

  ## Subset index for selected chromosomes.
  ## Note careful handshaking below to match up chr.sub.
  if(!is.null(chr)) {
    if(!is.numeric(chr))
      stop("chr must be numeric index to chromosomes")
    chr.sub <- match(gridone$chr, chr)
    chr.sub <- !is.na(chr.sub)
    if(!sum(chr.sub))
      stop(paste("no samples for chromosomes",
                 chr[sort(unique(chr.sub))],
                 collapse = ","))
    qbObject <- subset(qbObject, chr = chr)
    gridone <- pull.grid(qbObject, offset = TRUE, spacing = TRUE)
  }
  else {
    chr <- sort(unique(gridone$chr))
    chr.sub <- rep(TRUE, nrow(gridone))
  }

  nmap <- sum(chr.sub)

  ## Get indices into lod matrix.
  tmp <- unclass(ordered(paste(gridtwo[1,], gridtwo[2, ], sep = ":"),
                         paste(gridone$chr, gridone$map, sep = ":")))
  i.lower <- unclass(ordered(i.lower,
                             paste(gridone$chr, gridone$map, sep = ":")))
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
    type["upper"] <- type["lower"]
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
  niterone <- unclass(table(qb.inter(qbObject, gridone)))

  ## Smooth lod matrix by chromosome.
  lod <- qb.smoothtwo(gridone, nitertwo, niterone, lod, smooth)

  ## Make a scantwo object.
  lst <- list(lod = lod, map = gridone, scanoneX = NULL,
              niterone = niterone, nitertwo = nitertwo)
  attr(lst, "class") <- "scantwo"
  attr(lst, "method") <- type
  attr(lst, "type") <- attr(x, "type")
  attr(lst, "scan") <- scan
  invisible(lst)
}
###################################################################
plot.qb.scantwo <- function(x,
                          chr = NULL,
                          smooth = 0,
                          main = mains,
                          offset = offsets,
                          lower = "joint",
                          nodiag = all(diag(x2$lod) == 0),
                          slice = NULL,
                          show.locus = TRUE,
                          verbose = FALSE,
                          ...)
{
  chrs <- chr
  if(!is.null(slice))
    chrs <- c(chrs, slice[1])
  x2 <- qb.scantwo.smooth(x, chrs, smooth, ...)

  qbObject <- get(attr(x, "qb"))
  chr.qb <- attr(x, "chr")
  if(!is.null(chr.qb))
    qbObject <- subset(qbObject, chr = chr.qb)
  
  pheno.name <- attr(x, "pheno.name")
  type <- attr(x2, "method")
  scan <- attr(x2, "scan")
  min.iter <- attr(x, "min.iter")

  if(verbose) {
    cat(paste(c("\nupper:","lower:"), type[c("upper","lower")], "of",
              pheno.name, "for", scan[c("upper","lower")],
              collapse = "\n"),
        "\n")
    if(min.iter > 1)
      cat("Including only loci pairs with at least", min.iter, "samples.\n")
  }
  mains <- paste(type[c("upper","lower")], "of",
                 scan[c("upper","lower")], collapse = " / ")
  
  if(is.null(slice)) {
    ## Rescale values if type is estimate.
    tmpfn <- function(x) {
      max(x, -x, na.rm = TRUE)
    }
    offsets <- c(lower = 0, upper = 0)
    if(type["upper"] == "estimate")
      offsets["upper"] <- tmpfn(lod[row(lod) < col(lod)])
    if(type["lower"] == "estimate") {
      offsets["lower"] <- tmpfn(lod[row(lod) >= col(lod)])
    }
    if(is.null(names(offset)))
      names(offset) <- names(offsets)
    if(type["upper"] == "estimate") {
      lod[row(lod) < col(lod)] <- 
        1 + (lod[row(lod) < col(lod)] / offset["upper"])
    }
    if(type["lower"] == "estimate") {
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
    
    ## Plot scantwo object.
    plot(x2, lower = lower, nodiag = nodiag, main = main,
         incl.markers = TRUE, ...)
    if(verbose)
      cat("\n")
    invisible(x2)
  }
  else { ## 1-D slice through 2-D surface
    grid <- qb.scantwo.slice(x2, chr, slice, type, smooth)

    ## Plot slice.
    if(var(grid[[4]]) > 0 & show.locus) {
      tmpar <- par(mfrow=c(2,1), mar=c(2.1,4.1,0.1,0.1))
      on.exit(par(tmpar))
    }
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
                          offdiag = 0.5)
{
  if(smooth) {
    if(offdiag < 0)
      offdiag <- 0
    if(offdiag > 1)
       offdiag <- 1
    
    smoothtwo <- function(x, smooth, nitertwo) {
      n.map <- dim(x)
      if(min(n.map) > 3) {
        nr <- n.map[1]
        nc <- n.map[2]
        o <- (x != 0)
        for(i in seq(smooth)) {
          ## Set up numerator = weighted sum of xs.
          wt <- nitertwo * x
          x <- (wt[, c(1, seq(nc - 1))] +
                wt[, c(seq(2, nc), nc)] +
                wt[c(1, seq(nr - 1)), ] +
                wt[c(seq(2, nr), nr), ])
          if(offdiag)
            x <- x + offdiag *
              (wt[c(1, seq(nr - 1)), c(1, seq(nc - 1))] +
               wt[c(1, seq(nr - 1)), c(seq(2, nc), nc)] +
               wt[c(seq(2, nr), nr), c(seq(2, nc), nc)] +
               wt[c(seq(2, nr), nr), c(1, seq(nc - 1))])
            
          if(any(o))
            x[o] <- (x + 4 * (1 + offdiag) * wt)[o]
          
          ## Now get denominator = sum of weights.
          wt <- (nitertwo[, c(1, seq(nc - 1))] +
                 nitertwo[, c(seq(2, nc), nc)] +
                 nitertwo[c(1, seq(nr - 1)), ] +
                 nitertwo[c(seq(2, nr), nr), ])
          if(offdiag)
            wt <- wt + offdiag *
              (nitertwo[c(1, seq(nr - 1)), c(1, seq(nc - 1))] +
               nitertwo[c(1, seq(nr - 1)), c(seq(2, nc), nc)] +
               nitertwo[c(seq(2, nr), nr), c(seq(2, nc), nc)] +
               nitertwo[c(seq(2, nr), nr), c(1, seq(nc - 1))])
          ## Off-diagonal elements.
          if(any(o))
            wt[o] <- (wt + 4 * (1 + offdiag) * nitertwo)[o]
          x <- x / wt
          x[is.na(x)] <- 0
        }
      }
      x
    }
    smoothtwo.same <- function(x, smooth, nitertwo) {
      ## Smooth upper and lower half of x, leaving diagonal unchanged.
      is.upper <- row(x) > col(x)
      is.lower <- row(x) < col(x)

      ## Make mat symmetric using lower triangle.
      tmpfn <- function(x, smooth, nitertwo, is.upper) {
        tmpfn2 <- function(x) {
          mat <- x
          mat[is.upper] <- t(x)[is.upper]
          diagmat <- mat[row(mat) == 1 + col(mat)]
          ndiag <- length(diagmat)
          diag(mat) <- (diagmat[c(1, seq(ndiag))] +
                        diagmat[c(seq(ndiag), ndiag)]) / 2
          mat
        }
        smoothtwo(tmpfn2(x), smooth, tmpfn2(nitertwo))
      }
      x[is.lower] <- tmpfn(x, smooth, nitertwo, is.upper)[is.lower]
      x[is.upper] <- t(tmpfn(t(x), smooth, t(nitertwo), is.upper))[is.upper]
      x
    }
    chrs <- unique(grid$chr)
    n.chr <- length(chrs)
    if(n.chr == 1)
      x <- smoothtwo.same(x, smooth, nitertwo)
    else {
      for(i in seq(n.chr)) {
        ## process diagonal matrix.
        rows <- chrs[i] == grid$chr
        if(sum(rows)) {
          x[rows,rows] <- smoothtwo.same(x[rows,rows], smooth,
                                         nitertwo[rows,rows])
          ## Now off diagonal matrices.
          if(i < n.chr) for(j in seq(i + 1, n.chr)) {
            cols <- chrs[j] == grid$chr
            if(sum(cols)) {
              ## Lower triangle matrix.
              x[rows,cols] <- smoothtwo(x[rows,cols], smooth,
                                        nitertwo[rows,cols])
              ## Upper triangle matrix.
              x[cols,rows] <- smoothtwo(x[cols,rows], smooth,
                                        nitertwo[cols,rows])
            }
          }
        }
      }
    }
    diag(x) <- qb.smoothone(diag(x), grid, smooth, niterone)
  }
  x
}
