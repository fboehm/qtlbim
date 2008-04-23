#####################################################################
##
## $Id: cond.R,v 1.11.2.10 2006/12/12 19:23:28 byandell Exp $
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
qb.condloci <- function(qbObject, chr = 1, cutoff = 25, nqtl = uqtl[nuqtl],
                        use.qtl = FALSE, individual = TRUE)
{
  niter <- qb.niter(qbObject)
  
  mainloci <- qb.get(qbObject, "mainloci")
  mainloci <- mainloci[mainloci$chrom == chr, c("niter","nqtl","locus")]
  if(!nrow(mainloci))
    return(NULL)

  ## Reassign nqtl as number of QTL on chr.
  mainloci$nqtl <- c(table(mainloci$niter)[as.character(mainloci$niter)])
  
  ## Posterior probabilities for nqtl on chr.
  prob <- table(mainloci$nqtl)
  uqtl <- as.numeric(names(prob))
  prob <- prob / uqtl
  prob <- 100 * prob / niter

  ## Only include nqtl with posterior above cutoff.
  lt.cut <- rev(cumsum(rev(prob))) >=  cutoff
  nuqtl <- max(1, sum(lt.cut))
  ## Override nuqtl if nqtl is not missing.
  nuqtl <- which.min(abs(nqtl - uqtl))
  if(uqtl[nuqtl] != nqtl)
    stop(paste("No samples with", nqtl, "QTL"))
  
  ## If any extra, collapse them into one category.
  if(nuqtl < length(prob)) {
    prob <- c(prob[seq(nuqtl)], sum(prob[-seq(nuqtl)]))
    uqtl <- c(uqtl[seq(nuqtl)], nqtl + 1)
    mainloci$nqtl[mainloci$nqtl > nqtl] <- nqtl + 1
  }

  if(nqtl > 1) {
    ## Match QTL number with QTL of largest list using medians.
    tmpdata <- data.frame(locus = mainloci$locus[mainloci$nqtl == nqtl])
    tmpdata$qtl <- rep(seq(nqtl), nrow(tmpdata) / nqtl)
    ## Linear discriminant analysis.
    fit.lda <- lda(as.matrix(tmpdata$locus),tmpdata$qtl)
    if(individual) {
      ## Pick using LDA for each sample. Does not protect against ties.
      mainloci$qtl <- c(predict(fit.lda, as.matrix(mainloci$locus))$class)
    }
    else {
      ## Pick using LDA on mode for index conditional on nqtl.
      ## Rarely get ties except i > nuqtl.
      reassign.qtl <- function(fit.lda, locus, uqtli, tmpdata) {
        tmpqtl <- rep(seq(uqtli), nrow(tmpdata) / uqtli)
        mode.qtl <- tapply(locus, tmpqtl,
                           function(locus) {
                             f <- density(locus)
                             f$x[which.max(f$y)]
                           })
        c(predict(fit.lda, as.matrix(mode.qtl))$class[tmpqtl])
      }
      mainloci$qtl <- rep(1, nrow(mainloci))
      myseq <- seq(length(uqtl))
      if(use.qtl)
        myseq <- myseq[-nuqtl]
      for(i in myseq) {
        uqtli <- uqtl[i]
        mainloci$qtl[mainloci$nqtl == uqtli] <-
          reassign.qtl(fit.lda, mainloci$locus[mainloci$nqtl == uqtli],
                       uqtli, tmpdata)
      }
    }
    if(use.qtl)
      mainloci$qtl[mainloci$nqtl == nqtl] <- tmpdata$qtl
  }
  else
    mainloci$qtl <- rep(1, nrow(mainloci))

  class(mainloci) <- c("qb.condloci", "data.frame")
  attr(mainloci, "chr") <- chr
  attr(mainloci, "cutoff") <- cutoff
  attr(mainloci, "prob") <- prob
  attr(mainloci, "uqtl") <- uqtl
  attr(mainloci, "nqtl") <- nqtl
  attr(mainloci, "niter") <- niter
  attr(mainloci, "step") <- qb.get(qbObject, "step")
  mainloci
}
##############################################################################
print.qb.condloci <- function(x,...) print(summary(x, ...))
##############################################################################
summary.qb.condloci <- function(object, merge = TRUE, ...)
{
  nqtl <- attr(object, "nqtl")
  prob <- attr(object, "prob")
  niter <- attr(object, "niter")

  my.summary <- function(x) {
    tmp <- factor(paste("QTL", x$qtl))
    out <- t(sapply(tapply(x$locus,tmp,
                           function(x) {
                             c(summary(x),
                               Pct. = round(100 * length(x) / niter, 2))
                           }),
                    c))
    ties <- apply(table(x$niter, tmp) > 1, 2, sum) / niter
    cbind(out, Ties = round(100 * ties[rownames(out)], 2))
  }
  if(merge | length(prob) == 1)
    my.summary(object)
  else {
    tmp <- split(object, object$nqtl)
    tmp2 <- c("=",">=")[1 + (as.numeric(names(tmp)) > nqtl)]
    names(tmp) <- paste("nqtl", tmp2, names(tmp))
    lapply(tmp, my.summary)
  }
}
##############################################################################
plot.qb.condloci <- function(x, merge = TRUE, jitter = 0.5, ...)
{
  prob <- round(attr(x, "prob"))
  nqtl <- attr(x, "nqtl")
  niter <- attr(x, "niter")
  
  step <- attr(x, "step")
  chr <- attr(x, "chr")
  x$locus <- jitter(x$locus, amount = step * jitter)
  
  if(merge | length(prob) == 1) {
    ## Density plot averaged over number of QTL.
    tmp <- round(100 * table(x$qtl) / niter)
    densityplot(~locus, x, groups = x$qtl,
                main = paste("chr = ", chr, ", nqtl = ", nqtl,
                  " (", paste(tmp, collapse = ", "), "%)", sep = ""),
                ...)
  }
  else {
    ## Density plots conditioning on number of QTL.
    tmp <- c("=",">=")[1 + (x$nqtl > nqtl)]
    x$nqtl <- factor(paste("chr = ", chr, ", nqtl ", tmp, " ", x$nqtl,
                           " (", prob[x$nqtl], "%)", sep = ""))
    densityplot(~locus | nqtl, x, groups = x$qtl,
                layout = c(1,length(prob)))
  }
}
##############################################################################
qb.epimodes <- function(qbObject, cutoff = 1, nqtl = nqtl.est,
                        n.iter = qb.niter(qbObject),
                        pairloci = qb.get(qbObject, "pairloci"), ...)
{
  pairloci <- pairloci[pairloci$chrom1 == pairloci$chrom2, ]
  if(!nrow(pairloci))
    return(NULL)

  ## Subset on chromosomes.
  region <- qb.get(qbObject, "region")
  if(!is.null(region)) {
    tmp <- region$chr[region$start < region$end]
    geno.names <- geno.names[tmp]
  }

  ## Posterior probabilities on number of QTL per chromosome.
  nqtl.post <- tapply(pairloci$niter, pairloci$chrom1,
                      function(x) {
                        if(length(x)) {
                          ## prob = probability that number of QTL per chr is uqtl.
                          prob <- table(table(x)[as.character(x)])
                          uqtl <- as.numeric(names(prob))
                          prob <- prob / uqtl
                          prob <- 100 * prob / n.iter
                          tmp <- 100 - sum(prob)
                          if(tmp > 0.01)
                            prob <- c("0" = tmp, prob)
                          prob
                        }
                        else
                          0
                      },
                      simplify = FALSE)
  nqtl.post <- nqtl.post[geno.names]

  ## Estimate number of QTL per chromosome.
  nqtl.est <- sapply(nqtl.post, function(x) {
    if(sum(x)) {
      uqtl <- as.numeric(names(x))
      uqtl[max(1, sum(rev(cumsum(rev(x))) >=  cutoff))]
    }
    else
      0
  })
  nqtl.est <- nqtl.est[geno.names]
  
  if(!missing(nqtl)) {
    if(length(nqtl) != length(nqtl.est))
      stop("nqtl length must match chromosomes in qb object")
    nqtl <- as.array(nqtl)
    names(nqtl) <- names(nqtl.est)
  }
  get.mode <- function(x, min.mode = FALSE, f = density(x)) {
    ## Find location of smoothed maximum or minimum.
    r <- range(x)
    tmp <- f$x >= r[1] & f$x <= r[2]
    f$x[tmp][ifelse(min.mode, which.min(f$y[tmp]), which.max(f$y[tmp]))]
  }
  chrsplit <- function(pairloci, chr, nqtl) {
    all.loci <- stack(pairloci[pairloci$chrom1 == chr, c("locus1","locus2")])
    if(nqtl == 1) {
      list(peaks = get.mode(all.loci$values))
    }
    else {
      ## Linear discriminant analysis between elements of epistatic pairs.
      pred.lda <- c(predict(lda(as.matrix(all.loci$values), all.loci$ind),
                            as.matrix(all.loci$values))$class)
      
      ## Find peaks, then locate valleys between peaks.
      ## Use density across chromosome.
      peaks <- tapply(all.loci$values, pred.lda, get.mode, FALSE)
      if(length(peaks) > 1) {
        valleys <- tapply(all.loci$values,
                          cut(all.loci$values, peaks, labels = FALSE),
                          get.mode, TRUE)
        ## Re-estimate peaks given valleys.
        r <- c(min(all.loci$values) - 1,
               valleys,
               max(all.loci$values) + 1)
        peaks <- tapply(all.loci$values,
                        cut(all.loci$values, r, labels = FALSE),
                        get.mode, FALSE)
      }
      else
        valleys <- NULL
      list(peaks = peaks, valleys = valleys)
    }
  }
  
  peaks <- valleys <- list()
  for(i in names(nqtl)) {
    if(nqtl[i]) {
      tmp <- chrsplit(pairloci, i, nqtl[i])
      peaks[[i]] <- tmp$peaks
      if(nqtl[i] > 1) 
        valleys[[i]] <- tmp$valleys
    }
  }
  out <- list(nqtl.post = nqtl.post, nqtl.est = nqtl,
              peaks = peaks, valleys = valleys)
  class(out) <- c("qb.mainmodes", "list")
  out
}
##############################################################################
qb.mainmodes <- function(qbObject, cutoff = 25, nqtl = NULL,
                        n.iter = qb.niter(qbObject),
                        mainloci = qb.get(qbObject, "mainloci"), ...)
{
  geno.names <- qb.geno.names(qbObject)
  chrom.names <- ordered(geno.names[mainloci$chrom], geno.names)
  
  ## Subset on chromosomes.
  region <- qb.get(qbObject, "region")
  if(!is.null(region)) {
    tmp <- region$chr[region$start < region$end]
    geno.names <- geno.names[tmp]
  }

  ## Posterior probabilities on number of QTL per chromosome.
  nqtl.post <- tapply(mainloci$niter, chrom.names,
                      function(x) {
                        if(length(x)) {
                          ## prob = probability that number of QTL per chr is uqtl.
                          prob <- table(table(x)[as.character(x)])
                          uqtl <- as.numeric(names(prob))
                          prob <- prob / uqtl
                          prob <- 100 * prob / n.iter
                          tmp <- 100 - sum(prob)
                          if(tmp > 0.01)
                            prob <- c("0" = tmp, prob)
                          prob
                        }
                        else
                          0
                      }, simplify = FALSE)
  nqtl.post <- nqtl.post[geno.names]
  
  ## Estimate number of QTL per chromosome.
  nqtl.est <- sapply(nqtl.post, function(x) {
    if(sum(x)) {
      uqtl <- as.numeric(names(x))
      uqtl[max(1, sum(rev(cumsum(rev(x))) >=  cutoff))]
    }
    else
      0
  })
  nqtl.est <- nqtl.est[geno.names]
  
  if(is.null(nqtl))
    nqtl <- nqtl.est

  if(!missing(nqtl)) {
    n.chr <- length(nqtl.est)
    if(length(nqtl) != n.chr) {
      if(length(nqtl) == 1)
        nqtl <- rep(nqtl, n.chr)
      else
        stop(paste("nqtl length must be number of chromosomes (",
                   n.chr, ") in qb object", sep = ""))
    }
    nqtl <- as.array(nqtl)
    names(nqtl) <- names(nqtl.est)
  }
  get.mode <- function(x, min.mode = FALSE, f = density(x)) {
    ## Find location of smoothed maximum or minimum.
    r <- range(x)
    tmp <- f$x >= r[1] & f$x <= r[2]
    f$x[tmp][ifelse(min.mode, which.min(f$y[tmp]), which.max(f$y[tmp]))]
  }
  chrsplit <- function(mainloci, nqtl) {
    all.loci <- mainloci[, "locus"]
    if(nqtl == 1) {
      list(peaks = get.mode(all.loci))
    }
    else {
      ## all.nqtl = number of qtl for each sample.
      all.nqtl <- mainloci[, "niter"]
      all.nqtl <- table(all.nqtl)[as.character(all.nqtl)]

      ## nqtl.loci = loci with number of QTL = nqtl.
      nqtl.loci <- all.loci[all.nqtl == nqtl]
      ## qtl.id = ID for QTL: 1, 2, ..., nqtl.
      qtl.id <- rep(seq(nqtl), length(nqtl.loci) / nqtl)

      ## Linear discriminant analysis.
      if(max(tapply(nqtl.loci,qtl.id,function(x)diff(range(x))))) {
        pred.lda <- c(predict(lda(as.matrix(nqtl.loci), qtl.id),
                              as.matrix(all.loci))$class)
        
        ## Find peaks, then locate valleys between peaks.
        ## Use density across chromosome.
        peaks <- tapply(all.loci, pred.lda, get.mode, FALSE)
      }
      else ## nqtl.loci only take on a nqtl values. Assume one QTL.
        peaks <- get.mode(all.loci)
    
      if(length(peaks) > 1) {
        valleys <- tapply(all.loci, cut(all.loci, peaks, labels = FALSE),
                          get.mode, TRUE)
        ## Re-estimate peaks given valleys.
        r <- c(min(all.loci) - 1, valleys, max(all.loci) + 1)
        peaks <- tapply(all.loci, cut(all.loci, r, labels = FALSE),
                        get.mode, FALSE)
      }
      else
        valleys <- NULL
      
      list(peaks = peaks, valleys = valleys)
    }
  }

  peaks <- valleys <- list()
  for(chri in names(nqtl)) {
    if(nqtl[chri]) {
      tmp <- chrsplit(mainloci[chrom.names == chri, , drop = FALSE], nqtl[chri])
      peaks[[chri]] <- tmp$peaks
      if(nqtl[chri] > 1) 
        valleys[[chri]] <- tmp$valleys
    }
  }
  out <- list(nqtl.post = nqtl.post, nqtl.est = nqtl,
              peaks = peaks, valleys = valleys)
  class(out) <- c("qb.mainmodes", "list")
  out
}
##############################################################################
print.qb.mainmodes <- function(x, ...) print(summary(x, ...))
##############################################################################
summary.qb.mainmodes <- function(object, digits = 4, ...)
{
  n.chr <- length(object$nqtl.est)
  max.nqtl <- max(object$nqtl.est)

  probs <- sort(as.numeric(unique(unlist(lapply(object$nqtl.post, names)))))
  if(n.chr == 1)
    object$nqtl.post <- round(object$nqtl.post[[1]], 1)
  else {
    tmp <- matrix("", length(probs), n.chr)
    dimnames(tmp) <- list(as.character(probs), names(object$nqtl.post))
    for(i in names(object$nqtl.post)) {
      tmpi <- object$nqtl.post[[i]]
      tmp[names(tmpi), i] <- round(tmpi, 1)
    }
    object$nqtl.post <- tmp
  }

  ## Round off modes and arrange in truncated array.
  for(modes in c("peaks", "valleys")) {
    if(length(object[[modes]])) {
      if(n.chr == 1)
        object[[modes]] <- as.character(signif(object[[modes]][[1]], digits))
      else {
        object[[modes]] <- lapply(object[[modes]], signif, digits)
        n.rows <- max.nqtl - (modes == "valleys")
        tmp <- names(object[[modes]])
        object[[modes]] <- sapply(object[[modes]],
                                  function(x) {
                                    c(x, rep("", n.rows - length(x)))
                                  })
        if(n.rows == 1)
          names(object[[modes]]) <- tmp
      }
    }
  }
  class(object) <- c("summary.qb.mainmodes", "list")
  object
}
##############################################################################
print.summary.qb.mainmodes <- function(x, ...)
{
  ## Need to add geno.names.
  n.chr <- length(x$nqtl.est)

  cat("Posterior distribution on number of QTL")
  if(n.chr > 1) cat("\n")
  print(x$nqtl.post, quote = F)
  cat("\n")

  cat("Estimated number of QTL")
  if(n.chr == 1) {
    tmp <- c(x$nqtl.est)
    names(tmp) <- NULL
    cat(":", tmp, "\n")
  }
  else {
    cat(" per chromosome:\n")
    print(x$nqtl.est)
  }

  cat("\nPeaks:\n")
  print(x$peaks, quote = FALSE)

  if(length(x$valleys)) {
    cat("\nValleys:\n")
    print(x$valleys, quote = FALSE)
  }
}
##############################################################################
qb.chrsplit <- function(grid, mainloci, chrs = seq(geno.names), n.iter,
                        geno.names = levels(ordered(grid$chr)),
                        split.chr = qb.mainmodes(n.iter = n.iter,
                          mainloci = mainloci, ...)$valleys,
                        ...)
{
  ## Only use chromosomes in chrs vector.
  dont.use <- is.na(match(grid$chr, chrs))
  
  gchr <- geno.names[grid$chr]
  
  split.names <- names(split.chr)
  split.names <- split.names[match(geno.names[chrs], split.names, nomatch = 0)]
  for(i in split.names) {
    tmp <- geno.names[grid$chr] == i
    if(any(tmp)) { ## Just to make sure there are some.
      gname <- gchr[tmp][1]
      gchr[tmp] <- NA
      split.chri <- range(grid$pos[tmp])
      split.chri <- c(split.chri[1] - 1, split.chr[[i]], split.chri[2])
      for(j in seq(length(split.chri) - 1)) {
        tmp2 <- (grid$pos[tmp] > split.chri[j] &
                 grid$pos[tmp] <= split.chri[j+1])
        gchr[tmp][tmp2] <- paste(gname, j, sep = ".")
      }
    }
  }

  ## This is wasteful, as we assign, then make NA.
  gchr[dont.use] <- NA
  tmp <- !is.na(gchr)
  tmp2 <- !duplicated(gchr[tmp])
  unsplit <- grid$chr[tmp][tmp2]

  grid$chr <- ordered(gchr, unique(gchr[tmp]))
  attr(grid, "unsplit") <- unsplit
  grid
}


