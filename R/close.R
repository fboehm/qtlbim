#####################################################################
##
## $Id: close.R,v 1.0 2007/7/20 byandell Exp $
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
#####################################################################


##############################################################
qb.nulldist <- function(target = NULL, signed = FALSE,
                        score.type = c("sq.atten","attenuation",
                          "variance","recombination","distance"))
{
  ## Compare target to null
  score.type <- match.arg(score.type)

  if(score.type == "variance" | is.null(target)) {
    data.frame(score = 0)
  }
  else {
    ## For each target match sample chrom and compare locus.
    n.target <- nrow(target)

    if(score.type == "distance")
      score <- Inf
    else
      score <- ifelse(score.type == "recombination", 0.5 * n.target, 0)
    if(signed)
      score <- -score
    data.frame(score = score, n.union = n.target)
  }
}
##############################################################
qb.archdist <- function(sample, target = NULL, signed = FALSE,
                        score.type = c("sq.atten","attenuation",
                          "variance","recombination","distance"),
                        qtl = seq(nrow(target)))
{
  ## Proximity of sampled QTL to target QTL.
  ## attenuation = (1-2r), r = recombination fraction.
  ## Uses Haldane map function on cM distance.

  score.type <- match.arg(score.type)
  niter <- sample[, "niter"]

  ## Want to add LPD and 2logBF to this list. Take some care.
  ## Also want variance above null model of no QTL.
  ## variance = variance of model - null model variance.
  ## variance = total variance - envvar.
  ## LPD = log10(envvar/totvar)
  ## BF comes from qb.bf routine.
  if(score.type == "variance" | is.null(target)) {
    tmp <- unlist(tapply(sample[, "variance"], niter, sum))
    score <- tmp - sum(target[, "variance"])
    if(!signed)
      score <- abs(score)
    data.frame(score = score)
  }
  else {
    ## For each target match sample chrom and compare locus.
    n.target <- nrow(target)

    ## Accumulate information about extra sample loci (see below).
    iter.names <- unique(niter)
    n.iter <- length(iter.names)

    tmp <- 0
    if(score.type == "recombination")
      tmp <- 0.5
    else if(score.type == "distance")
      tmp <- Inf
    score <-  rep(tmp, n.iter)
    extra <- rep(0, n.iter)
    names(extra) <- names(score) <- iter.names
    
    distfn <- function(dist, signed) {
      ## Compute score for target.set using sample.set.
      diff.sign <- 1
      if(signed) {
        diff.sign <- sign(dist)
        if(any(diff.sign == 0))
          diff.sign[diff.sign == 0] <- 1
      }
      dist <- abs(dist)
      if(score.type != "distance") {
        dist <- exp(-0.02 * dist)
        if(score.type == "sq.atten")
          dist <- dist * dist
        else if(score.type == "recombination")
          dist <- (1 - dist) / 2
      }
      sum(dist, na.rm = TRUE)
    }

    target$chrom <- ordered(target$chrom, unique(as.character(target$chrom)))
    ## Need to consider matrix of rows[ii], cols[chr == chr[i]].
    ## If more than one, then need to sort out their interplay.
    for(chri in levels(target[, "chrom"])) {
      ii <- target[, "chrom"] == chri
      n.ii <- sum(ii)
      jj <- sample[, "chrom"] == chri
      if(sum(jj)) {
        jj <- seq(nrow(sample))[jj]
        iteri <- sample[jj, "niter"]
        
        ## Distance set = list of c(target.set, sample.set)
        dist.set <- tapply(sample[jj, "locus"], iteri,
                           qb.distset, target[ii, "locus"])
        
        ## Setup for extra and dist computations.
        seq.target <- seq(n.ii)
        chr.target <- seq(n.target)[ii]
        ## sample.set indexes sample matrix.
        sample.set <- unlist(lapply(dist.set,
                                    function(x, seq.target)
                                    x[-seq.target],
                                    seq.target))
        iterj <- iteri[sample.set]
        sample.set <- jj[sample.set]
        ## target.set indexes target matrix.
        target.set <- unlist(lapply(dist.set,
                                    function(x, seq.target, chr.target)
                                    chr.target[x[seq.target]],
                                    seq.target, chr.target))

        ## Only consider selected qtl. 
        is.qtl <- !is.na(match(target.set, qtl))

        if(any(is.qtl)) {
          ## Accumulate number of extra linked loci.
          iter.names <- names(dist.set)
          extra[iter.names] <- extra[iter.names] +
            pmax(0,
                 c(table(iteri)) - sum(!is.na(match(chr.target, qtl))))
          
          ## Score as difference in locus positions.
          iter.names <- as.character(unique(iterj[is.qtl]))
          score[iter.names] <- score[iter.names] +
            unlist(tapply(sample[sample.set[is.qtl], "locus"] -
                          target[target.set[is.qtl], "locus"],
                          iterj[is.qtl], distfn, signed))
        }
      }
    }
    ## Add in extra for unlinked QTL if recombination.
    if(score.type == "recombination")
      score <- score + extra * 0.5
    data.frame(score = score, n.union = n.target + extra)
  }
}
##############################################################
qb.distset <- function(sample, target)
{
  ## 
  ## Distance between architectures for specified chromosome.
  ## target.set = index to target matrix.
  ## sample.set = index to sample matrix.
  
  n.target <- length(target)
  n.sample <- length(sample)
  target.set <- rep(TRUE, n.target)
  sample.set <- rep(TRUE, n.sample)

  if(n.target == n.sample)
    return(c(target.set, sample.set))
  
  if(n.target == 1 | n.sample == 1) {
    ## Find minimal distance of sampled QTL to QTL on chri.
    tmp <- abs(sample - target)
    wh <- which.min(tmp)
    if(n.target > 1)
      target.set[-wh] <- FALSE
    else ## n.sample > 1
      sample.set[-wh] <- FALSE
  }
  else {
    ## Following assumes target and sample are already
    ## ordered by chrom, then locus.
    if(n.sample > n.target)
      sample.set <- qb.distlocus(sample, target)
    else ## n.sample < n.target
      target.set <- qb.distlocus(target, sample)
  }
  c(target.set, sample.set)
}
##############################################################
qb.distlocus <- function(sample, target)
{
  n.target <- length(target)
  n.sample <- length(sample)
  
  ## Generate matrix with all combinations of n.target
  ## subsets fo n.sample.
  sample.locus <- combn(sample, n.target)
  
  ## Find which k-tuple of sample.locus is closest to
  ## which k-tuple of target.locus.
  alldist <- as.matrix(apply(sample.locus, 2,
                             function(x, y) sum((x - y)^2),
                             target))

  ## Find which subset of sample is closest to target.
  wh <- which.min(alldist)
  tmp <- seq(n.sample)
  wh <- combn(tmp, n.target)[, col(alldist)[wh]]

  ## Return sample.set as logical vector.
  !is.na(match(tmp, wh))
}
#####################################################################33
qb.splititer <- function(mainloci, pairloci = NULL, splits)
{
  ## Get loci samples and split by iteration number.
  out <- mainloci[, c("niter","chrom","locus","varadd")]
  names(out)[4] <- "variance"
  var.names <- grep("^var", names(mainloci))
  out$variance <- apply(as.matrix(mainloci[, var.names]), 1, sum)
  chr.pos <- paste(out$chrom, out$locus, sep = ":")
  if(!is.null(pairloci)) if(nrow(pairloci)) {
    index <- paste(out$niter, chr.pos, sep = ":")
    var.names <- grep("^var", names(pairloci))
    pair.var <- apply(as.matrix(pairloci[, var.names]), 1, sum) / 2

    for(i in 1:2) {
      index1 <- paste(pairloci[, "niter"],
                      pairloci[, paste("chrom", i, sep = "")],
                      pairloci[, paste("locus", i, sep = "")],
                      sep = ":")
      tmp <- match(index1, index)
      out[tmp, "variance"] <- out[tmp, "variance"] + pair.var
    }
  }
  out$chrom <- ordered(splits[chr.pos], unique(splits))
  out
}
#####################################################################
qb.close <- function(qbObject, target = NULL, epistasis = TRUE,
                     signed = FALSE,
                     score.type = c("sq.atten","attenuation","variance",
                       "recombination", "distance"), ...)
{
  qb.exists(qbObject)
  
  score.type <- match.arg(score.type)

  if(is.null(target))
    score.type <- "variance"
  else {
    ## If chrom and locus missing from target, assume chr and pos.
    ## This happens if target = summary(qb.hpdone(qbObject)).
    if(inherits(target, "qb.hpdone") | inherits(target, "qb.scanone"))
      target <- summary(target)
    
    ## Make sure target is data.frame with columns: chrom, locus, variance.
    target <- as.data.frame(target)
    if(is.null(target$chrom))
      target$chrom <- target$chr
    if(is.null(target$locus))
      target$locus <- target$pos
    tmp <- match(c("chrom","locus","variance"), names(target), nomatch = 0)
    target <- target[, tmp, drop = FALSE]
    
    if(is.null(target$locus) | is.null(target$chrom))
      stop("target must have both chrom and locus columns")
    
    ## Give target 0 variance if missing.
    if(is.null(target$variance))
      target$variance <- rep(0, nrow(target))

    ## Silly business to ensure target in chrom.locus order.
    target <- target[order(target$chrom, target$locus), ]
  }
  
  ## Get loci samples and split by iteration number.
  mainloci <- qb.get(qbObject, "mainloci", ...)
  if(epistasis)
    pairloci <- qb.get(qbObject, "pairloci", ...)
  else
    pairloci <- NULL
  epistasis <- !is.null(pairloci)

  ## Score each sample against target.
  out <- qb.nulldist(target, signed, score.type)
  iterdiag <- qb.get(qbObject, "iterdiag", ...)
  n.iter <- qb.niter(qbObject)
  out <- out[rep(1, n.iter), ]

  if(is.null(dim(out)))
    out <- data.frame(score = out)
  row.names(out) <- iterdiag$niter

  splits <- find.splits(qbObject, mainloci)
  tmp <- qb.archdist(as.data.frame(qb.splititer(mainloci, pairloci, splits)),
                     target, signed, score.type)
  out[row.names(tmp), ] <- tmp

  ## Pattern of sampled QTL per iteration.
  out$pattern <- qb.makepattern(qbObject, epistasis, mainloci = mainloci, ...)
  
  ## Number of sampled QTL per iteration.
  out$n.qtl <- rep(0, n.iter)
  tmp <- c(table(mainloci$niter))
  out[names(tmp), "n.qtl"] <- tmp

  class(out) <- c("qb.close", "data.frame")
  attr(out, "target") <- target
  attr(out, "epistasis") <- epistasis
  attr(out, "signed") <- signed
  attr(out, "score.type") <- score.type
  out
} 
#######################################################################
print.qb.close <- function(x, ...) print(summary(x, ...))
#######################################################################
summary.qb.close <- function(object,
                             cutoff = ifelse(attr(object, "epistasis"),
                               0.25, 0.5),
                             digits = 0,
                             show = "score",
                             ...)
{
  if(is.na(match(show, names(object))))
    stop("show name must be in names of qb.close object")
  if(!is.numeric(object[[show]]))
    stop("show name must be numeric")

  out <- list()
  out[[1]] <- attr(object, "target")
  
  ## Summary by Number of QTL.
  tmpfn <- function(score, nqtl) {
    tmp <- tapply(score, nqtl, summary)
    tmp.row <- names(tmp)
    tmp.col <- names(tmp[[1]])
    tmp <- matrix(unlist(tmp), ncol = length(tmp.col), byrow = TRUE)
    dimnames(tmp) <- list(tmp.row, tmp.col)
    tmp
  }
  out[[2]] <- tmpfn(object[[show]], object$n.qtl)

  ## Find patterns above cutoff.  
  tmp2 <- paste(object$n.qtl, object$pattern, sep = "@")
  pct <- table(tmp2) * 100 / nrow(object)
  maxpct <- max(pct)

  ## Make sure cutoff allows for at least one entry.
  if(maxpct <= cutoff){
    warning(paste("best pattern cutoff =", cutoff, "too large: max percent =", signif(maxpct, 3)))
    cutoff <- maxpct - 1e-6
  }
  pct <- rev(sort(pct[pct > cutoff]))

  ## Summarize score on most common subset.
  tmp <- !is.na(match(tmp2, names(pct)))
  tmp2 <- tmp2[tmp]
  out[[3]] <- cbind(Percent = round(pct, 1),
                    tmpfn(object$score[tmp], tmp2)[names(pct), ])
  
  names(out) <- c(paste("target for score", attr(object, "score.type")),
                  paste(show, c("by sample number of qtl",
                              "by sample chromosome pattern")))
  class(out) <- c("summary.qb.close", "list")
  out
}
#######################################################################
print.summary.qb.close <- function(x, ...)
{
  for(i in names(x)) {
    cat("\n", i, "\n")
    print(x[[i]])
  }
  invisible()
}
#######################################################################
plot.qb.close <- function(x, category = c("pattern", "nqtl"),
                          xlab = attr(x, "score.type"),
                          cutoff = ifelse(attr(x, "epistasis"),
                            0.25, 0.5),
                          sort.pattern = c("percent", "score"),
                          ...)
{
  sort.pattern <- match.arg(sort.pattern)
  
  if(!missing(cutoff) & missing(category))
    category <- "pattern"
  else
    category <- match.arg(category)
  switch(category,
         nqtl = print(bwplot(ordered(n.qtl) ~ score, x, xlab = xlab, ...), ...),
         pattern = {
           ## For pattern, want to ideally aggregate I think.
           ## At very least, only consider more frequent QTL.
           ## Find patterns above cutoff.  
           tmp2 <- paste(x$n.qtl, x$pattern, sep = "@")
           pct <- table(tmp2) * 100 / nrow(x)
           maxpct <- max(pct)
           pct <- rev(sort(pct[pct > cutoff]))
           if(!length(pct))
             stop(paste("cutoff =", cutoff, "too large: max percent =",
                        signif(maxpct, 3)))
           tmp <- !is.na(match(tmp2, names(pct)))
           data <- data.frame(score = x$score[tmp], pattern = tmp2[tmp])
           tmp <- match(data$pattern, names(pct))
           tmp2 <- paste(names(pct), " (", round(pct, 1), "%)", sep = "")
           data$pattern <- if(sort.pattern == "percent")
             ordered(tmp2[tmp], tmp2)
           else {
             o <- order(- unlist(tapply(data$score,
                                        data$pattern, mean))[names(pct)])
             ordered(tmp2[tmp], tmp2[o])
           }
           print(bwplot(pattern ~ score, data, xlab = xlab, ...), ...)
         })
}
#######################################################################
qb.patterniter <- function(qbObject,
                           epistasis = TRUE,
                           category = "pattern",
                           cutoff = ifelse(epistasis, 0.25, 0.5),
                           mainloci = qb.get(qbObject, "mainloci", ...), ...)
{
  if(category == "pattern") {
    ## Pattern of sampled QTL per iteration.
    pattern <- qb.makepattern(qbObject, epistasis, mainloci = mainloci, ...)
  }
  else
    pattern <- qb.nqtl(qbObject, mainloci = mainloci, ...)

  ## Restrict to most common patterns. Reset cutoff if too large.
  pct <- table(pattern) * 100 / length(pattern)
  if(max(pct) < cutoff)
    cutoff <- max(pct)
  pct <- rev(sort(pct[pct >= cutoff]))

  n.pat <- length(pct)
  if(n.pat == 1)
    stop(paste("only one pattern with cutoff", cutoff))

  ## Number of QTL per iteration.
  nqtl <- rep(0, length(pattern))
  names(nqtl) <- names(pattern)
  tmp <- table(mainloci$niter)
  nqtl[names(tmp)] <- tmp

  list(pattern = pattern, nqtl = nqtl, pct = pct)
}
#######################################################################
qb.patternave <- function(qbObject, epistasis = TRUE, pattern, nqtl, pct,
                          mainloci = qb.get(qbObject, "mainloci", ...),
                          include = c("nested","all","exact"),
                          center = c("median","mean"),
                          level = 5, ...)
{
  center <- match.arg(center)
  if(level <= 0 | level >= 100)
    stop("level must be between 0 and 100")
  level <- c(level / 200, 1 - level / 200)

  restrict <- !is.na(match(pattern, names(pct)))

  if(epistasis)
    pairloci <- qb.get(qbObject, "pairloci", ...)
  else
    pairloci <- NULL

  include <- match.arg(include)
  if(include == "exact") {
    tmp <- rep(restrict, nqtl)
    mainloci <- mainloci[tmp, ]
    if(epistasis) {
      if(!is.null(pairloci))
        pairloci <- pairloci[!is.na(match(pairloci$niter,
                                          unique(mainloci$niter))), ]
    }
    ## Reduce to restricted list of patterns.
    pattern <- pattern[restrict]
    nqtl <- nqtl[restrict]
  }

  ## Accumulate variance by locus (halve epistatic).
  splits <- find.splits(qbObject, mainloci)
  iters <- qb.splititer(mainloci, pairloci, splits)

  ## Set up framework for summary.
  ## Only niter and chrom are correct.
  ## Need to model average to get locus, variance (below).
  if(include == "exact")
    sumpat <- iters[rep(!duplicated(pattern), nqtl), ]
  else
    sumpat <- iters[rep(!duplicated(pattern) & restrict, nqtl), ]

  pattern.nqtl <- rep(pattern, nqtl)

  conf <- matrix(0, nrow(sumpat), 5)
  dimnames(conf) <- list(NULL,
                         c("locus.LCL","locus.UCL","variance.LCL","variance.UCL","n.qtl"))

  if(center == "mean")
    myave <- mean
  else
    myave <- median
  
  ## Score each sample against target.
  if(include == "exact") {
    index <- paste(unlist(apply(as.matrix(nqtl),1,
                                function(x) seq(length = x))),
                   pattern.nqtl, sep = "@")
    index <- ordered(index, unique(index))

    sumpat[, "locus"] <- c(tapply(iters[, "locus"], index, myave,
                                  na.rm = TRUE))
    sumpat[, "variance"] <- c(tapply(iters[, "variance"], index, myave,
                                     na.rm = TRUE))
    targets <- unique(pattern)

    for(i in c("locus","variance")) {
      tmp <- tapply(iters[, i], index, quantile, level)
      conf[, paste(i, "LCL", sep = ".")] <-
        unlist(sapply(tmp, function(x) x[1]))
      conf[, paste(i, "UCL", sep = ".")] <-
        unlist(sapply(tmp, function(x) x[2]))
   }
    conf[, "n.qtl"] <- c(table(index))
  }
  else {

    ## Find patterns that match each target. Watch out for NULL.
    patterns <- unique(pattern)
    targets <- unique(pattern[restrict])
    tmp <- rep(0, length(targets))
    tmp2 <- targets == "NULL"
    if(any(!tmp2))
      tmp[!tmp2] <- table(sumpat[, "niter"])
    pattern.sumpat <- rep(targets, tmp)
   
    if(include == "nested") {
      matches <- qb.match.pattern(qbObject, targets = targets,
                                  exact = FALSE, patterns = patterns)
      tmpfn <- function(chrs, splits.pat, i, pattern.nqtl, patterns, matches)
        ((pattern.nqtl %in% patterns[matches[, i]]) & (splits.pat %in% chrs))
        
    }
    else {
      matches <- NULL
      tmpfn <- function(chrs, splits.pat, i, pattern.nqtl, patterns, matches)
        (splits.pat %in% chrs)
    }

    ## Following will NOT get multiples in targets properly. Too bad.
    ## You can finess this by first splitting chromosomes.

    for(i in targets) {
      chrs <- strsplit(i, ",", fixed = TRUE)[[1]]
      chrs <- chrs[chrs %in% splits]
      
      tmp <- tmpfn(chrs, iters$chrom, i, pattern.nqtl, patterns, matches)
      itersi <- iters[tmp, ]
      splitsi <- factor(iters$chrom[tmp], chrs)

      for(j in c("locus", "variance")) {
        tmp2 <- c(tapply(itersi[, j], splitsi, myave, na.rm = TRUE))
        tmp <- match(chrs, names(tmp2), nomatch = 0)
        sumpat[pattern.sumpat == i, j][tmp > 0] <- tmp2[tmp]

        tmp2 <- tapply(itersi[, j], splitsi, quantile, level)
        conf[pattern.sumpat == i, paste(j, "LCL", sep = ".")][tmp > 0] <-
          unlist(sapply(tmp2, function(x) x[1])[tmp])
        conf[pattern.sumpat == i, paste(j, "UCL", sep = ".")][tmp > 0] <-
          unlist(sapply(tmp2, function(x) x[2])[tmp])

        if(j == "locus") {
          tmp2 <- c(table(splitsi))
          conf[pattern.sumpat == i, "n.qtl"][tmp > 0] <- tmp2[tmp]
        }
      }
    }
  }
  conf[, "n.qtl"] <- conf[, "n.qtl"] / qb.niter(qbObject)
  
  sumpat <- as.data.frame(sumpat)
  tmp <- unlist(strsplit(targets, ",", fixed = TRUE))
  tmp <- tmp[tmp %in% splits]
  sumpat$chrom <- ordered(tmp, unique(splits))
  
  attr(sumpat, "nqtl") <- {
    if(include == "exact")
      nqtl[!duplicated(pattern)]
    else
      nqtl[!duplicated(pattern) & restrict]
  }
  attr(sumpat, "pattern") <- targets
  attr(sumpat, "conf") <- conf
  sumpat
}
#######################################################################
find.splits <- function(qbObject, mainloci = qb.get(qbObject, "mainloci", ...), ...)
{
  grid <- pull.grid(qbObject, offset = TRUE)
  new.chr <- qb.chrsplit(grid, mainloci,
                         split.chr = qb.get(qbObject,"split.chr"))$chr
  names(new.chr) <- paste(grid$chr, grid$pos, sep = ":")
  new.chr
}
#######################################################################
qb.best <- function(...) qb.BestPattern(...)
#######################################################################
qb.BestPattern <- function(qbObject,
                           epistasis = TRUE,
                           category = c("pattern", "nqtl"),
                           cutoff = ifelse(epistasis, 0.25, 0.5),
                           score.type = c("sq.atten","attenuation",
                             "variance","recombination","distance"),
                           include = c("nested","all","exact"),
                           center = c("median","mean"),
                           level = 5, ...)
{
  qb.exists(qbObject)
  
  include <- match.arg(include)
  center <- match.arg(center)

  category <- match.arg(category)
  score.type <- match.arg(score.type)
  signed <- FALSE

  mainloci <- qb.get(qbObject, "mainloci", ...)
  
  tmp <- qb.patterniter(qbObject, epistasis, category, cutoff, mainloci)

  sumpat <- qb.patternave(qbObject, epistasis, tmp$pattern, tmp$nqtl,
                          tmp$pct, mainloci, include, center, level, ...)

  patterns <- attr(sumpat, "pattern")
  n.pat <- length(patterns)
  nqtl <- attr(sumpat, "nqtl")
  pct <- tmp$pct[patterns]
  niter <- names(nqtl)

  ## Split patterns into list. Add attributes. Be carefull about NULL.
  model <- split(data.frame(sumpat), sumpat[,"niter"])
  names(model) <- patterns[patterns != "NULL"]
  tmp <- which.max(pct)[1]
  if(patterns[tmp] == "NULL")
    best <- NULL
  else
    best <- model[[which.max(pct[patterns != "NULL"])]]

  conf <- split(data.frame(attr(sumpat, "conf")), sumpat[, "niter"])
  names(conf) <- patterns[patterns != "NULL"]

  ## Add confidence intervals and drop niter.
  for(i in names(model)) {
    model[[i]] <- cbind(model[[i]][, -1], conf[[i]])
    model[[i]] <- model[[i]][, c("n.qtl","chrom",
                                 "locus","locus.LCL","locus.UCL",
                                 "variance","variance.LCL","variance.UCL")]
  }
  
  ## Compute distances among patterns.
  ## Want to first take care of NULL pattern.
  score <- rep(qb.nulldist(NULL, signed, score.type)$score,
               length(nqtl))
  names(score) <- names(nqtl)
  if(!is.null(best)) {
    tmp <- c(qb.archdist(sumpat, best, FALSE, score.type)$score)
    names(tmp) <- names(nqtl)[nqtl > 0]
    score[names(tmp)] <- tmp
    names(nqtl) <- names(score) <- patterns
  }

  ## Order model by score, then pct.
  model <- model[order(-score[patterns != "NULL"], -pct[patterns != "NULL"])]
  ## Everthing else has its own order.

  out <- list(sumpat = sumpat, patterns = patterns, nqtl = nqtl, pct = pct,
              niter = niter, score = score, model = model)
  attr(out, "signed") <- signed
  attr(out, "score.type") <- score.type
  attr(out, "category") <- category
  class(out) <- c("qb.BestPattern", "list")
  out
}
##########################################################################
qb.patterndist <- function(qbBestObject,
                           score.type = attr(qbBestObject, "score.type"),
                           max.qtl = n.qtl,
                           use = c("complete","pairwise"), ...)
{
  ## Need to use special care with NULL pattern.
  
  use <- match.arg(use)
  sumpat <- qbBestObject$sumpat
  signed <- attr(qbBestObject, "signed")
  
  if(use == "complete") {
    ## Find out maximum number of QTL per chromosome.
    tmp <- levels(sumpat[, "chrom"])
    n.chr <- length(tmp)
    allchr <- rep(0, n.chr)
    names(allchr) <- tmp
    n.qtl <- sum(apply(matrix(unlist(tapply(sumpat[, "chrom"],
                                            sumpat[, "niter"],
                                            function(x, allchr) {
                                              tmp <- table(x)
                                              allchr[names(tmp)] <- c(tmp)
                                              allchr
                                            },
                                            allchr)),
                              length(tmp)),
                       1, max))
  }
  
  ## Dissimilarity of patterns.
  n.pat <- length(qbBestObject$patterns)
  patdist <- rep(0, n.pat * (n.pat - 1) / 2)
  extra <- 0
  if(score.type == "recombination")
    extra <- 0.5
  else if(score.type == "distance" | score.type == "variance")
    extra <- NA

  i.pat <- 0
  for(i in seq(n.pat - 1)) {
    if(qbBestObject$nqtl[i] > 0) {
      tmp <- sumpat[, "niter"] == qbBestObject$niter[i]
      target <- sumpat[tmp,, drop = FALSE]
      if(any(!tmp))
        sumpat <- sumpat[!tmp,, drop = FALSE]
      else
        sumpat <- NULL
    }
    else
      target <- NULL

    if(is.null(sumpat)) {
      dist <- qb.archdist(target, NULL, signed, score.type)
      if(signed)
        dist$score <- - dist$score
    }
    else
      dist <- qb.archdist(sumpat, target, signed, score.type)

    ## Change sum to mean.
    if(score.type != "variance") {
      if(use == "complete") {
        ## Expand number of QTL to be the same for all.
        tmp <- max.qtl - dist$n.union
        if(any(tmp < 0))
          warning(paste("max.qtl too small; increase by", tmp))
        if(score.type == "recombination")
          dist$score <- dist$score + tmp * 0.5
        if(score.type != "distance")
          dist$n.union <- max.qtl
      }
      dist$score <- dist$score / dist$n.union
    }

    n.dist <- length(dist$score)
    patdist[i.pat + seq(n.dist)] <- dist$score
    i.pat <- i.pat + n.dist
  }

  if(score.type == "attenuation" | score.type == "sq.atten")
    patdist <- 1 - patdist

  attr(patdist, "Size") <- n.pat
  attr(patdist, "Labels") <- qbBestObject$patterns
  attr(patdist, "Diag") <- FALSE
  attr(patdist, "Upper") <- FALSE
  attr(patdist, "method") <- "haldane"
  class(patdist) <- "dist"
  out <- list(dist = patdist, nqtl = qbBestObject$nqtl,
              pct = qbBestObject$pct,
              score = qbBestObject$score,
              max.qtl = max.qtl,
              model = qbBestObject$model)
  attr(out, "signed") <- signed
  attr(out, "score.type") <- score.type
  attr(out, "category") <- attr(qbBestObject, "category")
  class(out) <- c("qb.BestPattern", "list")
  out
}
#######################################################################
plot.qb.BestPattern <- function(x, type = c("mds","hclust"),
                              main = paste(score.type, "for",
                                attr(x, "category")),
                              xlab = score.type,
                              method = "complete",
                              cluster = 3,
                              cexmax = 5, colmax = 75,
                              cex = cexs, col = cols,
                              symbol = c("pattern","nqtl","cluster","c@n","c@p","n@p","c@n@p"),
                              ...)
{
  x <- qb.patterndist(x, ...)

  type <- match.arg(type)
  score.type <- attr(x, "score.type")
  
  symbol <- match.arg(symbol)
  colmax <- max(0, min(100, colmax))
  tmp <- (x$pct - min(x$pct)) / (max(x$pct) - min(x$pct))
  cexs <- 1 + (cexmax - 1) * tmp
  cols <- paste("gray", colmax - round(colmax * sqrt(tmp)), sep = "")
  
  hc <- hclust(x$dist, method = method)
  cluster <- cutree(hc, k = cluster)
  symbols <- switch(symbol,
                    pattern = names(x$pct),
                    nqtl = x$nqtl,
                    cluster = cluster,
                    "c@n" = paste(cluster, x$nqtl, sep = "@"),
                    "c@p" = paste(cluster, names(x$pct), sep = "@"),
                    "n@p" = paste(x$nqtl, names(x$pct), sep = "@"),
                    "c@n@p" = paste(cluster, x$nqtl, names(x$pct), sep = "@"))

  switch(type,
         mds = {
           mds <- cmdscale(x$dist, eig = TRUE)
           xmds <- range(mds$points[,1] * (1 + 0.5 * cexmax))
           ymds <- range(mds$points[,2] * (1 + 0.01 * cexmax))
           o <- order(cex)
           if(score.type == "variance") {
             ymds <- range(x$score)
             ymds <- ymds + c(-0.01, 0.01) * cexmax * diff(ymds)
             plot(c(0,2), ymds, type = "n",
                  xlab = "",
                  ylab = "explained variance",
                  main = main, xaxt = "n",
                  ...)
             text(rep(1, length(o)), x$score[o], symbols[o],
                  cex = cex[o], col = col[o])
           }
           else {
             plot(xmds, ymds, type = "n",
                  xlab = paste("MDS Axis 1 (eig=", round(mds$eig[1], 2), ")",
                    sep = ""),
                  ylab = paste("MDS Axis 2 (eig=", round(mds$eig[2], 2), ")",
                    sep = ""),
                  main = main,
                  ...)
             text(mds$points[o, 1], mds$points[o, 2], symbols[o],
                  cex = cex[o], col = col[o])
           }
         },
         hclust = {
           plot(hc, main = main, labels = symbols, xlab = xlab, ...)
         })
}
#######################################################################
print.qb.BestPattern <- function(x, ...) print(summary(x, ...))
#######################################################################
summary.qb.BestPattern <- function(object, method = "complete",
                                 cluster = 3,
                                 n.best = 1, ...)
{
  object <- qb.patterndist(object, ...)
  
  out <- list(max.qtl = object$max.qtl)
  hc <- hclust(object$dist, method = method)

  out$summary <- data.frame(terms = object$nqtl, percent = object$pct,
                            score = object$score,
                            cluster = cutree(hc, k = cluster))[hc$order, ]
  out$summary <- out$summary[order(-out$summary$score, -out$summary$percent), ]

  ## Need to use [patterns != "NULL"] below to get this right.
  pattern.null <- names(object$pct) == "NULL"
  if(n.best == 1) {
    tmp <- which.max(object$score)[1]
    if(pattern.null[tmp])
      out$best <- NULL
    else {
      out$best <- object$model[[which.max(object$score[!pattern.null])]]

      ## Reduce to significant digits.
      out$best[, c("n.qtl","variance","variance.LCL","variance.UCL")] <-
        signif(out$best[, c("n.qtl","variance","variance.LCL","variance.UCL")], 3)
      out$best[, c("locus","locus.LCL","locus.UCL")] <-
        signif(out$best[, c("locus","locus.LCL","locus.UCL")], 5)
    }
  }
  else {
    ## If picking more than one, skip over NULL pattern.
    tmp <- order(-object$score[!pattern.null])[seq(min(n.best, sum(!pattern.null)))]
    out$best <- object$model[tmp]

    ## Reduce to significant digits.
    for(i in names(out$best)) {
      out$best[[i]][, c("n.qtl","variance","variance.LCL","variance.UCL")] <-
        signif(out$best[[i]][, c("n.qtl","variance","variance.LCL","variance.UCL")], 3)
      out$best[[i]][, c("locus","locus.LCL","locus.UCL")] <-
        signif(out$best[[i]][, c("locus","locus.LCL","locus.UCL")], 5)
    }
  }

  out$score.type <- attr(object, "score.type")

  class(out) <- c("summary.qb.BestPattern", "list")
  attr(out, "n.best") <- n.best
  out
}
#######################################################################
print.summary.qb.BestPattern <- function(x, ...)
{
  cat("Best pattern(s) by", x$score.type, "score\n")
  print(x$best)

  cat("\nSummary by better patterns\n")
  print(x$summary)

  cat("\nMaximum number of QTL in architecture:", x$max.qtl, "\n")

  invisible()
}
