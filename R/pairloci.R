#####################################################################
## This code has been written by Dr.Yandell.
## $Id: pairloci.R,v 1.7.2.2 2006/10/02 19:18:53 byandell Exp $
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
qb.pair.posterior <- function(qbObject, cutoff = 1, nmax = 15)
{
  pairloci <- qb.get(qbObject, "pairloci")
  if(is.null(pairloci)) {
    cat("no epistatic pairs\n")
    return(invisible(NULL))
  }
  geno.names <- names(qb.cross(qbObject)$geno)
  percent <- 100 *
    rev(sort(table(interaction(geno.names[pairloci[, "chrom1"]],
                               geno.names[pairloci[, "chrom2"]])))) /
                                 nrow(qb.get(qbObject, "iterdiag"))
  percent <- percent[ percent >  cutoff ]
  if(length(percent) > nmax)
    percent <- percent[seq(nmax)]
  round(percent)
}
##############################################################################
qb.pair.nqtl <- function(qbObject, cutoff = 1,
  pairloci = qb.get(qbObject, "pairloci"))
{
  if(is.null(pairloci)) {
    cat("no epistatic pairs\n")
    return(invisible(NULL))
  }
  pairs <- pairloci$n.epis[ !duplicated(pairloci$niter) ]
  m <- mean(pairs)
  l <- length(pairs)
  pairs <- 100 * table(pairs) / l
  round(pairs[ pairs > cutoff ])
}
##############################################################################
qb.pairloci <- function(qbObject, chr)
{
  pairloci <- qb.get(qbObject, "pairloci")
  if(is.null(pairloci)) {
    cat("no epistatic pairs\n")
    return(invisible(NULL))
  }

  ## Find epistatic pairs of loci for chr.
  chr <- sort(chr)
  if(length(chr) == 1)
    chr <- rep(chr, 2)
  else
    chr <- chr[1:2]

  d <- pairloci[pairloci$chrom1 == chr[1] & pairloci$chrom2 == chr[2],
                c("locus1", "locus2")]
  names(d) <- paste("chr", chr, sep = ".")
  d <- apply(d, 2, jitter, amount = 0.5)
  
  class(d) <- c("qb.pairloci", "matrix")
  attr(d, "chr") <- chr
  attr(d, "niter") <- nrow(qb.get(qbObject, "iterdiag"))
  attr(d, "map") <- pull.map(qb.cross(qbObject))[chr]
  attr(d, "post") <- qb.pair.posterior(qbObject)[paste(chr, collapse = ".")]
  d
}
##############################################################################
summary.qb.pairloci <- function(object, ...)
{
  rbind(apply(object,2,quantile, c(.25,.5,.75)),
        samples=nrow(object),
        percent = round(100 * nrow(object) / attr(object, "niter"), 3))
}
##############################################################################
print.qb.pairloci <- function(x, ...) print(summary(x, ...))
##############################################################################
plot.qb.pairloci <- function(x, main = mainchr, cex = 0.75, ...)
{
  map <- attr(x, "map")
  chr <- attr(x, "chr")
  post <- attr(x, "post")

  class(x) <- "matrix"
  ## qb starts chromosome from 0, but map may be positive
  rng <- list()
  for(r in 1:2) {
    rng[[r]] <- range(map[[r]])
    x[,r] <- x[,r] + rng[[r]][1]
  }
  
  plot(x, cex = cex, xlim = rng[[1]], ylim = rng[[2]])
  for(r in 1:2)
    axis(r, map[[r]], lwd = 3, labels = FALSE)
  abline(v = median(x[, 1]), h = median(x[, 2]),
         lty = 2, lwd = 3, col = "blue")

  ## Plots only in uppertriangle if chr are the same.
  if(chr[1] == chr[2])
    abline(0, 1)

  ## Add title.
  mainchr <- paste(chr[1], " by ", chr[2]," (", post, "%)", sep = "")
  title(main)
  invisible()
}
##############################################################################
qb.epistasis <- function(qbObject, effects = c("aa","ad","da","dd"),
                         cutoff = 1, maxpair = 5, pairs = names(post), ...)
{
  pairloci <- qb.get(qbObject, "pairloci")
  if(is.null(pairloci)) {
    cat("no epistatic pairs\n")
    return(invisible(NULL))
  }

  ## Identify pairs of chromosomes with interacting QTL.
  geno.names <- names(qb.cross(qbObject)$geno)
  inter <- interaction(geno.names[pairloci[, "chrom1"]],
                       geno.names[pairloci[, "chrom2"]])
  post <- qb.pair.posterior(qbObject, cutoff)
  if(length(post) > maxpair)
    post <- post[seq(maxpair)]
  if(!is.character(pairs))
    stop("pairs must be character")
  pairs <- pairs[match(names(post), pairs, nomatch = 0)]
  post <- post[pairs]
  inter.pairs <- !is.na(match(inter, pairs))

  ## Subset on desired Cockerham effects.
  effects <- effects[ match(names(pairloci), effects, nomatch = 0) ]
  pairloci <- as.data.frame(as.matrix(pairloci[inter.pairs, effects]))
  names(pairloci) <- effects

  ## Kludge inter to add posterior to name of pairs.
  inter <- ordered(as.character(inter[inter.pairs]), pairs)
  pairs.pct <- paste(pairs, "\n", post, "%", sep = "")

  pairloci$inter <- ordered(pairs.pct[unclass(inter)], pairs.pct)

  class(pairloci) <- c("qb.epistasis", "data.frame")
  attr(pairloci, "post") <- post
  pairloci
}
##############################################################################
summary.qb.epistasis <- function(object, ...)
{
  nc <- ncol(object)

  tmp <- as.matrix(object[, -nc])
  if(ncol(tmp) == 1)
    dimnames(tmp) <- list(NULL, names(object)[-nc])
  
  signif(cbind("%" =  attr(object, "post"),
               apply(tmp, 2,
                     function(x, y) tapply(x, y, median),
                     object[, nc])),
         3)
}
##############################################################################
print.qb.epistasis <- function(x, ...) print(summary(x, ...))
##############################################################################
plot.qb.epistasis <- function(x, effects = c("aa","ad","da","dd"),
                              cex = 0.5, ...)
{
  require("lattice")
  trellis.par.set(theme=col.whitebg(), warn = FALSE)

  effects <- effects[ match(names(x), effects, nomatch = 0) ]
  for(j in seq(along = effects)) {
    form <- formula(paste(effects[j], "inter", sep = "~"))
    print(bwplot(form, x, jitter = TRUE, factor = 1, cex = cex,
                 col = "gray75", lwd = 2,
                 xlab = "epistatic pair", ylab = "epistatic effect",
                 panel = function(x, y, ...) {
                   panel.abline(h = 0, lwd = 3, lty = 2, col = "red")
                   panel.stripplot(x,y,...)
                   lx <- levels(x)
                   if(sum(y != 0) > 100)
                     panel.bwplot(x, y, ..., do.out = FALSE)
                   for(i in seq(along = lx)) {
                     ii <- (x == lx[i])
                     if(any(ii))
                       panel.lines(i+c(-0.25,0.25),
                                   rep(median(y[ii], na.rm = TRUE), 2),
                                   lwd = 5, col = "blue")
                   }
                 },
                 horizontal = FALSE,
                 main = effects[j]),
          split = c(j, 1, length(effects), 1),
          more = (j < length(effects)))
  }
}
##############################################################################
qb.chrom <- function(qbObject)
{
  geno.names <- names(qb.cross(qbObject)$geno)
  chrom <- c(table(geno.names[qb.get(qbObject, "mainloci")$chrom]))
  maplen <- unlist(lapply(pull.map(qb.cross(qbObject)),
                           function(x) diff(range(x))))
  niter <- nrow(qb.get(qbObject, "iterdiag"))
  ## caution: posterior does not account for duplicate chromosomes
  assess <- data.frame(posterior = chrom / sum(chrom),
    prior = maplen[names(chrom)] / sum(maplen))
  assess$bf <- assess$posterior / assess$prior
  bf <- min(assess$bf)
  if(bf > 0 & min(assess$prior) > 0)
    assess$bf <- assess$bf / bf
  assess$bfse <- assess$bf * 
    sqrt((1 - assess$posterior) / (assess$posterior * niter))
  assess
}
##############################################################################
qb.pairs <- function(qbObject, cutoff = 1, nmax = 15)
{
  pairloci <- qb.get(qbObject, "pairloci")
  if(is.null(pairloci)) {
    cat("no epistatic pairs\n")
    return(invisible(NULL))
  }
  npair <- qb.pair.nqtl(qbObject, cutoff, pairloci)
  niter <- nrow(qb.get(qbObject, "iterdiag"))
  geno.names <- names(qb.cross(qbObject)$geno)
  inter <- interaction(geno.names[pairloci[, "chrom1"]],
                       geno.names[pairloci[, "chrom2"]])
  posterior <- rev(sort(table(inter))) / niter
  posterior <- posterior[ posterior > cutoff / 100 ]
  posterior[posterior > 1] <- 1
  if(!length(posterior))
    return(NULL)
  if(length(posterior) > nmax)
    posterior <- posterior[1:nmax]

  maplen <- unlist(lapply(pull.map(qb.cross(qbObject)),
                          function(x)diff(range(x))))
  i <- !duplicated(inter)
  prior <- maplen[pairloci$chrom1[i]] * maplen[pairloci$chrom2[i]]
  names(prior) <- inter[i]
  prior <- prior / sum(prior)
  prior <- prior[match(names(posterior), names(prior), nomatch = 0)]
  assess <- data.frame(posterior = posterior, prior = prior)
  assess$bf <- assess$posterior / assess$prior
  bf <- min(assess$bf)
  if(bf > 0 & min(assess$prior) > 0)
    assess$bf <- assess$bf / bf
  assess$bfse <- assess$bf * 
    sqrt((1 - assess$posterior) / (assess$posterior * niter))
  assess
}
##############################################################################
plot.qb.pattern <- function (x, bars = seq(x$posterior),
  labels = c("model index", "model posterior", "pattern posterior"),
  barlabels = names(x$posterior), 
  threshold = c(weak = 3, moderate = 10, strong = 30), units = 2, 
  rescale = TRUE, ...) 
{
  bar <- barplot(c(x$posterior), col = "white", names = bars, ...)
  tmp <- if (rescale) 
    x$prior * 0.99 * max(x$posterior) / max(x$prior)
  else x$prior
  lines(bar, tmp, type = "b", col = "blue", lwd = 2)
  if (!is.null(barlabels)) {
    cex <- 1
#    usr <- par("usr")
#    tmp <- as.numeric(2 * (x$posterior - usr[3]) >= diff(usr[3:4]))
#    for (i in 0:1) {
#      ii <- i == tmp
#      if (any(ii)) 
#        text(bar[ii], x$posterior[ii], barlabels[ii], 
#             srt = 90, adj = i, cex = cex)
#    }
    text(bar, 0, barlabels, srt = 90, adj = 0, cex = cex)
    
  }
  mtext(labels[1], 1, 2)
  mtext(labels[2], 2, 2)
  mtext(labels[3], 3, 0.5)
  x$bf[x$bf == 0 | x$prior == 0] <- NA
  plot(seq(bars), x$bf, log = "y", xaxt = "n",
       xlim = c(0.5, length(bars)), xlab = "", ylab = "", ...)
  mtext(labels[1], 1, 2)
  mtext("posterior / prior", 2, 2)
  mtext("Bayes factor ratios", 3, 0.5)
  axis(1, seq(bars), bars)
  if (!is.null(barlabels)) {
    cxy <- par("cxy")[1] * cex/2
    nqtl <- strsplit(barlabels, "")
    nqtl <- lapply(nqtl, function(x) {
      colon <- seq(x)[x == ":" | x == "*"][1]
      x <- if (is.na(colon)) 
        "1"
      else
        x[seq(colon - 1)]
      paste(x, collapse = "")
    })
    text(seq(barlabels) - cxy, x$bf, nqtl, srt = 90, cex = cex)
  }
  usr <- 10^par("usr")[3:4]
  for (i in seq(bars)) {
    if (x$bfse[i] > 0) {
      bfbar <- x$bf[i] + c(-units, units) * x$bfse[i]
      bfbar[1] <- max(usr[1], bfbar[1])
      bfbar[2] <- min(usr[2], bfbar[2])
    }
    else bfbar <- usr
    lines(rep(i, 2), bfbar)
  }
  if (length(threshold)) {
    bars <- floor(mean(seq(bars))/2) + 0.5
    maxusr <- usr[2]
    usr <- prod(usr^c(0.95, 0.05))
    lines(bars + c(-0.25, 0.25), rep(usr, 2), lwd = 3, col = "blue")
    texusr <- usr
    for (i in seq(length(threshold))) {
      sigusr <- min(maxusr, usr * threshold[i])
      if (texusr < maxusr) 
        text(bars + 0.5, sqrt(texusr * sigusr), names(threshold)[i], 
             col = "blue", adj = 0)
      arrows(bars, usr, bars, sigusr, 0.1, lwd = 3, col = "blue")
      texusr <- sigusr
    }
  }
  x
}
