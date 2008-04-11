#####################################################################
##
## $Id: findloci.R,v 1.3.2.2 2006/10/24 15:22:26 byandell Exp $
##
##     Copyright (C) 2002 Brian S. Yandell
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


##############################################################
qb.locus <- function(qbObject, chr = 1, n.iter = qb.niter(qbObject),
  mainloci = qb.get(qbObject, "mainloci"))
{
  mainloci <- qb.get(qbObject, "mainloci")
  count <- table(tapply(mainloci[, "chrom"], mainloci[, "niter"], 
    function(x, y) sum(x == y),
    chr))
  if(names(count)[1] != "0")
    count <- c("0" = n.iter - sum(count), count)
  count
}
##############################################################
qb.loci <- function(qbObject, chr = sort(unique(mainloci[, "chrom"])),
  threshold = 25)
{
  n.iter <- qb.niter(qbObject)
  mainloci <- qb.get(qbObject, "mainloci")
  counts <- matrix(0, length(chr), qb.get(qbObject, "max.nqtl"))
  chrs <- as.character(chr) ## could get names from cross$geno
  dimnames(counts) <- list(chrs, as.character(seq(ncol(counts)) - 1))
  for(i in seq(chr)) {
    tmp <- qb.locus(qbObject, chr[i], n.iter, mainloci)
    counts[chrs[i], names(tmp)] <- tmp
  }
  tmp <- apply(counts, 2, function(x) any(x > 0))
  tmp[1] <- TRUE ## always include 0 column
  counts <- round(100 * counts[, tmp] / n.iter)
  counts[apply(counts[, -1], 1, sum) >= threshold, ]
}
##############################################################
qb.multloci <- function(qbObject, chr = 1, cutoff = 25, nqtl = NULL, ...)
{
  qb.exists(qbObject)
  
  chr <- chr[1]

  ## Get main loci and extract modes.
  res <- qb.mainmodes(subset(qbObject, chr = chr), cutoff, nqtl, ...)

  ## Main loci on chromosome chr.
  tmp <- qb.get(qbObject, "mainloci")[, c("niter", "chrom", "locus")]
  tmp <- tmp[tmp[, "chrom"] == chr, c("niter", "locus")]
  res$mainloci <- tmp

  ## Epistatic pairs with both loci on chr.
  tmp <- qb.get(qbObject, "pairloci")
  tmp <- tmp[, c("niter", "chrom1", "chrom2", "locus1", "locus2")]
  tmp <- tmp[tmp[, "chrom1"] == chr & tmp[, "chrom2"] == chr,
             c("niter", "locus1", "locus2")]
  res$pairloci <- tmp
  class(res) <- c("qb.multloci", "qb.mainmodes", "list")
  attr(res, "chr") <- chr
  attr(res, "cex") <- qb.cex(qbObject)
  attr(res, "step") <- qb.get(qbObject, "step")
  attr(res, "pos") <- pull.map(qb.cross(qbObject, genoprob = FALSE))[[chr]]
  attr(res, "niter") <- qb.niter(qbObject)
  res
}
##############################################################
print.qb.multloci <- function(x, ...) print(summary(x, ...))
##############################################################
summary.qb.multloci <- function(object, merge = TRUE, ...)
{
  out <- object[1:4]
  niter <- attr(object, "niter")

  locus.more <- object$mainloci
  locus.more$qtl <- apply(outer(locus.more$locus, object$valleys[[1]],
                                function(x,y) x > y),
                          1, sum) + 1
  locus.more$nqtl <-
    c(table(object$mainloci$niter)[as.character(object$mainloci$niter)])
  my.summary <- function(x, niter) {
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
  
  out$qtl <- if(merge | length(object$nqtl.post[[1]]) == 1)
    my.summary(locus.more, niter)
  else {
    tmp <- split(locus.more, locus.more$nqtl)
    tmp2 <- c("=",">=")[1 + (as.numeric(names(tmp)) > c(object$nqtl.est))]
    names(tmp) <- paste("nqtl", tmp2, names(tmp))
    lapply(tmp, my.summary, niter)
  }

  class(out) <- c("summary.qb.multloci", "list")
  out
}
##############################################################
print.summary.qb.multloci <- function(x, ...)
{
  cat("Posterior Percent by Number of QTL\n")
  print(x$nqtl.post[[1]])
  cat("Estimated Number of QTL:", as.numeric(x$nqtl.est), "\n\n")
  cat("Peaks\n")
  print(x$peaks[[1]])
  cat("Valleys\n")
  print(x$valleys[[1]])
  cat("QTL Summaries\n")
  print(x$qtl)
  invisible()
}
##############################################################
plot.qb.multloci <- function(x, amount = .5, cex = attr(x, "cex"),
                             split = TRUE, contour = TRUE,
                             weight = TRUE,
                             merge = TRUE,
                             ...)
{
  chr <- attr(x, "chr")
  pos <- attr(x, "pos")
  require("lattice")
  trellis.par.set(theme = col.whitebg(), warn = FALSE) ## white background

  ## Split can be logical or numeric.
  ## If numeric, plot only those pages.
  if(is.logical(split)) {
    if(merge)
      show.plot <- 1:4
    else 
      show.plot <- 1
  }
  else {
    if(!is.numeric(split))
      stop("split must be logical or numeric")
    show.plot <- pmax(1, pmin(4, split))
    split <- FALSE
  }
  
  is.two <- as.numeric(x$nqtl.post[[1]])
  is.two <- max(is.two) > 1
  
  nqtl.more <- c(x$nqtl.est > 1)
  if(nqtl.more)
    valleys <- x$valleys[[1]]
  else
    valleys <- min(x$mainloci$locus) - 1
  amount <- amount * attr(x, "step")

  if(is.two) {
    table.niter <- table(x$mainloci$niter)
    if(!is.null(x$pairloci))
      mainloci.nqtl <- table.niter[as.character(x$pairloci$niter)]
  }

  locus.more <- data.frame(locus = x$mainloci$locus)
  if(any(show.plot == 1)) {
    if(is.two) {
      ## Density plot of QTL on chr; group by estimated QTL.
      locus.more$qtl <- apply(outer(locus.more$locus, valleys,
                                    function(x,y) x > y),
                              1, sum) + 1
      locus.more$locus <- jitter(locus.more$locus, amount = amount)
      locus.more$nqtl <- c(table.niter[as.character(x$mainloci$niter)])
      if(nqtl.more) {
        if(weight) {
          weights <- 1 / table.niter[as.character(x$mainloci$niter)]
          weights <- weights / sum(weights)
        dens <- density(locus.more$locus, weights = weights)
          darg <- list(weights = weights)
        }
        else {
          dens <- density(locus.more$locus)
          darg <- list()
        }
      }
      ## My version of density and weighted density, adding annotation.
      ## Assumes response = locus, group = qtl, number of QTL in nqtl.
      mydensityplot <- function(locus.more, dens, main) {
        if(weight) {
          densityplot(~locus, locus.more, cex = cex, groups = qtl, lty = 1,
                      main = main,
                      panel = function(x, ...) {
                        ## Get trellis plot info.
                        line.info <- trellis.par.get("superpose.line")
                        symbol.info <- trellis.par.get("superpose.symbol")
                        jitter.amount <-
                          0.01 * diff(current.panel.limits()$ylim)
                        uqtl <- sort(unique(locus.more$qtl))
                        for(i in seq(length(uqtl))) {
                          ii <- (uqtl[i] == locus.more$qtl)
                          weights <- 1 / locus.more$nqtl[ii]
                          weights <- weights / sum(weights)
                          panel.lines(density(locus.more$locus[ii],
                                              weights = weights),
                                      col = line.info$col[i],
                                      lty = 1)
                          panel.xyplot(x = locus.more$locus[ii],
                                       y = jitter(rep(0,
                                         length(locus.more$locus[ii])),
                                         amount = jitter.amount),
                                       pch = symbol.info$pch[i],
                                       col = symbol.info$col[i],
                                       ...)
                        }
                        if(nqtl.more) {
                          panel.abline(v = valleys, lty = 2,
                                       col = "darkgray")
                          panel.lines(dens$x,
                                      (1 + length(valleys)) * dens$y,
                                      col = "black", lty = 3)
                        }
                        panel.rug(x = pos, end = 0.02, quiet = TRUE,
                                  col = "black")
                      })
        }
        else {
          densityplot(~locus, locus.more, cex = cex, groups = qtl, lty = 1,
                      main = main,
                      panel = function(x, ...) {
                        panel.superpose(x, ...)
                        if(nqtl.more) {
                          panel.abline(v = valleys, lty = 2, col = "darkgray")
                          panel.lines(dens$x, (1 + length(valleys)) * dens$y,
                                      col = "black", lty = 3)
                        }
                        panel.rug(x = pos, end = 0.02, quiet = TRUE, col = "black")
                      },
                      panel.groups = "panel.densityplot")
        }
      }
      if(merge) {
        print(mydensityplot(locus.more, dens,
                            paste("all loci by QTL on chr", chr)),
              more = split & TRUE,
              split = if(split) c(2, 2, 2, 2) else c(1,1,1,1))
      }
      else {
        tmp <- c("=",">=")[1 + (locus.more$nqtl > c(x$nqtl.est))]
        print(densityplot(~locus | factor(paste("n.qtl", tmp, locus.more$nqtl)),
                          locus.more, cex = cex, groups = qtl, lty = 1,
                          panel = function(x, ...) {
                            panel.densityplot(x, ...)
                            panel.rug(x = pos, end = 0.02, quiet = TRUE, col = "black")
                          },
                          main = paste("all loci by QTL on chr", chr),
                          layout = c(1, length(unique(locus.more$nqtl)))))
      }
    }
    else {
      ## Density plot of QTL on chr; at most 1 QTL per iteration.
      locus.more$locus <- jitter(locus.more$locus, amount = amount)
      print(densityplot(~locus, locus.more, cex = cex,
                        main = paste("all loci on chr", chr),
                        panel = function(...) {
                          panel.densityplot(...)
                          panel.rug(x = pos, end = 0.02, quiet = TRUE,
                                    col = "black")
                        }),
            more = split & TRUE,
            split = if(split) c(1, 1, 2, 1) else c(1,1,1,1))
    }
  }
    
  ## Bar chart of number of QTL.
  if(any(show.plot == 2)) {
    print(barchart(Percent ~ nqtl,
                   data.frame(Percent = c(x$nqtl.post[[1]]),
                              nqtl = names(x$nqtl.post[[1]])),
                   main = paste("number of QTL on chr", chr)),
          more = split & is.two,
          split = if(split) c(2, 1, 2, 1 + is.two) else c(1,1,1,1))
  }
  if(is.two) {
    if(any(show.plot == 3) & !is.null(x$pairloci)) {
      ## Density plot of epistatic pairs; group by estimated QTL.
      locus.more <-
        data.frame(locus = unlist(x$pairloci[, c("locus1", "locus2")]))
      locus.more$qtl <- apply(outer(locus.more$locus, valleys,
                                    function(x,y) x > y),
                              1, sum) + 1
      locus.more$locus <- jitter(locus.more$locus, amount = amount)
      locus.more$nqtl <-
        table(x$pairloci$niter)[as.character(x$pairloci$niter)]
      if(nqtl.more)
        dens <- density(locus.more$locus)
      if(nqtl.more) {
        if(weight) {
          weights <- 1 / locus.more$nqtl
          weights <- weights / sum(weights)
          dens <- density(locus.more$locus, weights = weights)
        }
        else
          dens <- density(locus.more$locus)
      }
      print(mydensityplot(locus.more, dens,
                          paste("epi loci by QTL on chr", chr)),
            more = split & TRUE,
            split = if(split) c(1, 1, 2, 2) else c(1,1,1,1))
    }
    if(any(show.plot == 4)) {
      ## Scatter plot of pairs: lower triangle = all pairs.
      locus.more <- tapply(x$mainloci[, "locus"], x$mainloci[, "niter"],
                           function(x) {
                             nx <- length(x)
                             if(nx == 1)
                               NULL
                             else {
                               rx <- matrix(seq(nx), nx, nx)
                               cx <- col(rx)
                               rbind(x[c(rx[rx > cx])],
                                     x[c(cx[rx > cx])],
                                     nx)
                             }
                           })
      locus.more <- matrix(unlist(locus.more), 3)
      locus.more <- as.data.frame(t(locus.more))
      names(locus.more) <- c("locus1", "locus2","nqtl")
      loci <- paste(locus.more$locus1, locus.more$locus2, sep = ":")
      locus.more$point <- table(loci)[loci]
      
      if(!is.null(x$pairloci)) {
        ## Scatter plot of pairs: upper triangle = epistatic pairs.
        locus.more <- rbind(locus.more,
                            cbind(x$pairloci[, c("locus1", "locus2")],
                                  nqtl = mainloci.nqtl, point = 0))
        rm(mainloci.nqtl)
      }
      
      ## Add diagonal with solo QTL.
      table.niter <- table.niter[as.character(x$mainloci$niter)]
      locus.one <- x$mainloci$locus[table.niter == 1]
      rm(table.niter)
      gc()
      
      ## Print scatter plot.
      if(contour) {
        ## Overlay contour lines.
        loci <- sort(unique(c(locus.more$locus1, locus.more$locus2)))
        grid <- expand.grid(list(locus1 = loci, locus2 = loci))
        if(weight) {
          weights <- 1 / (locus.more$nqtl - 1)
          weights <- weights / sum(weights)
          contours <- loess(point ~ locus1 * locus2, locus.more, weights)
        }
        else
          contours <- loess(point ~ locus1 * locus2, locus.more)
        contours <- c(predict(contours, grid))
        ## Force zero on epistasis side.
        contours[grid$locus1 <= grid$locus2] <- 0
        contours <- contourLines(loci, loci,
                                 matrix(contours, length(loci), length(loci)),
                                 levels = pretty(max(contours, na.rm = TRUE) * c(0.6, 1)))
        rm(loci)
      }
        
      ## Jitter loci across pseudomarker interval.
      for(i in c("locus1","locus2"))
        locus.more[[i]] <- jitter(locus.more[[i]], amount = amount)
      
      ## XY plot with optional contour lines.
      tmp <- range(locus.one, locus.more$locus1, locus.more$locus2)
      print(xyplot(locus2 ~ locus1, locus.more, cex = cex, groups = nqtl,
                   xlim = tmp, ylim = tmp,
                   panel = function(x,y,...) {
                     panel.abline(0,1, col = "gray")
                     panel.abline(v = valleys, h = valleys,
                                  lty = 2, col = "darkgray")
                     if(length(locus.one))
                       panel.points(jitter(locus.one, amount = amount),
                                    jitter(locus.one, amount = amount),
                                    col = "black", cex = cex)
                     panel.superpose(x,y, ...)
                     ## Could use contourLines here instead.
                     panel.rug(x = pos, end = 0.02, quiet = TRUE)
                     panel.rug(y = pos, end = 0.02, quiet = TRUE)
                     if(contour)
                       lapply(contours, function(l)
                              panel.lines(l$x, l$y, col = "black"))
                   },
                   panel.groups = "panel.xyplot",
                   main = paste("epi/all pairs by nqtl on chr", chr)),
            split = if(split) c(1, 2, 2, 2) else c(1,1,1,1))
    }
  }
  invisible()
}
