#####################################################################
##
## $Id: hpd.R,v 1.7.2.4 2006/12/01 19:59:09 byandell Exp $
##
##     Copyright (C) 2003 Brian S. Yandell
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
qb.hpdHeight <- function(prob, levels = seq(0.5, 0.95, by = 0.05) )
{
  ## HPD region across genome: find probability critical values
  o <- order(-prob)
  p <- cumsum(prob[o])
  p <- p / max(p)

  hpd <- numeric( length( levels ))
  names( hpd ) <- levels
  for (i in seq(along = levels)) {
    tmp <- p <= levels[i]
    if( !any( tmp ))
      tmp <- 1
    hpd[i] <- min(prob[o][tmp])
  }
  hpd
}
##############################################################################
qb.hpdone <- function(qbObject, level = 0.5, profile = "2logBF",
                      effects = "cellmean", scan = "sum",
                      chr = chrs, smooth = 3, ...)
{
  qb.exists(qbObject)
  
  ## Goal: identify regions with posterior above level.

  qbName <- deparse(substitute(qbObject))
  
  if(level < 0 | level > 1)
    stop("level must be between 0 and 1")

  ## Get grid of pseudomarkers.
  grid <- pull.grid(qbObject, offset = TRUE)
  niter <- unclass(table(qb.inter(qbObject, grid, ...)))

  ## Reduce to subset of chromosomes if given in ...
  tmp <- list(...)
  if(!is.null(tmp$chr)) {
    tmp <- !is.na(match(grid$chr, tmp$chr))
    niter <- niter[tmp]
    grid <- grid[tmp, ]
  }
  scan <- scan[1]
  one <- qb.scanone(qbObject, ..., type.scan = "posterior")
  tmp <- names(one)[-(1:2)]
  if(is.na(match(scan, tmp)))
    scan <- tmp[1]

  grid[[scan]] <- qb.smoothone(one[, scan], grid, smooth, niter)

  ## Find HPD height and reduce grid to those at or above height.
  hpd <- qb.hpdHeight(grid[[scan]], level)
  grid <- grid[grid[[scan]] >= hpd, ]
  wh <- tapply(grid[[scan]], grid$chr, which.max)
  len <- cumsum(table(grid$chr))
  wh <- wh + c(0, len[-length(len)])
  pos <- grid$pos[wh]
  peak <- grid[[scan]][wh]

  ## Set up matrix with limits.
  if(nrow(grid)) {
    prob <- unlist(tapply(grid[[scan]], grid$chr, sum))
    grid <- tapply(grid$pos, grid$chr, range)
    gridnames <- list(names(grid),
                      c("chr", "n.qtl", "pos",
                        paste(c("lo", "hi"), ".",
                              round(100 * level), "%", sep = "")))
    chrs <- as.numeric(names(grid))
    grid <- cbind(chrs, prob, pos,
                  matrix(unlist(grid), ncol = 2, byrow = TRUE))
    dimnames(grid) <- gridnames

    ## Subset to selected chr.
    chr <- chr[!is.na(match(chr, row.names(grid)))]
    chr.hpd <- as.character(chr)
    grid <- grid[chr.hpd, ]
    if(length(chr) == 1) {
      grid <- matrix(grid, 1)
      dimnames(grid) <- list(chr.hpd, gridnames[[2]])
    }
    geno.names <- qb.geno.names(qbObject)
    grid <- data.frame(grid)
    grid$chr <- ordered(geno.names[grid$chr], geno.names)
    rownames(grid) <- as.character(grid$chr)
    out <- list(hpd.region = grid)

    out$profile <- qb.scanone(qbObject, type.scan = profile, chr = chr, ...)
    attr(out$profile, "qb") <- qbName
    out$effects <- if(effects == "cellmean")
      qb.scanone(qbObject, type.scan = "cellmean", chr = chr, ...)
    else
      qb.scanone(qbObject, type.scan = "estimate",
                    scan = "main", aggregate = FALSE, chr = chr, ...)
    attr(out$effects, "qb") <- qbName

    class(out) <- c("qb.hpdone", "matrix")
    attr(out, "profile") <- profile
    attr(out, "effects") <- effects
    attr(out, "hpd.height") <- hpd
    attr(out, "hpd.level") <- level
    attr(out, "map.range") <- lapply(pull.map(qb.cross(qbObject, genoprob = FALSE)), range)
    out

    ## Might consider separate lines at epistasis and main?
    ## Would need to write qb.get.main to complement qb.get.epis
    ## in summary.qb.scanone. Best left for later iteration.
  }
  else
    NULL
}
###################################################################
print.qb.hpdone <- function(x, ...) print(summary(x, ...))
###################################################################
summary.qb.hpdone <- function(object, chr = chrs, digits = 3, ...)
{
  if(is.null(object))
    stop("No HPD regions identified (level to small?)")

  ## The selected chr is subset of row.names of hpd.region.
  chrs <- row.names(object$hpd.region)
  chr <- chr[!is.na(match(chr, row.names(object$hpd.region)))]
  chr.hpd <- as.character(chr)

  tmp <- object$hpd.region[chr.hpd,, drop = FALSE]

  ## Add peak of profile.
  tmp2 <- summary(object$profile)
  tmp2 <- tmp2[match(chr, tmp2[, "chr"]), "sum", drop = FALSE]

  ## Add effects estimates.
  tmp3 <- summary(object$effects)
  tmp3 <- tmp3[match(chr, tmp3[, "chr"]), -(1:4), drop = FALSE]

  tmp <- cbind(tmp, tmp2, tmp3)
  tmprow <- dimnames(tmp3)[[2]]
  dimnames(tmp) <- list(chr,
                        c(dimnames(object$hpd.region)[[2]],
                          attr(object, "profile"),
                          tmprow))
  class(tmp) <- c("summary.qb.hpdone", "data.frame")
  attr(tmp, "digits") <- digits
  tmp
}
###################################################################
print.summary.qb.hpdone <- function(x, ...)
{
  x[, -1] <- round(x[, -1], attr(x, "digits"))
  class(x) <- "data.frame"
  NextMethod()
}
###################################################################
plot.qb.hpdone <- function(x, chr = chrs, ...)
{
  if(is.null(x))
    stop("No HPD regions identified (level to small?)")

  ## The selected chr is subset of row.names of hpd.region.
  chrs <- x$hpd.region$chr
  chr <- x$hpd.region$chr[match(chr,
                                {
                                  if(is.numeric(chr))
                                    unclass(x$hpd.region$chr)
                                  else
                                    x$hpd.region$chr
                                },
                                nomatch = 0)]
  if(!length(chr))
    stop("no chromosomes to plot")
  chr.hpd <- as.character(chr)

  if(length(chr) > 1) {
    ## Expand map limits to add lines to scanone plot.
    map <- attr(x, "map.range")[chr]
    minmap <- unlist(lapply(map, min))
    maxmap <- unlist(lapply(map, max))
    map <- c(-25, maxmap - minmap)
    map <- cumsum(25 + map)[-length(map)] - minmap
  }
  else
    ## With one chr, scanone starts at 0 now.
    map <- 0

  tmpar <- par(mfrow = c(2,1), mar = c(4.1,4.1,3.1,0.1))
  on.exit(par(tmpar))

  ## NB: scan takes numeric chr for now.
  plot(x$profile, chr = chr, ...)
  tmp <- pmatch(c("lo.","hi."), dimnames(x$hpd.region)[[2]])
  apply(cbind(x$hpd.region[chr.hpd, tmp[1]] + map,
              x$hpd.region[chr.hpd, tmp[2]] + map), 1,
        function(x) lines(x, rep(0,2), col = "red", lwd = 3))

  ## Vertical line at peak.
  apply(as.matrix(map + x$hpd.region[chr.hpd,"pos"]), 1, function(x)
        abline(v = x, lty = 2, lwd = 2, col = "gray"))

  plot(x$effects, chr = chr, ...)
  apply(as.matrix(map + x$hpd.region[chr.hpd,"pos"]), 1, function(x)
        abline(v = x, lty = 2, lwd = 2, col = "gray"))
  invisible()
}
#######################################################################################3
qb.hpdchr <- function(qbObject, level = 0.5, height = hpd.height,
                      hpd = qb.hpdone(qbObject, level = level, profile = "post"), 
                      chr = hpd$hpd.region$chr,
                      smooth = 3)
{
  if(!missing(height)) {
    level <- 1
  }
  if(attr(hpd, "profile") != "post")
    stop("qb.hpdone must be run with profile = \"post\"")
  hpd.height <- attr(hpd,"hpd.height")
  out <- array(NA, length(chr))
  names(out) <- chr
  grid <- pull.grid(subset(qbObject, chr = unclass(chr)))
  for(i in seq(length(chr))) {
    tmp <- hpd$profile[grid$chr == unclass(chr[i]), "sum"]
    tmp <- qb.smoothchr(tmp, smooth, tmp)
    tmp2 <- tmp > height
    if(any(tmp2))
      out[as.character(chr[i])] <- sum(tmp[tmp2]) / sum(tmp)
  }
  out <- 100 * out[!is.na(out)]
  list(hpd.height = height, chr.posterior = out)
}
