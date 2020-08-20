#####################################################################
##
## $Id: gbye.R,v 1.7.2.8 2006/12/06 15:26:31 byandell Exp $
##
##     Copyright (C) 2008 Brian S. Yandell
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
qb.gbye.posterior <- function(qbObject, covar = get.covar, cutoff = 1, nmax = 5,
                              gbye = qb.get(qbObject, "gbye", ...), ...)
{
  if(is.null(gbye))
    return(invisible(NULL))

  cross <- qb.cross(qbObject, genoprob = FALSE)
  geno.names <- qb.geno.names(qbObject, cross)
  get.covar <- qb.get(qbObject,
                      "nfixcov")[as.logical(qb.get(qbObject, "intcov"))]
  covar.names <- names(cross$pheno)[covar]

  percent <- table(covar.names[gbye$covar], geno.names[gbye$chrom])
  np <- dimnames(percent)
  np <- outer(np[[2]], np[[1]], paste, sep = ":")
  percent <- c(percent)
  o <- order(-percent)
  percent <- 100 * percent[o] / qb.niter(qbObject)
  names(percent) <- np[o]
  percent <- percent[ percent >  cutoff ]
  if(length(percent) > nmax)
    percent <- percent[seq(nmax)]
  round(percent)
}
##############################################################################
qb.intcov <- function(qbObject, covar = get.covar, effects = c("add","dom"),
                         cutoff = 1, nmax = 5, cov.chr = names(post), ...)
{
  qb.exists(qbObject)
  
  ## Use qb.gbye.posterior to find interesting pairs.
  ## Set up object as in qb.epistasis.
  ## Set up plot or use plot.qb.epistasis.

  gbye <- qb.get(qbObject, "gbye", ...)
  if(is.null(gbye)) {
    cat("no GxE\n")
    return(invisible(NULL))
  }

  ## Identify covariates and chromosomes with interacting QTL.
  cross <- qb.cross(qbObject, genoprob = FALSE)
  geno.names <- qb.geno.names(qbObject, cross)
  get.covar <- qb.get(qbObject,
                      "nfixcov")[as.logical(qb.get(qbObject, "intcov"))]
  if(!length(covar))
    return(invisible(NULL))
  if(all(covar <= 0))
    return(invisible(NULL))

  covar.names <- names(cross$pheno)[covar]

  inter <- paste(geno.names[gbye[, "chrom"]],
                 covar.names[gbye[, "covar"]], sep = ":")
  post <- qb.gbye.posterior(qbObject, covar, cutoff, nmax, gbye = gbye)
  if(!is.character(cov.chr))
    stop("cov.chr must be character")
  cov.chr <- cov.chr[match(names(post), cov.chr, nomatch = 0)]
  post <- post[cov.chr]
  inter.cov.chr <- !is.na(match(inter, cov.chr))

  ## Subset on desired Cockerham effects.
  effects <- effects[match(names(gbye), effects, nomatch = 0)]
  gbye <- as.data.frame(as.matrix(gbye[inter.cov.chr, effects]))
  names(gbye) <- effects

  ## Kludge inter to add posterior to name of cov.chr.
  inter <- ordered(as.character(inter[inter.cov.chr]), cov.chr)
  cov.chr.pct <- paste(cov.chr, "\n", post, "%", sep = "")

  gbye$inter <- ordered(cov.chr.pct[unclass(inter)], cov.chr.pct)

  class(gbye) <- c("qb.epistasis", "data.frame")
  attr(gbye, "post") <- post
  gbye
}

