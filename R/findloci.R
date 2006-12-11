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
qb.locus <- function(qbObject, chr = 1,
  iterdiag = qb.get(qbObject, "iterdiag"),
  mainloci = qb.get(qbObject, "mainloci"))
{
  n.iter <- nrow(iterdiag)
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
  iterdiag <- qb.get(qbObject, "iterdiag")
  mainloci <- qb.get(qbObject, "mainloci")
  counts <- matrix(0, length(chr), qb.get(qbObject, "max.nqtl"))
  chrs <- as.character(chr) ## could get names from cross$geno
  dimnames(counts) <- list(chrs, as.character(seq(ncol(counts)) - 1))
  for(i in seq(chr)) {
    tmp <- qb.locus(qbObject, chr[i], iterdiag, mainloci)
    counts[chrs[i], names(tmp)] <- tmp
  }
  tmp <- apply(counts, 2, function(x) any(x > 0))
  tmp[1] <- TRUE ## always include 0 column
  counts <- round(100 * counts[, tmp] / nrow(iterdiag))
  counts[apply(counts[, -1], 1, sum) >= threshold, ]
}
##############################################################
qb.multloci <- function(qbObject, chr = 1)
{
  qb.exists(qbObject)
  
  chr <- chr[1]
  tmp <- qb.locus(qbObject, chr)
  tmp <- round(100 * tmp / sum(tmp), 1)
  res <- list(nqtl = tmp)
  tmp <- qb.get(qbObject, "mainloci")[, c("niter", "chrom", "locus")]
  tmp <- tmp[tmp[, "chrom"] == chr, c("niter", "locus")]
  res$mainloci <- tmp
  tmp <- qb.get(qbObject, "pairloci")[, c("chrom1", "chrom2", "locus1", "locus2")]
  tmp <- tmp[tmp[, "chrom1"] == chr & tmp[, "chrom2"] == chr, c("locus1", "locus2")]
  res$pairloci <- tmp
  class(res) <- c("qb.multloci", "list")
  attr(res, "chr") <- chr
  attr(res, "cex") <- qb.cex(qbObject)
  res
}
##############################################################
summary.qb.multloci <- function(object, ...)
{
  object$nqtl
}
##############################################################
print.qb.multloci <- function(x, ...) print(summary(x, ...))
##############################################################
plot.qb.multloci <- function(x, amount = .75,
                             cex = attr(x, "cex"), ...)
{
  chr <- attr(x, "chr")
  require("lattice")
  trellis.par.set(theme = col.whitebg(), warn = FALSE) ## white background

  is.two <- length(x$nqtl) > 2
  print(densityplot(~locus, x$mainloci, cex = cex,
                    main = paste("all loci on chr", chr)),
        more = TRUE, split = c(1 + is.two, 1 + is.two, 2, 1 + is.two))

  tmp <- data.frame(Percent = c(x$nqtl), nqtl = names(x$nqtl))
  print(barchart(Percent ~ nqtl, tmp,
                 main = paste("number of QTL on chr", chr)),
        more = is.two, split = c(2, 1, 2, 1 + is.two))
  
  if(is.two) {
    tmp <- data.frame(locus = unlist(x$pairloci))
    print(densityplot(~locus, tmp, cex = cex, horizontal = TRUE,
                      main = paste("epistatic loci on chr", chr)),
          more = TRUE, split=c(1, 1, 2, 2))
    tmp <- tapply(x$mainloci[, "locus"], x$mainloci[, "niter"],
      function(x) {
        nx <- length(x)
        if(nx == 1)
          NULL
        else {
          rx <- matrix(seq(nx), nx, nx)
          cx <- col(rx)
          rbind(x[c(rx[rx > cx])], x[c(cx[rx > cx])])
        }
      })
    tmp <- matrix(unlist(tmp), 2)
    tmp <- apply(tmp, 1, jitter, amount = amount)
    tmp <- as.data.frame(tmp)

    tmp2 <- apply(x$pairloci, 2, jitter, amount = amount)
    tmp2 <- as.data.frame(tmp2)
    names(tmp) <- names(tmp2) <- c("locus1","locus2")
    tmp <- rbind(tmp,tmp2)
    print(xyplot(locus2 ~ locus1, tmp, cex = cex,
                 panel = function(x,y,...) {
                   panel.abline(0,1)
                   panel.xyplot(x,y,...)
                   },
                 main = paste("epistatic/all pairs on chr", chr)),
          split = c(1,2,2,2))
  }
  invisible()
}
