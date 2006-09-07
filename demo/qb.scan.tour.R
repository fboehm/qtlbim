#####################################################################
##
## $Id: qb.scan.tour.R,v 1.2.2.4 2006/09/07 01:55:21 byandell Exp $
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
## create 1-D marginal scan
## diagnostics for each locus conditional on all other QTL
library(qtlbim)

qb.load(cross, qbExample)

one <- qb.scanone(qbExample, type = "LOD")

summary(one)
sum.one <- summary(one, threshold = c(main = 5, epistasis = 5, GxE = 5),
                   order = "sum")
sum.one
chrs <- as.vector(sum.one[,"chr"])

## ask before plot routines
tmpar <- par(ask = dev.interactive())
on.exit({
  par(tmpar)
})

## show conditional contribution to LOD across whole genome
plot(one, smooth = 3)
## focus on heritability for 
plot(one, chr = chrs, smooth = 3)
## show only main effects LOD
one <- qb.scanone(qbExample, type = "LOD", aggregate = FALSE)
plot(one, chr = chrs, scan = "main", smooth = 3,
     col=c(add.fix.cov="turquoise", dom.fix.cov="magenta"))

## 2-D scan conditional on all other QTL
two <- qb.scantwo(qbExample, chr = chrs, type = "LOD")

sum.two <- summary(two, threshold = c(upper=5), sort = "upper")
sum.two
chr2 <- as.matrix(sum.two[,1:2])
## nearest neighbor smoothing of order 3
plot(two, smooth = 3)

plot(two, chr=chr2[1,], smooth = 3)

####################################################3
## Subset on regions of chromosomes 1,3,5.

qb2 <- subset(qbExample, chr = chr2[1,],
               region = data.frame(
                 chr = chr2[1,],
                 start = rep(0,2),
                 end = rep(40,2)))

two <- qb.scantwo(qb2, type = "LOD")
plot(two, smooth = 3)

####################################################################
## Variance components: add and aa.
two <- qb.scantwo(qb2, type = "var",
                  scan = list(lower = "add", upper = "aa"))
plot(two, smooth = 3)
## Heritability: main and ad.
two <- qb.scantwo(qb2, type = "her",
                  scan = list(lower = "main", upper = "ad"))
plot(two, smooth = 3)

plot(qb.scanone(qbExample, chr = chrs, type = "detect"),
     smooth = 3)
