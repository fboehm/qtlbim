#####################################################################
##
## $Id: qb.scan.tour.R,v 1.2.2.6 2006/12/05 20:12:23 byandell Exp $
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

data(qbExample)

## Conditional contribution to LPD across whole genome.
one <- qb.scanone(qbExample, type = "LPD")
summary(one)
plot(one)

## Show only main effects LPD.
one <- qb.scanone(qbExample, type = "LPD", aggregate = FALSE)
plot(one, scan = "add")

## 2-D scan conditional on all other QTL
two <- qb.scantwo(qbExample, type = "LPD")
plot(two)

####################################################3
## Subset on regions of chromosomes 1,2.
qbSubset <- subset(qbExample, chr = c(1,2),
  region = data.frame(chr = c(1,2), start = c(35,2), end = c(55,22)))

two <- qb.scantwo(qbSubset, type = "LPD")
plot(two)

####################################################################
## Variance components: add and aa.
two <- qb.scantwo(qbSubset, type = "var",
                  scan = list(lower = "add", upper = "aa"))
plot(two)
## Heritability: main and ad.
two <- qb.scantwo(qbSubset, type = "her",
                  scan = list(lower = "main", upper = "ad"))
plot(two)

plot(qb.scanone(qbExample, type = "detect"))
