#####################################################################
##
## $Id: qb.hyper.tour.R,v 1.3.2.4 2006/09/07 01:55:21 byandell Exp $
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
## A tour of the hyper data with possibly some new insights.
library(qtlbim)

data(qbHyper)

## Standard summaries.
summary(qbHyper)
plot(qbHyper)

########################################################
## 1-D scans.

one <- qb.scanone(qbHyper, type = "LPD")
summary(one)
plot(one)

sum.one <- summary(one, sort = "sum",
                   threshold = c(sum = 4, epistasis = 2))
sum.one
chrs <- sort(unique(sum.one$chr))
chrs

## Look at LPD on key chromosomes.
plot(one, chr = chrs)

## Look at cell means.
onemean <- qb.scanone(qbHyper, type = "cellmean")
plot(onemean, chr = chrs)

########################################################
## 2-D scans.

two <- qb.scantwo(qbHyper, chr = chrs, type = "LPD")
summary(two, sort = "upper")
## Focus on epistatic chromosome pairs.
plot(two, chr = c(4,6,7,15))

## Slice for chromosome 15.
plot(two, chr = c(4,6,7), slice = 15)

## 2-D estimates of main and epistatic effects.
twoest <- qb.scantwo(qbHyper, chr = chrs, type = "estimate",
  scan = list(lower="add",upper="aa"))

## Slice for epistatic effect for chr 15.
plot(twoest, chr = c(4,6,7), slice = 15)

## Slice for cell mean.
slice <- qb.sliceone(qbHyper, type = "cellmean", chr = c(4,6,7), slice = 15)
summary(slice)
plot(slice, chr = 6:7)
plot(slice, chr = 6:7, scan = "slice")
