#####################################################################
##
## $Id: qb.covar.tour.R,v 1.3.2.8 2006/12/06 15:14:39 byandell Exp $
##
##     Copyright (C) 2006 Brian S. Yandell
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
## Make sure MCMC samples exist first.
library(qtlbim)

data(qbExample)
cross <- qb.cross(qbExample)

## Box and whiskers plot of covariate.
bwplot(pheno.normal~fix.cov,cross$pheno,
       horizontal = FALSE, jitter.data = TRUE)

## XY plot of random covariate.
gbye.geno <- factor(pull.geno(cross)[, find.marker(cross, 2, 12)])
xyplot(pheno.normal~random.cov,cross$pheno, type = c("p","r"),
       group = gbye.geno,
       key = simpleKey(levels(gbye.geno), space = "right"))

## Show grand mean vs. covariate estimates with 5 percentiles.
## Default for covar is all covariates.
tmp <- qb.meancomp(qbExample)
summary(tmp)
plot(tmp)

## Show main effects vs. GxE effect estimates.
## Default for covar is first covariate.
tmp <- qb.covar(qbExample)
plot(tmp)
summary(tmp)

## Show correlation of covariate with markers.
## Default for covar is first covariate.
tmp <- qb.confound(qbExample)
summary(tmp)
plot(tmp)

## Variance components, including random covariate.
tmp <- qb.varcomp(qbExample)
summary(tmp)
plot(tmp)

## Effects for interacting covariates.
tmp <- qb.intcov(qbExample)
summary(tmp)
plot(tmp)

## Cell mean scan for interacting covariates.
tmpar <- par(mfrow=c(2,1))
tmp <- qb.scanone(qbExample, type = "cellmean")
plot(tmp)
title("\n\nfemales")
## Lines at QTLs.
abline(v=c(15,45,60+25+12,2*(60+25)+15),lty=2)
tmp <- qb.scanone(qbExample, type = "cellmean", adjust.covar = 1)
plot(tmp)
title("\n\nmales")
## Lines at QTLs.
abline(v=c(15,45,60+25+12,2*(60+25)+15),lty=2)
par(tmpar)
