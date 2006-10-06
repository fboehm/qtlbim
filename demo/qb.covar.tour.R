#####################################################################
##
## $Id: qb.covar.tour.R,v 1.3.2.5 2006/10/06 15:42:45 byandell Exp $
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

qb.load(cross, qbExample)

## Set white background for lattice graphics.
trellis.par.set(theme = col.whitebg(), warn = FALSE)
## Prompt for trellis (same as par(ask=TRUE) for plot).
tmpgrid <- grid::grid.prompt(dev.interactive())
on.exit({
  grid::grid.prompt(tmpgrid)
})

## Box and whiskers plot of covariate.
bwplot(pheno.normal~fix.cov,cross$pheno,
       horizontal = FALSE, jitter.data = TRUE)

## XY plot of random covariate.
gbye.geno <- factor(pull.geno(cross)[, find.marker(cross, 12, 35)])
xyplot(pheno.normal~random.cov,cross$pheno, type = c("p","r"),
       group = gbye.geno,
       key = simpleKey(levels(gbye.geno), space = "right"))

## Show grand mean vs. covariate estimates with 5 percentiles.
## Default for covar is all covariates.
tmp <- qb.meancomp(qbExample)
summary(tmp)
plot(tmp)

if(!exists("qbsub"))
  qbsub <- subset(qbExample, pattern = c(1,3,5,12), chr = c(1,3,5,12))

## Show main effects vs. GxE effect estimates.
## Default for covar is first covariate.
tmp <- qb.covar(qbsub)
plot(tmp)
summary(tmp)
plot(qb.covar(qbsub, element = "dom"))

## Show correlation of covariate with markers.
## Default for covar is first covariate.
tmp <- qb.confound(qbsub)
summary(tmp)
plot(tmp)

## Variance components, including random covariate.
tmp <- qb.varcomp(qbsub)
summary(tmp)
plot(tmp)

## Effects for interacting covariates.
tmp <- qb.intcov(qbsub)
summary(tmp)
plot(tmp)

## Cell mean scan for interacting covariates.
tmpar <- par(mfrow=c(3,1))
tmp <- qb.scanone(qbsub, type = "cellmean")
plot(tmp)
abline(v=35,lty=2)
tmp <- qb.scanone(qbsub, type = "cellmean", adjust.covar = 1)
plot(tmp)
abline(v=35,lty=2)
tmp <- qb.scanone(qbsub, type = "cellmean", adjust.covar = 1)
plot(tmp)
abline(v=35,lty=2)
par(tmpar)
