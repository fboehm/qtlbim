#####################################################################
##
## $Id: qb.plot.tour.R,v 1.5.2.5 2006/10/23 17:11:54 byandell Exp $
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

## summary of MCMC run
summary(qbExample)

## Default plots.
plot(qbExample)

## Model selection via Bayes factor.
plot(qb.BayesFactor(qbExample))

## Show effects and HPD region.
plot(qb.hpdone(qbExample))

