#####################################################################
##
## $Id: qb.coda.tour.R,v 1.3.2.8 2006/12/06 15:14:39 byandell Exp $
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
library(qtlbim)

data(qbExample)

x = qb.coda(qbExample, "iterdiag")
summary(x)

HPDinterval(x,prob=0.95)
effectiveSize(x)
heidel.diag(x)
autocorr.plot(x)
crosscorr.plot(x)
cumuplot(x)
acfplot(x)

  
