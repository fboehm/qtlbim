####################################################################################

##     Copyright (C) 2006 Nengjun Yi and Tapan Mehta
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

####################################################################################


### NOTE: Depends on R/qtl 1.03.
qb.genoprob <- function(cross, map.function = map.functions, step = 2,
                        tolerance = 1e-6, stepwidth = "variable", ...)
{
  ## This uses calc.genoprob directly.
  ## Call sequence maintained for legacy user code.
  map.functions <- unlist(as.list(formals(calc.genoprob)$map.function)[-1])
  map.function <- match.arg(tolower(map.function[1]), map.function)

  ## First make sure QTL are not too close.
  if(any(sapply(pull.map(cross), function(x)
                ifelse(length(x) > 1, max(diff(x)) < tolerance, FALSE))))
    cross <- jittermap(cross, tolerance)
  
  ## Call R/qtl routine calc.genoprob for calculations.
  cross <- calc.genoprob(cross, step, map.function = map.function,
                         stepwidth = stepwidth,
                         ...)

  ## Add tolerance as attribute.
  attr(cross$geno[[1]]$prob, "tolerance") <- tolerance
  
  return(cross)
}

