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
qb.genoprob <- function(cross, map.func = c("Haldane","Kosambi"), step = 2,
                        tolerance = 1e-6, ...)
{
  ## This uses new variable-width create.map function
  ## within call to calc.genoprob. That is, we no longer
  ## need to call our internal C code.
  map.functions <- c("haldane","kosambi","c-f","morgan")
  map.func <- pmatch(tolower(map.func), map.functions, nomatch = 1)
  map.func <- map.functions[map.func[1]]
  
  ## First make sure QTL are not too close.
  if(any(sapply(pull.map(cross), function(x) max(diff(x))) < tolerance))
    cross <- jittermap(cross, tolerance)
  
  ## Call R/qtl routine calc.genoprob for calculations.
    cross <- calc.genoprob(cross, step, map.function = map.func,
                           stepwidth = "variable",
                           ...)
  return(cross)
}

