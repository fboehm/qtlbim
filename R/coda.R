####################################################################################

##     Copyright (C) 2006 Nengjun Yi 
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


qb.coda <- function(object,
                    element = c("iterdiag","mainloci","pairloci","covariates","gbye"),
                    variables = variable.list[[element]])
{
  qb.exists(object)
  
   require("coda")
   variable.list <- list(iterdiag = c("nqtl","mean","envvar","var"),
                         mainloci = c("chrom","add","dom"),
                         pairloci = c("chrom","aa","ad","da","dd"),
                         covariates = "cov1",
                         gbye = c("n.gbye","chrom","add","dom"))
   element <- match.arg(element)
   tmp <- qb.get(object, element)
   variables <- names(tmp)[match(variables, names(tmp), nomatch = 0)]
   if(length(variables))
      mcmc(tmp[, variables])
}



