#####################################################################
##
## $Id: zzz.R,v 1.3 2006/08/08 02:34:33 byandell Exp $
##
## .onLoad is run when the package is loaded with library(qtlbim)
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

.onLoad <- function(lib,pkg)
  {
  # We no longer need "library.dynam("qtlbim",pkg,lib) because this
  # is taken care of by the NAMESPACE directive useDynLib(qtlbim).

  # The R documentation suggests that "require(qtl)" should be
  # replaced by an "import(qtl)" directive in the NAMESPACE file.
  # However package::qtl does not define a name space.
  require(qtl)
  
  }


# end of zzz.R

