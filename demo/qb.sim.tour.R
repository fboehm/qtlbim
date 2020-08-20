#####################################################################
##
## $Id: qb.sim.tour.R,v 1.3.2.3 2006/09/06 18:43:28 byandell Exp $
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
library(qtlbim)

sim.compare <- function(n.ind = 100, effect = 1, n.qtl = 2,
                        n.mark = 200, bymark = 1, pos = n.mark / 2,
                        ... )
{
  ## Simulate cross data with 1 QTL.
  markers <- seq(0, n.mark, by = bymark)
  names(markers) <- paste("M",markers,sep="")
  sim <- sim.cross(list(ch1 = markers),
                   matrix(c(1, pos, effect), 1, 3),
                   n.ind,
                   "bc" )

  ## Classical R/qtl.
  sim <- calc.genoprob(sim, 1, error = 0)
  sim.scan <- scanone(sim)

  ## Bayesian R/qtlbim.
  sim <- qb.genoprob(sim, step=1)
  assign("sim", sim, 1)

  qbSim <- qb.mcmc(sim, chr.nqtl = n.qtl, main.nqtl = 1, mean.nqtl = 1,
                   max.nqtl = n.qtl, epistasis = FALSE)

  bone.lod <- qb.scanone(qbSim, type = "LPD")
  bone.add <- qb.scanone(qbSim, type="estimate", scan="add")
  ## Kludge to match effectscan.
#  bone.add$add <- -bone.add$add
  
  cat("Classical R/qtl summary\n")
  print(summary(sim.scan))
  cat("Bayesian R/qtlbim summary\n")
  print(summary(bone.lod))
  print(summary(bone.add))

  ## Idealized confounded effect.
  poss <- seq(0, n.mark, by = 1)
  expect <- effect * exp(-2 * abs(poss - pos) / 100)
  rss <- (n.ind - 2) * (1 + ((effect)^2 - expect^2) / 4)
  lod <- n.ind * log10(((n.ind - 1) * (1 + (effect)^2 / 4)) / rss) / 2
 
  ## LOD scans.
  cat("Plot of LOD scans\n")
  tmpar <- par( mfrow = c(2,1), mar = c(5.1,4.1,2.1,0.1) )
  on.exit(par(tmpar))

  ylim <- range(bone.lod$add, sim.scan$lod, lod)
  plot(sim.scan, ylim = ylim, lty = 2, ...)
  mtext("black dash = qtl; red solid = qtlbim; blue dot = ideal 1 QTL",
        1, 2)
  plot(bone.lod, add = TRUE, col = "red")
  lines(poss, lod, lty = 3, lwd = 2, col = "blue")

  ## Substitution effect scans.
  ## Classical one QTL.
  cat("Plot of effect scans\n")
  
  sim <- sim.geno(sim, n.draws = 1)
  tmp <- effectscan(sim, draw = FALSE, ...)
  ylim <- range(bone.add$add, 0, effect, tmp$a)
  plot(tmp, ylim = ylim, lty = 2, lwd = 2, add.legend = FALSE, ...)
  ## Want line color to be black, so redraw.
  lines(tmp$pos, tmp$a, lty = 2, lwd = 2, col = "black")
  axis( 1, at = markers, labels = FALSE )
  mtext("black dash = qtl; red solid = qtlbim; blue dot = ideal 1 QTL",
        1, 2)

  ## Bayesian 1-2 QTL.
  plot(bone.add, add = TRUE, col = "red")

  lines(poss, expect, lty = 3, lwd = 2, col = "blue")
  invisible(list(sim.data = sim, sim.qb = qbSim))
}
#########################################################

## Large effect simulation.
sim.compare(n.ind = 100, effect = 1, n.qtl = 2)
## Smaller effect simulation.
sim.compare(n.ind = 200, effect = 0.5, n.qtl = 2)
