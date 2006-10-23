library(qtlbim)

qb.load(cross, qbExample)

x = qb.coda( qbExample, "iterdiag")
summary(x)

## ask before plot routines
tmpar <- par(ask = dev.interactive())
## ask before lattice plots
tmpgrid <- grid::grid.prompt(dev.interactive())
on.exit({
  par(tmpar)
  grid::grid.prompt(tmpgrid)
})

HPDinterval(x,prob=0.95)
effectiveSize(x)
heidel.diag(x)
autocorr.plot(x)
crosscorr.plot(x)
cumuplot(x)
acfplot(x)

  
