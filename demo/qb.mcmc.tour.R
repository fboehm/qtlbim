library(qtlbim)

## F2 example.
cross <- qb.sim.cross(len = rep(100,20), n.mar = 11, eq.spacing = F,
                      n.ind = 500, type = "f2", ordinal = c(0.3,0.3,0.2,0.2),
                      missing.geno = 0.03, missing.pheno = 0.03,
                      qtl.pos = rbind(qtl.1 = c(chr=1,pos=15),
                        qtl.2 = c(1,45), qtl.3 = c(3,12),
                        qtl.4 = c(5,15), qtl.5 = c(7,15),
                        qtl.6 = c(10,15), qtl.7 = c(12,35),
                        qtl.8 = c(19,15)),
                      qtl.main = rbind(main.1 = c(qtl=1,add=0.5,dom=0),
                        main.2 = c(2,0,0.7), main3 = c(3,-0.5,0),
                        main4 = c(4,0.5,-0.5)),
                      qtl.epis = rbind(epis1 = c(qtl.a=4,qtl.b=5,aa=-0.7,
                                         ad=0,da=0,dd=0),
                        epis2 = c(6,8,0,1.2,0,0)),
                      covariate = c(fix.cov=0.5,ran.cov = 0.07),
                      gbye = rbind(GxE.1 = c(qtl=7,add=0.8,dom=0)) ) 

## Summary of simulation information.
summary(cross$qtl)

## Compute genotype probabilities and run MCMC.
cross <- qb.genoprob(cross, step=2) 

## Prepare data for MCMC.
qbData <- qb.data(cross, pheno.col = 3, rancov = 2, fixcov = 1)

## Set up interacting QTL model for qb.mcmc
qbModel <- qb.model(cross, chr.nqtl = rep(3,nchr(cross)),
                    prop = c(0.5, 1, 0.5), intcov = 1)

## Create MCMC samples.
qbExample <- qb.mcmc(cross, data = qbData, model = qbModel,
                     n.iter = 3000, n.thin = 10, n.burnin = 0,
                     genoupdate = F)
