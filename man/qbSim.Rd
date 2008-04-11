\name{qbSim}
\alias{qbSim}
\alias{qbSimMain}
\alias{qbSimEpi}
\docType{data}
\title{MCMC samples used in scan.pdf document}
\description{
  Retrieve or recreate MCMC samples used in scan.pdf document.
}
\usage{
data(qbSimMain)
data(qbSimEpi)
}
\details{
  Both calls to \code{data} create \code{qb} objects names \code{qbSim}.
  See vignette \code{scan.pdf} or see \code{scan.Rnw} in doc folder of package.
}
\seealso{
\code{\link{qb.genoprob}}, \code{\link{qb.mcmc}}, \code{\link{qb.sim.cross}}
}
\examples{
data(qbSimMain)
summary(qbSim)

data(qbSimEpi)
summary(qbSim)

\dontrun{
## Setup for Simulated Data used in scan.pdf.
n.ind <- 100		## number of individuals
n.mark <- 200		## number of markers
by.mark <- 1		## cM spacing between markers
qtl.positions <- n.mark / 2	## position of QTL
markers <- seq(0, n.mark, by = by.mark)
names(markers) <- paste("M", markers, sep = "")
sim.map <- list(ch1 = markers)
sim.model <- matrix(c(1, qtl.positions, qtl.effect / 2), 1, 3)
colnames(sim.model) <- c("chromosome","qtl-position","effect-size")
n.iter <- 1000		## number of iterations for MCMC
qb.random.seed <- 1626  ## random seed for MCMC

## Genetic architecture for scan simulations: 3 QTL.
qtl.positions <- rbind(qtl1 = c(chromosome = 1, locus = 5), 
                       qtl2 = c(chromosome = 1, locus = 50), 
                       qtl3 = c(chromosome = 2, locus = 33) )

qtl.positions
qtl.main.model <-
  rbind(qtl1.main.effect = c(qtl = 1, main.effect.size = 0),
        qtl2.main.effect = c(qtl = 2, main.effect.size = 0),
        qtl3.main.effect = c(qtl = 3, main.effect.size = 0))

qtl.main.model
qtl.epi.model <- rbind(qtl1.and.qtl3.epi.effect =
  c(qtl1 = 1, qtl2 = 3, epi.effect.size = 10))
qtl.epi.model

## SimEpi
set.seed(1234)
sim <- qb.sim.cross(len = rep(100, 2), n.mar = 10, eq.spacing = TRUE, 
                    n.ind = 100, type = "bc", missing.geno = 0.03,
                    qtl.pos = qtl.positions, 
                    qtl.main = qtl.main.model,
                    qtl.epis = qtl.epi.model)
sim <- qb.genoprob(sim)
qbSim <- qb.mcmc(sim, n.iter = n.iter, verbose = FALSE, n.thin = 40,
                 seed = qb.random.seed)

## The next line saves qbSim as an external binary file.
save("qbSim", file = "qbSimEpi.RData")

## SimMain
qtl.main.model[2, "main.effect.size"] = 10
set.seed(1234)
sim <- qb.sim.cross(len = rep(100, 2), n.mar = 10, eq.spacing = TRUE, 
                    n.ind = 100, type = "bc", missing.geno = 0.03,
                    qtl.pos = qtl.positions, 
                    qtl.main = qtl.main.model,
                    qtl.epis = NULL)

## After the data is simulated call qb.genoprob to fill in 
## missing data.
sim <- qb.genoprob(sim, step = 2)

## Call qb.mcmc and then analysis code.
qbSim <- qb.mcmc(sim, n.iter = n.iter, verbose = FALSE, n.thin = 40,
                 seed = qb.random.seed)

## The next line saves qbSim as an external binary file.
save("qbSim", file = "qbSimMain.RData")
}
}
\keyword{datasets}
