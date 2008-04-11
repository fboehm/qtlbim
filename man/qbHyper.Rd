\name{qbHyper}
\alias{qbHyper}
\docType{data}
\title{MCMC samples for hyper dataset from R/qtl}
\description{
Retrieve MCMC samples for hyper dataset from R/qtl. Code provided to
recreate qbHyper if desired. These samples are used in vignettes and demos.
}
\usage{data(qbHyper)}
\seealso{
\code{\link[qtl]{hyper}}, \code{\link{qb.genoprob}}, \code{\link{qb.mcmc}}
}
\examples{
data(qbHyper)
summary(qbHyper)
\dontrun{
## Here is code to generate qbHyper.

## Get data from R/qtl.
data(hyper)

## Restrict to autosomes.
hyper <- subset(hyper, chr = 1:19)

## Calculate genotype probabilities.
hyper <- qb.genoprob(hyper, step=2)

## Create MCMC samples.
qbHyper <- qb.mcmc(hyper, n.thin = 40, seed = 1616)

## The next line saves qbHyper as an external binary file.
save("qbHyper", file = "qbHyper.RData")
}
}
\keyword{datasets}
