\name{qb.remove}
\alias{qb.remove}
\alias{qb.recover}
\title{Remove or recover qb object and associated MCMC samples}
\description{
This removes the object and the directory that contains the MCMC
samples. It is not enough to remove the R object.
}
\usage{
qb.remove(qbObject, verbose = TRUE)
qb.recover(cross, traitName, output.dir, n.thin = 40, n.burnin,
  algorithm = "M-H", genoupdate = FALSE, \dots)
}
\arguments{
  \item{qbObject}{Object of class \code{qb} (see \code{\link{qb.mcmc}}).}
  \item{verbose}{Print warning if \code{TRUE}.}
  \item{cross}{Object of class \code{cross} (see \code{\link[qtl]{read.cross}}).}
  \item{traitName}{Character string name of trait to recover.}
  \item{output.dir}{Character string with name of output directory
    (inferred if missing).}
  \item{n.thin}{Thining of MCMC chain used in \code{\link{qb.mcmc}}.}
  \item{n.burnin}{Burnin of MCMC chain used in \code{\link{qb.mcmc}}.}
  \item{algorithm}{Algorithm of MCMC chain used in \code{\link{qb.mcmc}}.}
  \item{genoupdate}{Genotype update flag for MCMC chain used in \code{\link{qb.mcmc}}.}
  \item{\dots}{Options passed to \code{\link{qb.data}} and
    \code{\link{qb.model}}.}
}
\details{
At the present time, \code{\link{qb.mcmc}} stores MCMC samples in
external files located in directory \code{output.dir}, whose name is
typically the \code{traitName} followed by the date. \code{qb.remove}
removes this directory along with the
\code{qbObject}. \code{qb.recover} attempts to recover the use of an
orphaned \code{output.dir} after a crash of R. These are fragile routines.
}
\references{http://www.qtlbim.org}
\author{Brian S. Yandell}
\seealso{\code{\link{qb.mcmc}}, \code{\link{qb.genoprob}}}
\examples{
\dontrun{
## Recover qbExample for trait "bp" of cross "hyper" using default output.dir.
qbExample <- qb.recover(hyper, "bp")

## Remove internal qbExample and external output.dir.
qb.remove(qbExample)
}
}
\keyword{utilities}
