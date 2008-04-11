\name{qb.get}
\alias{covar.mean}
\alias{pull.grid}
\alias{pull.loci}
\alias{qb.get}
\alias{qb.cross}
\alias{qb.cross.class}
\alias{qb.demo}
\alias{qb.load}
\alias{qb.save}
\alias{qb.reorder}
\title{Internal qtlbim routines}
\description{
  These are internal qtlbim routines that are made visible in the
  namespace for technical use.
}
\usage{
covar.mean(qbObject, adjust.covar, verbose = FALSE)
pull.grid(qbObject, offset, spacing, mask.region, cross, step,
  off.end, stepwidth, \dots)
pull.loci(cross, step, off.end, stepwidth, region)
qb.cross(qbObject, genoprob = TRUE, \dots)
qb.cross.class(qbObject)
qb.demo()
qb.get(qbObject, element, pheno.id = 1, warn = TRUE, \dots)
qb.load(cross, qbObject, dir, file)
qb.save(cross, qbObject, dir, Name)
qb.reorder(qbObject, warn = FALSE)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{adjust.covar}{Adjustments to covariates. Default is
    \code{NA}, which adjusts by covariate mean values. Values are
    assumed to be in order of fixed covariates.}
  \item{verbose}{Verbose mode if \code{TRUE}.}
  \item{cross}{Object of class \code{cross} (see \code{\link[qtl]{read.cross}}).}
  \item{offset}{Offset by first marker if \code{TRUE}.}
  \item{step}{Argument to \code{\link[qtl]{create.map}}.}
  \item{off.end}{Argument to \code{\link[qtl]{create.map}}.}
  \item{stepwidth}{Argument to \code{\link[qtl]{create.map}}.}
  \item{region}{Argument to \code{\link[qtl]{create.map}}.}
  \item{genoprob}{Compute genotype probabilities with
    \code{\link{qb.genoprob}} if \code{TRUE}.}
  \item{spacing}{Add columns for \code{map}, \code{eq.spacing} and
    \code{xchr} if \code{TRUE}. This corresponds to \code{map} element
    of a \code{\link[qtl]{scantwo}} object.}
  \item{mask.region}{Subset genome regions if \code{TRUE} (see
    \code{\link{subset.qb}}).}
  \item{element}{Character string for element of \code{qbObject} to
    get. Typically this is a parameter to \code{\link{qb.data}},
    \code{\link{qb.model}} or \code{\link{qb.mcmc}}, or it is one of the
    MCMC sample files in \code{output.dir}, from \code{c("iterdiag",
      "mainloci", "pairloci", "covariates", "gbye")}.}
  \item{pheno.id}{Numeric or character identifier for phenotype. Useful
    eventually for multiple traits.}
  \item{warn}{Warn user if \code{qbObject} is legacy format if \code{TRUE}.}
  \item{dir}{Character string name of directory for load if \code{qbObject}
    does not exist.}
  \item{file}{Character string name of file for load if \code{qbObject}
    does not exist.}
  \item{Name}{Character string name of suffix for save if \code{qbObject}
    exists.}
  \item{\dots}{Additional arguments.}
}
\details{
  These are all internal routines. But some may be useful beyond.

  \code{qb.demo} is called in \code{demo(qb.tour)} and provides an
  interactive selection of the \code{R/qtlbim} demos.

  \code{qb.cross} extracts the \code{cross} object associated with
  \code{qbObject}. \code{qb.get} is the internal main routine for
  extracting information from a \code{qbObject}. As stated elsewhere,
  currently \code{qbObject} refers to objects that are critical to it
  but not part of it: the \code{cross} object used to create it and the
  MCMC samples in files in \code{output.dir}.

  \code{covar.mean} finds covariate means or adjusts them to
  user-supplied values.
  
  \code{pull.grid} pulls the grid of pseudomarkers from the \code{cross}
  object associated with \code{qbObject}. The option \code{spacing}
  determines whether this is in a format similar to
  \code{\link[qtl]{scanone}} (\code{FALSE}) or
  \code{\link[qtl]{scantwo}} (\code{TRUE}). It is used \code{qb.get}
  when accessing external MCMC sample files and by several other
  routines that require pseudomarker information, notably genotype
  probabilities.

  \code{qb.reorder} is called by \code{\link{qb.mcmc}} to create
  pointers to reorder the MCMC samples so that chromosome numbers and
  positions within chromosomes are in increasing order. It creates the
  \code{subset} element of a \code{qb} object.

  \code{qb.save} and \code{qb.load} are depricated. They were used to
  save and retrieve MCMC 
  samples used in examples and vignettes. We now use
  \code{\link[base]{save}} and \code{\link[base]{data}}.
}
\examples{
data(qbExample)

covar.mean(qbExample)
qb.get(qbExample, "output.dir")
summary(qb.cross(qbExample))
temp <- qb.get(qbExample, "iterdiag")
dim(temp)
names(temp)

\dontrun{
## The following should have no effect.
qbExample <- qb.reorder(qbExampleb)

## You can call the following rather than demo() to get a tour.
qb.demo()
}
}
\keyword{utilities}
