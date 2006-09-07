\name{qb.get}
\alias{pull.grid}
\alias{qb.get}
\alias{qb.cross}
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
pull.grid(qbObject, offset, spacing, mask.region)
qb.cross(qbObject)
qb.demo()
qb.get(qbObject, element, sub)
qb.load(cross, qbObject, dir, file)
qb.save(cross, qbObject, dir, file)
qb.reorder(qbObject)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{cross}{Object of class \code{cross} (see \code{\link[qtl]{read.cross}}).}
  \item{offset}{Offset by first marker if \code{TRUE}.}
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
  \item{sub}{Character string for subelement of \code{qbObject} to
    get.}
  \item{dir}{Character string name of directory for load if \code{qbObject}
    does not exist.}
  \item{file}{Character string name of file for load if \code{qbObject}
    does not exist.}
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

  \code{qb.save} and \code{qb.load} are used to speed up examples and
  vignettes by saving and loading an external workspace containing a
  \code{cross} and a \code{qb} object. \code{qb.save} returns
  \code{TRUE} if the save was successful. \code{qb.load} returns
  \code{TRUE} if the objects exist already or they were loaded properly.
}
\examples{
\dontshow{
qb.load(cross, qbExample)


}
qb.get(qbExample, "output.dir")
summary(qb.cross(qbExample))
temp <- qb.get(qbExample, "iterdiag")
dim(temp)
names(temp)

\dontrun{
## The following should have no effect.
qbExample <- qb.reorder(qbExampleb}

## You can call the following rather than demo() to get a tour.
qb.demo()
}
}
\keyword{utilities}
