\name{qb.sweave}
\alias{qb.sweave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run sweave to automate QTL search with MCMC samples.}
\description{
This routine runs a separate Sweave file (*.Rnw) of commands, making
substitutions for the user-supplied data and thresholds. It can be used
to automate the search for genetic architecture.
}
\usage{
qb.sweave(cross, pheno.col = 1, n.iter = 3000, n.draws = 64,
  scan.type = "2logBF", hpd.level = 0.5,
  upper.threshold, SweaveFile, SweaveExtra, PDFDir, remove.qb = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{Object of class \code{cross}.}
  \item{pheno.col}{Phenotype column in object \code{cross}.}
  \item{n.iter}{Number of MCMC iterations to be stored.}
  \item{n.draws}{Number of MC draws to use for \code{fitqtl}.}
  \item{scan.type}{Type of 1-D and 2-D scan to perform; see
    \code{\link{qb.scanone}}.} 
  \item{hpd.level}{Highest probability density level for scan; see
    \code{\link{qb.hpdone}}.} 
  \item{upper.threshold}{Threshold for upper triangle (epistasis) in 2-D
    scan; see \code{\link{qb.scantwo}}.}
  \item{SweaveFile}{Name of \code{\link[utils]{Sweave}} file (default is
    \code{system.file("doc", "hyperslide.Rnw", package = "qtlbim")}).}
  \item{SweaveExtra}{Name of user-supplied extra
    \code{\link[utils]{Sweave}} file (default is \code{NULL}).} 
  \item{PDFDir}{Name of directory to store PDF files (default is
    \code{phenoPDF}, where \code{pheno} is the name associated with
    phenotype \code{pheno.col}).}
  \item{remove.qb}{Remove constucted objects if \code{TRUE}.}
}
\details{
This is a simple shell around the \code{\link[utils]{Sweave}} routine to
create customized documents with embedded QTL analysis. The default file
\code{system.file("doc", "hyperslide.Rnw", package = "qtlbim")} creates a
"beamer" style PDF slide show. An alternative file
\code{system.file("doc", "hyperpaper.Rnw", package = "qtlbim")} creates a
preprint document. Both require post-processing with \code{pdflatex}.

A user-defined section can be added to the automated documents, using
the \code{SweaveExtra} option. We have provided
\code{system.file("external", "hyperslideextra.Rnw", package = "qtlbim")}
for the slide version and
\code{system.file("external", "hyperpaperextra.Rnw", package = "qtlbim")}
for the preprint version.
}
\references{ \url{http://www.soph.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link[utils]{Sweave}}}
\examples{
\dontrun{
data(hyper)

## Create default slide show LaTeX source without extra section.
qb.sweave(hyper)

## Turn LaTeX into PDF. Run twice to get outline correct.
## Need pdflatex on your system.
system("pdflatex hyperslide")
system("pdflatex hyperslide")

## Create document form, with extra section.
qb.sweave(hyper,
  SweaveFile = system.file("doc", "hyperpaper.Rnw", package = "qtlbim"))
  SweaveExtra = system.file("external", "hyperpaperextra.Rnw",
    package = "qtlbim"))
system("pdflatex hyperpaper")
system("pdflatex hyperpaper")
}
}
\keyword{hplot}
