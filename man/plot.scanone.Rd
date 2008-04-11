\name{plot.qb.scanone}
\alias{plot.qb.scanone}
\alias{print.qb.scanone}
\alias{print.summary.qb.scanone}
\title{Plot or print qb.scanone object.}
\description{
Plot or print marginal diagnostics of effects from a qb.scanone object.
}
\usage{
\method{print}{summary.qb.scanone}(x, digits = 3, \dots)
\method{print}{qb.scanone}(x, digits = 3, \dots)
\method{plot}{qb.scanone}(x, chr, scan, ylim, scan.name, \dots)
}

\arguments{
  \item{x}{An object of class \code{qb.scanone}.}
  \item{digits}{Significant digits to round with print.}
  \item{chr}{Vector of chromosomes to plot. Must be integer.}
  \item{scan}{The model effects to include, the default is all
    those included in the \code{\link[qtl]{scanone}} object \code{x}.}
  \item{ylim}{Vector of length 2 with vertical limits.}
  \item{scan.name}{Name of \code{scan} for plot; default is
    \code{"effects"} or comma-separated catenation of \code{scan}.}
  \item{\dots}{Other values passed to the generic plot function.}
}

\value{
  Data frame with colors and line types used in plots.
}
\details{
This plot method uses \code{\link[qtl]{plot.scanone}} as the engine to plot
marginal posterior diagnostics created with \code{\link{qb.scanone}}. When
there are multiple effects in \code{x}, these may be organized into one
or several stacked plots using \code{scan}. The default for most
diagnostics except counts is \code{scan} = c("sum", "main",
"epis"). Counts and posterior diagnostics are typically plotted in two
stacked plots. Individual columns from the \code{x} object can be
plotted by specifying their names as a vector to option \code{scan}.

Arguments \code{col} and \code{lty} use keywords to names in \code{scan}
argument. Default \code{main} and \code{sub} arguments are provided.
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link{qb.scanone}}, \code{\link{summary.qb.scanone}},
  \code{\link[qtl]{plot.scanone}}}
\examples{
example(qb.scanone)
}
\keyword{models}
\keyword{hplot}
