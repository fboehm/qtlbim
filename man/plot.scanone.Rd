\name{plot.qb.scanone}
\alias{plot.qb.scanone}
\alias{print.qb.scanone}
\title{Plot or print qb.scanone object.}
\description{
Plot or print marginal diagnostics of effects from a qb.scanone object.
}
\usage{
\method{print}{qb.scanone}(x, \dots)
\method{plot}{qb.scanone}(x, chr, smooth = 3, scan, ylim,
  scan.name, col, main, verbose = FALSE, \dots)
}

\arguments{
  \item{x}{An object of class \code{qb.scanone}.}
  \item{chr}{Vector of chromosomes to plot. Must be integer.}
  \item{smooth}{Perform smoothing if > 0 using weighted average over
    \code{smooth} adjacent points.}
  \item{scan}{The model effects to include, the default is all
    those included in the \code{\link[qtl]{scanone}} object \code{x}.}
  \item{ylim}{Limits for vertical (y) axis; default uses data limits.}
  \item{scan.name}{Name used in automatically generated main title.}
  \item{col}{Named vector of colors for plot. Names of colors correspond
    to effects to be plotted. Unnamed colors will be made \code{"black"}.}
  \item{main}{Main title for the plot}
  \item{verbose}{Give verbose feedback if \code{TRUE}.}
  \item{\dots}{Other values passed to the generic plot function.}
}

}
\value{Colors used in plots as character vector.
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
}
\references{ \url{http://www.ssg.uab.edu/qtlbim} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link{qb.scanone}}, \code{\link{summary.qb.scanone}},
  \code{\link[qtl]{plot.scanone}}}
\examples{
example(qb.scanone)
}
\keyword{models}
\keyword{hplot}
