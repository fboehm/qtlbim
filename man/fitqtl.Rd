\name{step.fitqtl}
\alias{step.fitqtl}
\alias{anova.step.fitqtl}
\title{Stepwise backward elimination and anova comparison.}
\description{
These functions mimic \code{step} and \code{anova} but have reduced
functionality. They are not truly methods, but can help study
qtl model fits.
}
\usage{
step.fitqtl(cross, qtl, pheno.col = 1, arch, cutoff = 0.05,
  trace = 1, steps = 100)
\method{anova}{step.fitqtl}(object, object2, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{Object of class \code{cross}.}
  \item{qtl}{Object of class \code{qtl}, as output of
    \code{\link[qtl]{makeqtl}}.}
  \item{pheno.col}{Column of phenotype (numeric).}
  \item{arch}{Object of class \code{qb.arch} from \code{\link{qb.arch}}.}
  \item{cutoff}{Significance cutoff for dropping terms.}
  \item{trace}{If positive, information is printed during the
    run. Values 1, 2, 3 give gradually more detailed information.}
  \item{steps}{Maximum number of steps to be considered.}
  \item{object}{Object of class \code{step.fitqtl} from \code{step.fitqtl}.}
  \item{object2}{Object of class \code{step.fitqtl} from \code{step.fitqtl}.}
  \item{\dots}{Currently not used.}
}
\details{
\code{step.fitqtl} is analogous to \code{\link[stats]{step}} applied to
analysis with \code{\link[qtl]{fitqtl}}. \code{anova.step.fitqtl} is an
  S3 method for \code{\link[stats]{anova}}. \code{anova.step.fitqtl}
  with one argument calls \code{\link[qtl]{summary.fitqtl}}; with two
  arguments it attempts to conduct a general F comparison of anova fits.
}
\value{
  \code{step.fitqtl} returns an object of class \code{step.fitqtl} with
  \item{fit}{Object of class \code{\link[qtl]{fitqtl}}.}
  \item{arch}{Object of class \code{\link{qb.arch}}.}
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{qb.arch}}, \code{\link[qtl]{fitqtl}},
  \code{\link[qtl]{summary.fitqtl}},
  \code{\link[qtl]{makeqtl}}}
\examples{
\dontshow{
if(!exists("cross.arch") | !exists("qbExample"))
  example(qb.arch)
  cross <- qb.cross(qbExample)
}
cross <- sim.geno(cross, n.draws = 8, step = 2, error = 0.01)
qtl <- makeqtl(cross, chr = c(1,1,2,3), pos = c(15,45,12,15))
cross.step <- step.fitqtl(cross, qtl, pheno.col = 3, arch = cross.arch)
anova(cross.step)
cross.step$arch
}
\keyword{models}
