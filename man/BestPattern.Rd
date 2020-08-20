\name{qb.BestPattern}
\alias{qb.best}
\alias{qb.BestPattern}
\alias{plot.qb.BestPattern}
\alias{summary.qb.BestPattern}
\title{Proximity of common genetic architecture patterns.}
\description{
Multidimensional scaling and hierarchical clustering of most common
patterns of genetic architecture.
}
\usage{
qb.best(\dots)
qb.BestPattern(qbObject, epistasis = TRUE,
  category = c("pattern", "nqtl"), cutoff, score.type =
  c("sq.atten","attenuation","variance","recombination","distance"),
  include = c("nested","all","exact"),
  center = c("median","mean"), level = 5, \dots)
\method{plot}{qb.BestPattern}(x, type = c("mds", "hclust"),
  main, xlab, method = "complete", cluster = 3, cexmax = 5,
  colmax = 75, cex, col,
  symbol = c("pattern","nqtl","cluster","c@n","c@p","n@p","c@n@p"), \dots)
\method{summary}{qb.BestPattern}(object, method = "complete", 
  cluster = 3, n.best = 1, \dots)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{x,object}{Object of class \code{qb.BestPattern}.}
  \item{epistasis}{Include epistasis in patterns if \code{TRUE}.}
  \item{category}{Distances indexed by \code{nqtl} or \code{pattern}.}
 \item{cutoff}{Percent cutoff for pattern inclusion in model
    selection. Default is 0.25 (0.5) if \code{epistasis} is \code{TRUE}
    (\code{FALSE}).}
  \item{score.type}{Type of score to use as distance. See
    \code{\link{qb.close}}.}
  \item{type}{Plot dendrogram for \code{hclust} or 2-D multidimensional
    scaling projection for \code{mds}.}
  \item{main}{Main plot title as character string.}
  \item{xlab}{Character string for horizontal (x) axis.}
  \item{method}{Method for hierarchical clustering.}
  \item{cluster}{Number of clusters desired.}
  \item{n.best}{Number of better models to display.}
  \item{cexmax}{Maximum font size (minimum is set to 1); patterns are
    displayed in \code{mds} plot proportional to their posterior
    probability.}
  \item{colmax}{Maximum number of colors.}
  \item{cex}{Manual override of font size for \code{mds} plot; should be
    length 1 or the number of patterns exceeding \code{cutoff}.}
  \item{col}{Colors for plotting.}
  \item{symbol}{Plot symbol for \code{mds} plot. Shorthand using at sign
    \code{@} signifies catenation of two or more symbols into one.}
  \item{include}{Action for model averaging of chromosome-specific locus
    and explained variance:
    use \code{all} MCMC samples that match the chromosome;
    use only MCMC samples for patterns that have the target pattern
    \code{nested} within them;
    or use only MCMC samples with the \code{exact} same target pattern.}
  \item{center}{Method of estimating the center for \code{locus} and
    explained \code{variance}.}
  \item{level}{Confidence level as percent between 0 and 100 for loci
    and variance contributions.}
  \item{\dots}{Parameters to methods.}
}
\details{
This uses the closeness measure from \code{\link{qb.close}} to
compute a similarity matrix among patterns whose posterior probabilities
exceed \code{cutoff}. Distance = 1 - similarity is used for hierarchical
clustering or multidimensional scaling.

The \code{best} pattern is chosen as the one with highest posterior
mean; all other patterns are compared to that pattern in terms of the
\code{score.type}. This best pattern is a natural \code{target} for \code{\link{qb.close}}.
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{qb.close}}}

\examples{
data(qbExample)

temp <- qb.BestPattern(qbExample)
summary(temp, n.best = 3)
plot(temp, type = "hclust")
plot(temp)
plot(temp, symbol = "c@n")

best <- summary(temp)$best
temp <- qb.close(qbExample, best)
summary(temp)
plot(temp)
}
\keyword{models}
