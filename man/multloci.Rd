\name{qb.multloci}
\alias{qb.multloci}
\alias{plot.qb.multloci}
\alias{summary.qb.multloci}
\title{Summaries of multiple loci on a chromosome.}
\description{
Summaries and up to four plots showing loci found in MCMC
samples for a chromosome.
}
\usage{
qb.multloci(qbObject, chr = 1, cutoff = 25, nqtl, \dots)
\method{plot}{qb.multloci}(x, amount = 0.5, cex,
  split = TRUE, contour = TRUE, weight = TRUE, merge = TRUE, \dots)
\method{summary}{qb.multloci}(object, merge = TRUE, \dots)
}
\arguments{
  \item{qbObject}{Object of class \code{qb}.}
  \item{x,object}{Object of class \code{qb.multloci}.}
  \item{chr}{Identifier for one chromosome.}
  \item{cutoff}{Smallest posterior probability for \code{nqtl} (ignored
    if \code{nqtl} provided).}
  \item{nqtl}{Number of QTL on chromosome (inferred by default).}
  \item{amount}{Amount to \code{\link[base]{jitter}} points.}
  \item{cex}{Character expansion of plot symbols.}
  \item{split}{Split plots into four panels on one page if
    \code{TRUE}. Otherwise plot each panel separately. The \code{split}
    may be a numeric vector with values 1:4 signifying which panels to
    show. See details.}
  \item{contour}{Contour plot overlaid on pairs if \code{TRUE}.}
  \item{weight}{Inversely weight loci in density plot by number of QTL
    if \code{TRUE}.}
  \item{merge}{Merge across number of QTL if \code{TRUE}. Otherwise,
    show separate summary or plot by number of QTL. See details.}
  \item{\dots}{Parameters to \code{\link{qb.mainmodes}} or to methods.}
}
\details{
  Find multiple loci in MCMC samples for chromosome \code{chr}. The
  number of QTL, \code{nqtl} is inferred from the histogram as the
  largest number of QTL above the percent \code{cutoff}.

  The generic \code{plot} command produces the following plots: (1)
  density plot of main QTL grouped by QTL; (2) histogram of number of
  QTL; (3) density plot of epistatic pairs; (4) scatter plot of pairs of
  QTL. The density plots are divided into \code{nqtl} groups. The
  scatter plot shows pairs of main loci below diagonal and epistatic
  pairs above using codes corresponding to the number of QTL per sample;
  note   that 3 QTL have 3 pairs, 4 QTL have 6, etc., and that solitary
  QTL are   displayed along the diagonal.

  \code{split} and \code{merge} control the manner of plotting. Setting
  \code{merge} to \code{FALSE} yields only density plots for main loci
  conditioned on the number of QTL per sample. Setting \code{split} to
  \code{FALSE} or to numbers between 1 and 4 yields plots on separate
  pages.
}
\references{ \url{http://www.qtlbim.org} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.qb}}, \code{\link{qb.scantwo}}, \code{\link{qb.mainmodes}}}

\examples{
data(qbExample)

temp <- qb.multloci(qbExample, 1)
summary(temp)
plot(temp)
plot(temp, merge = FALSE)
summary(temp, merge = FALSE)
}
\keyword{models}
