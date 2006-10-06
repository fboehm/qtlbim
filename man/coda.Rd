\name{qb.coda}
\alias{qb.coda}
\title{ Coerce to an MCMC object for use with the coda package.}
\description{
This function creates an object of class \code{\link[coda]{mcmc}} from
an object of class \code{qb} produced by the package \code{R/qb}.
}
\usage{
qb.coda( object, element,
  variables = c("nqtl","mean","envvar","var"))                            
}

\arguments{
  \item{object}{ An object of class \code{qb} returned by calling the function
                \code{\link{qb.mcmc}}.
                 }
  \item{element}{ A character string which has to one of \code{"iterdiag"}, \code{"mainloci"}, 
                \code{"pairloci"}, \code{"covariates"}, or
		\code{"gbye"}; default is \code{"iterdiag"}. }
  \item{variables}{ A vector of integers specifying the columns or
    column names of \code{element} to be considered. Details about the
    columns can be found in \code{\link{qb.mcmc}}}
}
\details{
This package requires the package \code{coda}.
}
\value{ An object of class \code{\link[coda]{mcmc}}. This object could be
used to analyze the MCMC output using R/coda.}
}


\author{Dr. Nengjun Yi, et al.,  nyi@ms.ssg.uab.edu}
 
 \references{
  \url{http://www.ssg.uab.edu/qtlbim}
  }

\seealso{ \code{\link[coda]{mcmc}}, \code{\link{qb.mcmc}},
  \code{\link{plot.qb}}}
\examples{     
\dontshow{
qb.load(cross, qbExample)
}
## Default plots for iteration diagnostics "iterdiag".
temp <- qb.coda(qbExample)
plot(temp)

## Summaries for some "mainloci" elements.
temp <- qb.coda(qbExample, "mainloci", c("chrom","locus","add","dom"))
plot(temp)
}
\keyword{manip}
