#####################################################################
##
## $Id: qb.R,v 1.16.2.9 2006/10/24 15:22:26 byandell Exp $
##
##     Copyright (C) 2002 Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
##############################################################################


##############################################################################
## summary.qb
##     Summary method for objects of class qb.
##
## arguments
##     object          An object of class qb (from run.bmapqtl),
##                     qb.BayesFactor (from qb.BayesFactor), or qb.qtl (from qb.qtl
##                     or qb.effects)
##
##     ...             Additional arguments to be passed to the generic
##                     method base::summary.
## returns
##     None


##############################################################################
print.qb <- function(x, ...) summary(x)
##############################################################################
summary.qb <- function(object, cutoff = 1, ...)
{
  cat("Bayesian model selection QTL mapping object",
      substitute(object), "on cross object",
      qb.get(object, "cross.name"), "\n")
  cat("had",
      qb.get(object, "n.iter"), "iterations recorded at each",
      qb.get(object, "n.thin"), "steps\n")
  cat("with",
      qb.get(object, "n.burnin"), "burn-in steps.\n")

  if(!is.legacy(object))
    cat("MCMC runs saved in qb object.\n")
  else
    cat("MCMC runs saved in temporary directory\n",
        qb.get(object, "output.dir"),
        "\n(use qb.remove to remove).\n")

  cross <- qb.cross(object, genoprob = FALSE)
  pheno.col <- qb.get(object, "pheno.col")
  cat(paste("Trait", ifelse(length(pheno.col) == 1, "", "s"), sep = ""),
      paste(names(cross$pheno)[pheno.col], collapse = ","),
      "(", paste(pheno.col, collapse = ","), ") treated as",
      qb.get(object, "trait"), ".\n")
  
  nfixcov <- qb.get(object, "nfixcov")
  nrancov <- qb.get(object, "nrancov")
  covar <- qb.get(object, "covar")
  if(nfixcov)
    cat(paste("Fixed covariate", "s"[nfixcov > 1], ":", sep = ""),
        paste(names(cross$pheno)[covar[seq(nfixcov)]], collapse = ", "),
        "\n")
  if(nrancov)
    cat(paste("Random covariate", "s"[nrancov > 1], ":", sep = ""),
        paste(names(cross$pheno)[covar[nfixcov+seq(nrancov)]],
              collapse = ", "),
        "\n")

  ## Partially sure on the following.
  cat("Trait was", ifelse(qb.get(object, "standardize"), "", "not"),
      "standardized.\n")
  cat(paste("Epistasis was",
            ifelse(qb.get(object, "epistasis"), "", " not"),
            " allowed.\n", sep = ""))
      
  cat("Prior number of QTL:",
      qb.get(object, "main.nqtl"), "main,",
      qb.get(object, "mean.nqtl"), "total, with",
      qb.get(object, "max.nqtl"), "maximum.\n")
  cat("Minimum distance between QTL:\n")
  tmp <- signif(qb.get(object, "interval"), 3)
  print(tmp)
  cat("Maximum number of QTL:\n")
  tmp2 <- qb.get(object, "chr.nqtl")
  names(tmp2) <- names(tmp)
  print(tmp2)
  
  cat(paste("QTL by environment",
            ifelse(qb.get(object, "qtl_envi"), "", " not"),
            " allowed.\n", sep = ""))
  cat("Interacting covariates:", qb.get(object, "intcov"))
##  cat("Contrast:", qb.get(object, "contrast"), "\n")

  qb.exists(object)
  
  ## Summaries of results.
  iterdiag <- qb.get(object, "iterdiag", ...)
  cat("\nDiagnostic summaries:\n")
  print(apply(iterdiag[ , -1], 2, summary))

  cat("\nPercentages for number of QTL detected:")
  print(round(100 * table(iterdiag$nqtl) / qb.niter(object)))

  if(qb.get(object, "epistasis")) {
    cat("\nPercentages for number of epistatic pairs detected:\n")
    print(qb.pair.nqtl(object, cutoff, ...))
    cat("\nPercentages for common epistatic pairs:")
    print(qb.pair.posterior(object, cutoff, 15, ...))
  }
  ## Need to add covariate stuff here.
  invisible()
}


##############################################################################
## qb.prior
##    This is a private method, not visible to package end-users.  Use this
## method to find the prior (distribution on the number of qtl) over a range
## of values.
##
## arguments
##     qbObject       An object of class qb.
##
##     range           A vector of nonnegative intesters.
##
## returns
##    Prior probabilities for each entry in the range argument arranged in
## the same order and with names corresonding to range entries.
##

qb.prior <- function(qbObject, range, mean = qb.get(qbObject, "mean.nqtl"))
{
  ## Extract prior distribution's name and mean.
  prior <- "poisson"

  ## Compute values of prior distribution over the range.
  pr <- switch(prior,
               geometric = (1 / mean) ^ range / (1 - (1 / mean)),
               poisson   = dpois(range, mean),
               uniform   = rep(1 / (1 + mean), length(range)))

  if(is.null(pr))
    stop("only geometric, poisson, and uniform priors recognized\n")

  ## Make sure that the values given in the range argument are used as
  ## the keys in the return value (pr).
  names(pr) <- range
  pr
}

##############################################################################
## plot.qb
##       This is a generic function for plotting qb objects.
## arguments
##       qbObject               An object of type qb.  
##       
##       
##

plot.qb <- function(x, ask = dev.interactive(), verbose = TRUE, ...)
{
  nqtl <- qb.nqtl(x, ...)
  if(max(nqtl) == 0)
    stop("no QTL found in MCMC runs")
  
  qb.exists(x)
  
  ## Now get a series of plots after prompts.

  cat("\nTime series of mcmc runs\n")
  tmp <- qb.coda(x, ...)
  if(verbose)
    print(summary(tmp, ...))
  plot(tmp, ...)
  
  cat("\nJittered plot of quantitative trait loci by chromosome...\n")
  tmp <- qb.loci(x, ...)
  if(verbose)
    print(summary(tmp, ...))
  plot(tmp, ...)
  
  cat("\nBayes Factor selection plots...\n")
  tmp <- qb.BayesFactor(x, ...)
  if(verbose)
    print(summary(tmp, ...))
  plot(tmp, ...)
  
  cat("\nHPD regions and best estimates...\n")
  tmp <- qb.hpdone(x, ...)
  attr(tmp$profile, "qb") <- deparse(substitute(x))
  attr(tmp$effects, "qb") <- deparse(substitute(x))
  
  if(verbose)
    print(summary(tmp, ...))
  plot(tmp, ...)
  
  cat("\nEpistatic effects...\n")
  tmp <- qb.epistasis(x, ...)
  if(verbose)
    print(summary(tmp, ...))
  plot(tmp, ...)

  cat("\nSummary diagnostics as histograms and boxplots by number of QTL\n")
  tmp <- qb.diag(x, ...)
  if(verbose)
    print(summary(tmp, ...))
  plot(tmp, ...)
  invisible()
}
##############################################################################
qb.smooth <- function(x, y)
{
  ux <- unique(x)
  if(length(ux) < 50) {
    ##    smo <- list(x = sort(ux),
    ##                y = rep(mean(y), length(ux)))
    ##    smo$sd <- rep(mad(y), length(smo$x))
    lmy <- lm(y ~ x)
    smo <- list(x = sort(ux),
                y = predict(lmy, data.frame(x = sort(ux))),
                sd = rep(sqrt(sum(resid(lmy) ^ 2) / lmy$df.resid),
                  length(ux)))
  }
  else {
    smo <- smooth.spline(x, y)
    smo$sd <- sqrt(pmax(0,
                         smooth.spline(x, (y - predict(smo, x)$y) ^ 2)$y))
  }
  smo  
}
##############################################################################
qb.loci <- function(qbObject, loci = c("main", "epistasis", "GxE"),
                    covar = get.covar, ...)
{
  qb.exists(qbObject)
  
  locis <- c("all", "main", "epistasis", "GxE")
  loci <- locis[pmatch(tolower(loci), tolower(locis), nomatch = 1)]
  loci <- unique(loci)

  elements <- c("mainloci", "mainloci", "pairloci", "gbye")
  names(elements) <- locis

  out.list <- list()
  for(element in loci) {
    out <- qb.get(qbObject, elements[element], ...)
    if(is.null(out))
      break
    if(0 == nrow(out))
      break
    switch(element,
           all = {
             out.list[[element]] <- out[, c("chrom", "locus")]
           },
           main = {
             tmp <- apply(as.matrix(out[, grep("^var", names(out))]), 1, sum) > 0
             if(sum(tmp) > 0)
               out.list[[element]] <- out[tmp, c("chrom", "locus")]
           },
           epistasis = {
             ## Could get fancy here or in plot.qb.loci
           ## by coloring most frequent pairs (using qb.pair.posterior).
             tmp <- names(out)
             find.chrom <- grep("chrom", tmp)
             find.locus <- grep("locus", tmp)
             out.list[[element]] <- data.frame(chrom = c(t(out[, find.chrom])),
                                               locus = c(t(out[, find.locus])))
           },
           GxE = {
             ## Could add covar option to subset.
             get.covar <- qb.get(qbObject,
                                 "nfixcov")[as.logical(qb.get(qbObject,
                                                              "intcov"))]
             tmp <- !is.na(match(out$covar, covar))
             if(sum(tmp) > 0) {
               out <- out[tmp, c("chrom", "locus"), drop = FALSE]
               out.list[[element]] <- out
             }
           })
  }
  class(out.list) <- c("qb.loci", "list")
  attr(out.list, "map") <- pull.map(qb.cross(qbObject, genoprob = FALSE))

  ## get mean spacing between grid points
  grid <- diff(pull.grid(qbObject)$pos)
  grid <- mean(grid[grid > 0])
  attr(out.list, "grid") <- grid
  
  attr(out.list, "n.iter") <- qb.get(qbObject, "n.iter")
  attr(out.list, "cex") <- qb.cex(qbObject)
  out.list
}
##############################################################################
plot.qb.loci <- function(x, loci = names(x), labels = FALSE, amount = .35,
                         cex = attr(x, "cex"),
                         col = c(all = "gray50", main = "blue", epistasis = "purple",
                           GxE = "darkred", marker = "green"),
                         ...)
{
  amount <- pmax(0, pmin(0.45, amount))
  amount <- array(amount, 2)

  map <- attr(x, "map")
  nmap <- names(map)

  if(is.null(names(col))) {
    col <- array(col, length(loci))
    names(col) <- loci
  }
  if(is.na(match("marker", names(col))))
    col <- c(col, marker = "green")

  ## Jittered plot of sampled loci.
  rlocus <- range(unlist(lapply(x, function(x) range(x$locus)))) +
    c(-1,1) * amount[2] * attr(x, "grid")
  uchrom <- sort(unique(unlist(lapply(x, function(x) unique(x$chrom)))))
  nchrom <- length(uchrom)
  ochrom <- rep(0, max(uchrom))
  ochrom[uchrom] <- seq(nchrom)
  plot(c(1,nchrom) + c(-.5,.5), range(0,rlocus),
       type = "n", xaxt = "n", xlab = "", ylab = "", ...)
  axis(1, seq(nchrom), nmap[uchrom])
  mtext("chromosome", 1, 2)
  mtext(paste(loci, col[loci], sep = "=", collapse = ", "), 1, 3)
  mtext("MCMC sampled loci", 2, 2)
  cxy <- par("cxy")[2] / 4

  ## Add lines for markers (and names if markers = TRUE).
  for(i in seq(nchrom)) {
    ii <- uchrom[i]
    tmp <- map[[nmap[ii]]]
    for(j in tmp)
      lines(i + c(-.475,.475), rep(j, 2), col = col["marker"],
            lwd = 2)
  }

  ## Jitter points vertically and horizontally.
  ## Side by side columns for loci at each chr.
  nloci <- length(loci)
  for(i in seq(nloci)) {
    points(jitter(ochrom[x[[loci[i]]]$chrom] +
                  (2 * i - nloci - 1) * amount[1] / nloci,
                  , amount[1] / nloci),
           jitter(x[[loci[i]]]$locus, , amount[2] * attr(x, "grid")),
           cex = cex,
           col = col[loci[i]])
  }
  ## Could add points from other loci elements here?
  ## But it may be too busy.

  ## Add marker names if markers = TRUE.
  if(labels) 
    for(i in seq(nchrom)) {
      ii <- uchrom[i]
      tmp <- map[[nmap[ii]]]
      text(rep(i-.5, length(tmp)), cxy + tmp, names(tmp),
           adj = 0, cex = 0.5, col = col["marker"])
    }
  invisible()
}
##############################################################################
summary.qb.loci <- function(object, digit = 1, ...)
{
  n.iter <- attr(object, "n.iter")
  nmap <- names(attr(object, "map"))
  lapply(object,
         function(object, digit) {
           tmp <- round(t(sapply(tapply(object$locus, nmap[object$chrom],
                                        function(x) {
                                          c(n.qtl = length(x) / n.iter,
                                            quantile(x, c(.25,.5,.75)))
                                        }),
                                 c)),
                        digit)
           tmp <- tmp[order(-tmp[, "n.qtl"]), ]
           tmp <- tmp[tmp[, "n.qtl"] > 0 | seq(nrow(tmp)) < 3, ]
         }, 
         digit)
}
##############################################################################
print.qb.loci <- function(x, ...) invisible(print(summary(x, ...)))
##############################################################################
qb.numqtl <- function(qbObject, ...)
{
  iterdiag <- qb.get(qbObject, "iterdiag", ...)
  n.iter <- qb.niter(qbObject)

  ## Posterior number of QTL.
  posterior <- table(iterdiag$nqtl)
  rnames <- names(posterior)
  posterior <- as.numeric(posterior) / n.iter
  
  ## Prior number of QTL.
  prior <- qb.prior(qbObject, as.numeric(rnames))

  ## Bayes factor ratios = posterior / prior.
  bf <- posterior / prior
  bf <- bf / min(bf, na.rm = TRUE)
  
  ## bfse = approximate Monte Carlo standard error for bf
  ## (actually binomial error)
  ## note that this is rescaled since bf[1] forced to be 1
  tmp <- data.frame(posterior = posterior, prior = prior, bf = bf,
       bfse = sqrt((1 - posterior) / (posterior * n.iter)) * bf)
  row.names(tmp) <- rnames
  tmp
}
##############################################################################
qb.pattern <- function(qbObject, cutoff = 1, nmax = 15, epistasis = TRUE, ...)
{
  iterdiag <- qb.get(qbObject, "iterdiag", ...)
  n.iter <- qb.niter(qbObject)

  mainloci <- qb.get(qbObject, "mainloci", ...)
  pairloci <- qb.get(qbObject, "pairloci", ...)
  pattern <- qb.makepattern(qbObject, epistasis, iterdiag = iterdiag,
                            mainloci = mainloci, pairloci = pairloci)
                           
  posterior <- rev(sort(table(pattern)))
  posterior <- posterior / sum(posterior)
  tmp <- posterior >= cutoff / 100
  if(sum(tmp))
    posterior <- posterior[tmp]
  else {
    cat("warning: pattern posterior cutoff", cutoff,
        "is too large and is ignored\n",
        "posterior range is", range(posterior), "\n")
  }
  if(length(posterior) > nmax)
    posterior <- posterior[1:nmax]
  ucount <- match(names(posterior), pattern)

  ## prior for pattern
  rng <- max(iterdiag$nqtl)
  pr <- qb.prior(qbObject, 0:rng)
  bf <- posterior
  map <- pull.map(qb.cross(qbObject, genoprob = FALSE))
  chrlen <- unlist(lapply(map, max))
  nchrom <- length(chrlen)
  chrlen <- chrlen / sum(chrlen)
  
  fact <- rep(1, rng)
  for(i in 2:(rng+1)) 
    fact[i] <- fact[i-1] * i

  ## New plan. Use only subset of mainloci matching ucount's.
  ## bundle table and for loop into one.

  ## Set up prior using null.
  prior <- rep(pr[1], length(posterior))
  names(prior) <- names(posterior)

  ## Find prior proportional to qb.prior and lengths of chromosomes.
  ## Use factorial adjustments for multiple linked QTL.
  tmpfn <- function(x, chrlen, fact) {
    ct <- c(table(x))
    prod(chrlen[names(ct)] ^ ct) * fact[sum(ct)] / prod(fact[ct])
  }

  ## Subset on non-null iterations corresponding to posterior.
  sub.post <- iterdiag$niter[ucount]
  is.depen <- qb.get(qbObject, "depen")
  if(epistasis & is.depen)
    tmp <- !is.na(match(pairloci$niter, sub.post))
  sub.post <- !is.na(match(mainloci$niter, sub.post))

  ## Adjust for epistatic effects.
  ## For now only considering hierarchical model case.
  ## That is epistasis only if main effects.
  ## In order to handle epistasis with 1 or no main effect
  ## we would need to compute probabilities for each iteration.
  ## This will be a lot more work!
  if(epistasis) {
    if(is.depen) {
      ## Epistatic priors depend on number of main loci.
      ## Not correct yet when c1 or c0 > 0.
      ## Find number of possible epistatic pairs of main loci.
      tbl.main <- c(table(mainloci[sub.post, "niter"]))
      tmp2 <- tbl.main * (tbl.main - 1) / 2
      
      tmp <- c(table(pairloci[tmp,"niter"]))
      prop <- qb.get(qbObject, "prop")
      if(sum(prop[-1]) > 0)
        warning("Bayes factor computations assume all QTL are main (may not be correct)")

      ## Epistasis with two main QTL.
      tmp2[names(tmp)] <- (prop[1] ^ tmp) *
        ((1 - prop[1]) ^ (tmp2[names(tmp)] - tmp))
      tmp2[is.na(match(names(tmp2), names(tmp)))] <-
        (1 - prop[1]) ^ tmp2[is.na(match(names(tmp2), names(tmp)))]
      tbl.main <- pr[1 + tbl.main] * tmp2
    }
    else {
      ## Independent prior for epistasis.
      ## Product of prior for main QTL * prior for epis-only QTL.

      ## Find mainloci entries matching patterns with high posterior.
      ## Kludge to catch iterations with 0 QTL.
      tmp <- rep(0, nrow(iterdiag))
      names(tmp) <- iterdiag$niter
      tmp2 <- c(table(mainloci[, "niter"]))
      tmp[names(tmp2)] <- tmp2
      tmp2 <- rep(!is.na(match(pattern, names(posterior))), tmp)

      ## Sum up main loci variance components and tally those not zero.
      tmp <- mainloci[tmp2, "varadd"]
      is.bc <- (qb.cross.class(qbObject) == "bc")
      if(!is.bc)
        tmp <- tmp + mainloci[tmp2, "vardom"]
      tmp <- tapply(tmp, mainloci[tmp2, "niter"], function(x) sum(x > 0))
      ## Table number of total QTL per iteration.
      tbl.main <- c(table(mainloci[tmp2, "niter"]))

      ## Product of priors for main and epis-only QTL.
      main.nqtl <- qb.get(qbObject, "main.nqtl")
      mean.nqtl <- qb.get(qbObject, "mean.nqtl")
      ## Average product of poisson probabilities by pattern.
      pr.main <- qb.prior(qbObject, seq(0, max(tmp)), main.nqtl)
      pr <- qb.prior(qbObject, seq(0, max(tbl.main - tmp)), mean.nqtl - main.nqtl)
      tbl.main <- tapply(pr.main[1 + tmp] * pr[1 + tbl.main - tmp],
                         pattern[names(tbl.main)], mean)
      ## Rearrange in order to match calculations below.
      tbl.main <- tbl.main[pattern[as.character(sort(iterdiag$niter[ucount]))]]
      names(tbl.main) <- sort(iterdiag$niter[ucount])
      tbl.main <- tbl.main[!is.na(tbl.main)]
    }
  }
  else
    tbl.main <- pr[1 + c(table(mainloci[sub.post, "niter"]))]

  geno.names <- names(map)
  tmp <- c(tapply(ordered(geno.names[mainloci[sub.post, "chrom"]], geno.names),
                  mainloci[sub.post, "niter"],
                  tmpfn, chrlen, fact))
  tmp <- tmp * tbl.main
       
  prior[pattern[names(tmp)]] <- tmp

  ## Bayes factor ratio = rescaled version of posterior / prior.
  bf <- posterior / prior
  ## rescale bf so smallest value is 1 (avoid NA, 0/0)
  minbf <- bf[!is.na(bf)]
  if(length(minbf))
    bf <- bf / min(minbf)

  ## bfse = approximate Monte Carlo standard error for bf
  ## (actually binomial error)
  ## note that this is rescaled since bf[1] forced to be 1
  tmp <- names(posterior)
  nqtl <- sapply(strsplit(tmp, ","),length)
  names(nqtl) <- tmp
  data.frame(nqtl = nqtl,
             posterior = posterior, prior = prior, bf = bf,
             bfse = sqrt((1 - posterior) / (posterior * n.iter)) * bf)
}
##############################################################################
qb.bf <- function(...) qb.BayesFactor(...)
##############################################################################
qb.BayesFactor <- function(qbObject,
                           items = c("nqtl","pattern","chrom","pairs"),
                           cutoff.pattern = ifelse(epistasis, 0.25, 0.5),
                           cutoff.pairs = 1, nmax = 15,
                           epistasis = TRUE, ...)
{
  qb.exists(qbObject)
  
  assess <- list()

  if(any(pmatch(tolower(items), "nqtl", nomatch = 0)))
    assess$nqtl <- qb.numqtl(qbObject, ...)

  if(any(pmatch(tolower(items), "pattern", nomatch = 0)))
    assess$pattern <- qb.pattern(qbObject, cutoff.pattern, nmax, epistasis, ...)

  if(any(pmatch(tolower(items), "chrom", nomatch = 0)))
    assess$chrom <- qb.chrom(qbObject, ...)

  if(any(pmatch(tolower(items), "pairs", nomatch = 0)))
    assess$pairs <- qb.pairs(qbObject, cutoff.pairs, nmax, ...)

  class(assess) <- c("qb.BayesFactor", "list")
  assess
}
##############################################################################
summary.qb.BayesFactor <- function(object, sort = TRUE, digits = 3, ...)
{
  if(sort)
    object <- lapply(object, function(x) x[order(-x$bf), ])
  lapply(object, signif, digits)
}
##############################################################################
print.qb.BayesFactor <- function(x, ...) print(summary(x, ...))
##############################################################################
plot.qb.BayesFactor <- function(x, ...)
{
  tmpar <- par(mfrow = c(length(x), 2), mar = c(3.1,3.1,2.1,0.1))
  on.exit(par(tmpar))

  if(!is.null(x$nqtl))
    plot.qb.pattern(x$nqtl, row.names(x$nqtl),
                     c("number of QTL","QTL posterior","QTL posterior"),
                     NULL, ...)

  if(!is.null(x$pattern))
    plot.qb.pattern(x$pattern, ,
                    c("pattern of QTL","pattern posterior","pattern posterior"),
                    ...)

  if(!is.null(x$chrom))
    plot.qb.pattern(x$chrom, row.names(x$chrom),
                     c("chromosome","chrom posterior","chrom posterior"),
                     NULL,...)

  if(!is.null(x$pairs))
    plot.qb.pattern(x$pairs, row.names(x$pairs),
                     c("pairs", "pairs posterior", "pairs posterior"),
                     NULL, ...)

  invisible()
}
##############################################################################
### marginal histograms
##############################################################################
qb.diag <- function(qbObject, items= c("mean","envvar","var","herit"), ...)
{
  qb.exists(qbObject)
  
  iterdiag <- qb.get(qbObject, "iterdiag", ...)
  nhist <- length(items)
  diag <- matrix(NA, qb.niter(qbObject), nhist + 1)
  dimnames(diag) <- list(NULL, c(items, "nqtl"))

  for(i in seq(nhist)) {
    ## marginal histogram
    if(items[i] == "herit")
      diag[,i] <- iterdiag$var / (iterdiag$envvar + iterdiag$var)
    else
      diag[,i] <- iterdiag[[items[i]]]
  }
  diag[, 1 + nhist] <- iterdiag$nqtl
  class(diag) <- c("qb.diag", "matrix")
  diag
}
##############################################################################
summary.qb.diag <- function(object, digits = 5, ...)
{
  nc <- ncol(object)
  apply(object[, -nc], 2,
        function(x, y) {
          signif(c(median = median(x, na.rm = TRUE),
                  unlist(tapply(x, y, median, na.rm = TRUE))), 5)
          },
        object[, nc])
}
##############################################################################
print.qb.diag <- function(x, ...) print(summary(x, ...))
##############################################################################
plot.qb.diag <- function(x, ...)
{
  nx <- ncol(x)
  nhist <- nx - 1

  tmpar <- par(mfrow = c(nhist, 2), mar=c(2.1,2.1,0.1,0.1), oma = c(0,0,3,0))
  on.exit(par(tmpar))

  for(i in dimnames(x)[[2]][-nx]) {
    diag <- x[, i]
    diag <- diag[!is.na(diag)]
    
    plot(density(diag), main = "", xlab = "", ylab = "", yaxt = "n", ...)
    mtext(i, 2, 0.5)
    ## hand-made boxplot on its side
    b <- boxplot(diag, plot = FALSE)
    tmp <- par("usr")[4] / 8
    up <- tmp * 1.5
    polygon(b$stats[c(2,4,4,2)], up+c(0,0,tmp,tmp))
    lines(rep(b$stats[3],2), up+c(0,tmp))
    lines(b$stats[1:2], up+rep(tmp/2,2), lty = 2)
    lines(b$stats[4:5], up+rep(tmp/2,2), lty = 2)
    lines(rep(b$stats[1],2), up+tmp*c(1,3)/4)
    lines(rep(b$stats[5],2), up+tmp*c(1,3)/4)
    
    ## conditional boxplots
    tmp <- split(diag, x[!is.na(x[, i]), nx])
    boxplot(tmp)
  }
}
##############################################################################
qb.demo <- function()
{
  demo.list <- c("mcmc","plot","scan","coda","covar","hyper","sim")
  repeat {
    current.list <- demo.list
    selection <- readline(paste("Next demo (", paste(current.list, collapse = ","), "): ",
                   sep = ""))
    if(!nchar(selection)) break

    stops <- c("exit","quit")
    helps <- c("help","?")

    ## Selection is not empty, so proceed to find unique entry.
    while(nchar(selection)) {
      selection <- grep(paste("^", tolower(selection), ".*", sep = ""),
                        c(current.list, stops), value = TRUE)
      if(length(selection) <= 1) break
      current.list <- selection
      selection <- readline(paste("Next demo (", paste(current.list, collapse = ","), "):",
                      sep = ""))
    }
    if(length(selection)) {
      if(match(selection, stops, nomatch = 0)) break
      if(match(selection, helps, nomatch = 0))
        demo(package = "qtlbim")
      else {
        if(nchar(selection))
          demo(paste("qb", selection, "tour", sep = "."), "qtlbim",
               character.only = TRUE)
      }
    }
  }
  invisible()
}
