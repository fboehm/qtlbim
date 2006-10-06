qb.gbye.posterior <- function(qbObject, covar = get.covar, cutoff = 1, nmax = 5,
                              gbye = qb.get(qbObject, "gbye"))
{
  if(is.null(gbye))
    return(invisible(NULL))

  cross <- qb.cross(qbObject)
  geno.names <- names(cross$geno)
  get.covar <- qb.get(qbObject,
                      "nfixcov")[as.logical(qb.get(qbObject, "intcov"))]
  covar.names <- names(cross$pheno)[covar]

  percent <- table(covar.names[gbye$covar], geno.names[gbye$chrom])
  np <- dimnames(percent)
  np <- outer(np[[2]], np[[1]], paste, sep = ".")
  percent <- c(percent)
  o <- order(-percent)
  percent <- 100 * percent[o] / nrow(qb.get(qbObject, "iterdiag"))
  names(percent) <- np[o]
  percent <- percent[ percent >  cutoff ]
  if(length(percent) > nmax)
    percent <- percent[seq(nmax)]
  round(percent)
}
qb.intcov <- function(qbObject, covar = get.covar, effects = c("add","dom"),
                         cutoff = 1, nmax = 5, cov.chr = names(post), ...)
{
  ## Use qb.gbye.posterior to find interesting pairs.
  ## Set up object as in qb.epistasis.
  ## Set up plot or use plot.qb.epistasis.

  gbye <- qb.get(qbObject, "gbye")
  if(is.null(gbye)) {
    cat("no GxE\n")
    return(invisible(NULL))
  }

  ## Identify covariates and chromosomes with interacting QTL.
  geno.names <- names(qb.cross(qbObject)$geno)
  get.covar <- qb.get(qbObject,
                      "nfixcov")[as.logical(qb.get(qbObject, "intcov"))]
  if(!length(covar))
    return(invisible(NULL))
  if(all(covar <= 0))
    return(invisible(NULL))

  covar.names <- names(cross$pheno)[covar]

  inter <- interaction(geno.names[gbye[, "chrom"]],
                       covar.names[gbye[, "covar"]])
  post <- qb.gbye.posterior(qbObject, covar, cutoff, nmax, gbye = gbye)
  if(!is.character(cov.chr))
    stop("cov.chr must be character")
  cov.chr <- cov.chr[match(names(post), cov.chr, nomatch = 0)]
  post <- post[cov.chr]
  inter.cov.chr <- !is.na(match(inter, cov.chr))

  ## Subset on desired Cockerham effects.
  effects <- effects[match(names(gbye), effects, nomatch = 0)]
  gbye <- as.data.frame(as.matrix(gbye[inter.cov.chr, effects]))
  names(gbye) <- effects

  ## Kludge inter to add posterior to name of cov.chr.
  inter <- ordered(as.character(inter[inter.cov.chr]), cov.chr)
  cov.chr.pct <- paste(cov.chr, "\n", post, "%", sep = "")

  gbye$inter <- ordered(cov.chr.pct[unclass(inter)], cov.chr.pct)

  class(gbye) <- c("qb.epistasis", "data.frame")
  attr(gbye, "post") <- post
  gbye
}

