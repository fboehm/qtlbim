#####################################################################
##
## $Id: qb.r,v 1.6 2005/11/28 byandell@wisc.edu
##
##     Copyright (C) 2005 Brian S. Yandell
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
qb.reorder <- function(qbObject)
{
  ## Need only do this once (however, it is idempotent).
  if(!is.null(qbObject$subset))
    return(qbObject)

  ## Sequence of iterations in iterdiag.
  iterdiag <- qb.get(qbObject, "iterdiag", sub = NULL)
  subs <- list(iterdiag = seq(nrow(iterdiag)))
  
  ## Need to reorder QTL in mainloci by niter, then chrom, then locus.
  mainloci <- qb.get(qbObject, "mainloci", NULL)
  if(!is.null(mainloci))
    subs$mainloci <- order(mainloci$niter, mainloci$chrom, mainloci$locus)

  ## Subset for GxE interaction: order by niter, covar, chrom, locus.
  gbye <- qb.get(qbObject, "gbye", NULL)
  if(!is.null(gbye))
    subs$gbye <- order(gbye$niter, gbye$covar, gbye$chrom, gbye$locus)
  
  ## Need to reorder QTL1, QTL2 in pairloci so that:
  ##      chrom1 < chrom2
  ##      locus1 < locus2 if chrom1 == chrom2
  p <- qb.get(qbObject, "pairloci", sub = NULL)
  if(is.null(p))
    subs$pairloci <- list(order = 0, flip = 0, left = 0, right = 0)
  else {
    s <- list(order = seq(nrow(p)))
    s$flip <- (p$chrom1 > p$chrom2) |
      ((p$chrom1 == p$chrom2) & (p$locus1 > p$locus2))
    ## order chromosome/locus sets
    if(any(s$flip)) {
      s$left <- c("chrom1","locus1")
      s$right <- c("chrom2","locus2")
      if(!is.null(p$ad)) {
        ## f2
        s$left <- c(s$left,"ad")
        s$right <- c(s$right,"da")
      }
    }
    else {
      s$left <- s$right <- 0
    }
    subs$pairloci <- s
  }
  ## Subset regions for chromosomes.
  map <- pull.map(qb.cross(qbObject))
  subs$region <- as.data.frame(lapply(map, range))
  subs$region <- data.frame(chr = seq(ncol(subs$region)),
                               start = unlist(subs$region[1,]),
                               end = unlist(subs$region[2,]))
  row.names(subs$region) <- names(map)

  ## Now attach subset to qb object.
  qbObject$subset <- subs
  qbObject
}
##############################################################################
qb.get <- function(qbObject, element, sub = qbObject$subset[[element]])
{

  ## Get mcmc output.
  renamedElement<-paste(element,".dat",sep="")
  x <- qbObject[[element]]
  if(!is.null(x))
    return(x)
  filename <- file.path(qbObject$output.dir, renamedElement) 
  if(!file.exists(filename))
    return(NULL)
  tmp <- scan(filename, n = 1, quiet = TRUE)
  if(length(tmp))
     x <- as.data.frame(read.table(filename))
  if(is.null(x)) 
    return(NULL)
  
  ## Specify variable names.
  is.bc <- (qbObject$cross == "bc")
  eff1 <- "add"
  eff2 <- "aa"
  var1 <- "varadd"
  var2 <- "varaa"
  if(!is.bc) {
      eff1 <- c(eff1,"dom")
      eff2 <- c(eff2,"ad","da","dd")
      var1 <- c(var1,"vardom")
      var2 <- c(var2,"varad","varda","vardd")
  }

  v <- var1
  if(qbObject$epistasis) v <- c(v,var2)
  if(qbObject$qtl_envi) {
     v <- c(v,"envadd")
     if(!is.bc) v <- c(v,"envdom")
  }
  if(qbObject$envi) v <- c(v,"varenv")
  v <- c(v,"var")

  xnames <- switch(element,
      covariates = {
        if(is.null(x))
          c("cov1")
        else
          paste("cov", seq(ncol(x)), sep = "")
      },
      gbye = c("niter","n.gbye","covar","chrom","locus",eff1,var1),
      mainloci = c("niter","nqtl","chrom","locus",eff1,var1),
      pairloci = c("niter","n.epis","chrom1","locus1","chrom2","locus2",eff2,var2),
      iterdiag = c("niter","nqtl","mean","envvar",v))

## assign the names to X
  if(is.null(x))
  {
      if(!is.null(xnames)) {
        x <- data.frame(matrix(NA, 0, length(xnames)))
        names(x) <- xnames
      }
   }
   else {
     if(element == "covariates")
        element <- "iterdiag"
     names(x) <- xnames
     if(!is.null(sub)) {
        if(element == "pairloci") {
          x[sub$flip, c(sub$left,sub$right)] =
            x[sub$flip, c(sub$right,sub$left)]
          x <- x[sub$order, ]
        }
      else
        x <- x[sub, ]
     }
   }
  grid <- pull.grid(qbObject, offset = TRUE, mask.region = FALSE)
  chrpos <- paste(grid$chr,
                  unlist(apply(as.matrix(table(grid$chr)), 1,
                               function(x) {
                                 seq(0, length = x)
                               })),
                  sep = ":")
  if(element == "mainloci" | element == "gbye") {
    x$locus <- grid$pos[match(paste(x$chrom, x$locus, sep = ":"), chrpos)]
  }
  if(element == "pairloci") {
    x$locus1 <- grid$pos[match(paste(x$chrom1, x$locus1, sep = ":"), chrpos)]
    x$locus2 <- grid$pos[match(paste(x$chrom2, x$locus2, sep = ":"), chrpos)]
  } 
  x
}

##############################################################################
qb.remove <- function(qbObject, verbose = TRUE)
{
  qbName <- deparse(substitute(qbObject))
  tmp <- qb.get(qbObject, "output.dir")
  if(dirname(tmp) != system.file("external", package = "qtlbim")) {
    if(verbose)
      warning(paste("Removing internal", qbName, "and external directory",
                    tmp),
              call. = FALSE, immediate. = TRUE)
    unlink(tmp, recursive = TRUE)
  }
  else {
    if(verbose)
      warning(paste("Removing internal", qbName),
              call. = FALSE, immediate. = TRUE)
  }
  remove(list = qbName, pos = 1)
}
##############################################################################
qb.exists <- function(qbObject)
{
  tmp <- qb.get(qbObject, "output.dir")
  if(!file.exists(tmp)) {
    cat("\n")
    stop(paste("Object contains no MCMC samples in", tmp), call. = FALSE)
  }
  invisible()
}
##############################################################################
qb.recover <- function(cross, traitName,
                        ## Find existing directory name if it exists.
                        output.dir = system(paste("ls -d ./", traitName, "_*",
                          sep = ""),
                          intern = TRUE),
                        ## Guess at other parameters.
                        n.thin = 40,
                        n.burnin = 0.01 * n.iter * n.thin,
                        algorithm = "M-H",
                        genoupdate = FALSE,
                        ## Arguments for qb.data and qb.model (except pheno.col).
                        ...)
{
  if(length(output.dir) == 0)
    stop("No trait MCMC directories exist to be recovered")
  if(length(output.dir) > 1) {
    print(output.dir)
    stop("Multiple trait MCMC directories exist: select one as output.dir")
  }

  step <- attr(cross$geno[[1]]$prob, "step")
  if(is.null(step))
    step <- 2
  
  qbObject = list(cross.name = deparse(substitute(cross)),
    cross = class(cross)[1],
    output.dir = output.dir,
    n.iter = 0,
    n.thin = n.thin,
    n.burnin = 0,
    algorithm = algorithm,
    genoupdate = genoupdate,
    step = step,
    verbose = TRUE)

  ## Now get the defaults (or supplied values via ...) for data and model.
  qbObject <- c(qbObject,
    qb.data(cross, pheno.col = find.pheno(cross, traitName), ...),
    qb.model(cross, ...))

  ## Get the actual number of iterations.
  ## Need to first have qbObject with some settings before qb.get call.
  qbObject$n.iter <- n.iter <- nrow(qb.get(qbObject, "iterdiag"))
  qbObject$n.burnin <- n.burnin

  ## Object is of class "qb".
  class(qbObject) <- "qb"

  ## Set up subset reordering.
  qbObject$subset <- NULL
  qb.reorder(qbObject)
}
##############################################################################
subset.qb <- function(x, nqtl = 1, pattern = NULL, exact = FALSE, chr,
                       region, offset = TRUE, restrict.pair = TRUE, ...)
{
  cross <- qb.cross(x)
  
  osub <- x$subset
  if(is.null(osub))
    x <- qb.reorder(x)
  if(missing(nqtl) & missing(pattern) & missing(exact) & missing(chr) &
     missing(region))
    return(x)
  nqt <- nqtl[1]
  if(!is.null(pattern)) {
    if(is.character(pattern))
      stop("pattern must be numeric, not character")
    nqtl <- max(nqtl, length(pattern))
  }
  
  ## Get elements (except covariate, which matches iterdiag)
  ## Note that previous subsetting is done here (see Incorporate below).
  iterdiag <- qb.get(x, "iterdiag")
  mainloci <- qb.get(x, "mainloci")
  pairloci <- qb.get(x, "pairloci")
  gbye <- qb.get(x, "gbye")

  ## And set up subset list on nqtl.
  ## These are TRUE/FALSE indicators except for region.
  sub <- list()

  ## Subset on number of QTL (led by iterdiag).
  if(exact) {
    ## Only exactly nqtl QTL (see also pattern below).
    sub$iterdiag <- iterdiag$nqtl == nqtl
    if(!sum(sub$iterdiag))
      stop(paste("empty object: no iterations with number of QTL =", nqtl))
  }
  else {
    ## At least nqtl QTL (see also pattern below).
    sub$iterdiag <- iterdiag$nqtl >= nqtl
    if(!sum(sub$iterdiag))
      stop(paste("empty object: no iterations with number of QTL >=", nqtl))
  }
  ## Now adjust other MCMC elements based on iterdiag.
  iters <- iterdiag[sub$iterdiag,"niter"]
  sub$mainloci <- !is.na(match(mainloci$niter, iters))
  sub$pairloci <- !is.na(match(pairloci$niter, iters))
  sub$gbye <- !is.na(match(gbye$niter, iters))

  ## Subset on pattern of chromosomes.
  if(!is.null(pattern)) {
    mypat <- table(pattern)
    if(exact) {
      ## Create function to only exactly this pattern retained in subset.
      mypat <- c(mypat, extra = 0)
      patfn <- function(x, mypat, yourpat) {
        tbl <- table(x)
        tmp <- match(names(tbl), names(mypat), nomatch = length(mypat))
        yourpat[tmp] <- tbl[tmp > 0]
        all(yourpat == mypat)
      }
    }
    else {
      ## Create function to retain all pattern sets that include this pattern.
      patfn <- function(x, mypat, yourpat) {
        tbl <- table(x)
        tmp <-  match(names(tbl), names(mypat), nomatch = 0)
        yourpat[tmp] <- tbl[tmp > 0]
        all(yourpat >= mypat)
      }
    }
    ## Now adjust all MCMC elements based on pattern in mainloci.
    blank <- rep(0, length(mypat))
    names(blank) <- names(mypat)
    iters <- unlist(tapply(mainloci$chrom, mainloci$niter, patfn,
                            mypat, blank, simplify = FALSE))
    sub$iterdiag <- sub$iterdiag & iters
    if(!sum(sub$iterdiag))
      stop(paste("empty object: no patterns like", pattern))
    iters <- iterdiag[sub$iterdiag, "niter"]
    sub$mainloci <- sub$mainloci & !is.na(match(mainloci$niter, iters))
    sub$gbye <- sub$gbye & !is.na(match(gbye$niter, iters))
    sub$pairloci <- sub$pairloci & !is.na(match(pairloci$niter, iters))
  }
  
  ## Subset of regions in chromosomes.
  sub$region <- x$subset$region
  if(!missing(region)) {
    region <- as.list(region)
    iters <- rep(TRUE, length(unique(mainloci$niter)))
    if(!offset)
      cross.map <- pull.map(cross)
    for(i in seq(length(region$chr))) {
      if(!offset) {
        region$start[i] <- region$start[i] + cross.map[[region$chr[i]]][1]
        region$end[i] <- region$end[i] + cross.map[[region$chr[i]]][1]
      }
      ## Score positive if a QTL in every region for an iteration.
      iters <- iters & tapply(mainloci$chrom == region$chr[i] &
                              (mainloci$locus >= region$start[i] &
                               mainloci$locus <= region$end[i]),
                              mainloci$niter, any)
      sub$region$start[region$chr[i]] <-
        max(sub$region$start[region$chr[i]], region$start[i])
      sub$region$end[region$chr[i]] <- 
        min(sub$region$end[region$chr[i]], region$end[i])
    }
    ## Note that some iterations now have 0 QTL; need careful handling.
    tmp <- match(names(iters), iterdiag$niter)
    sub$iterdiag[tmp] <- sub$iterdiag[tmp] & iters
    if(!sum(sub$iterdiag))
      stop(paste("empty object: no regions like", as.data.frame(region)))
    iters <- iterdiag[sub$iterdiag, "niter"]
    sub$mainloci <- sub$mainloci & !is.na(match(mainloci$niter, iters))
    sub$gbye <- sub$gbye & !is.na(match(gbye$niter, iters))
    sub$pairloci <- sub$pairloci & !is.na(match(pairloci$niter, iters))

    ## Drop linked QTL outside of region.
    if(!is.null(mainloci)) {
      tmp <- rep(FALSE, nrow(mainloci))
      for(i in seq(length(region$chr))) {
        tmp <- tmp | (mainloci$chrom == region$chr[i] &
                      (mainloci$locus < region$start[i] |
                       mainloci$locus > region$end[i]))
      }
      sub$mainloci <- sub$mainloci & !tmp
    }
    if(!is.null(gbye)) {
      tmp <- rep(FALSE, nrow(gbye))
      for(i in seq(length(region$chr))) {
        tmp <- tmp | (gbye$chrom == region$chr[i] &
                      (gbye$locus < region$start[i] |
                       gbye$locus > region$end[i]))
      }
      sub$gbye <- sub$gbye & !tmp
    }
    if(!is.null(pairloci)) {
      tmp <- rep(FALSE, nrow(pairloci))
      for(i in seq(length(region$chr))) {
        tmp <- tmp | ((pairloci$chrom1 == region$chr[i] &
                       (pairloci$locus1 < region$start[i] |
                        pairloci$locus1 > region$end[i])) |
                      (pairloci$chrom2 == region$chr[i] &
                       (pairloci$locus2 < region$start[i] |
                        pairloci$locus2 > region$end[i])))
      }
      sub$pairloci <- sub$pairloci & !tmp
    }
  }

  ## Subset of chromosomes.
  if(!missing(chr)) {
    ## First find out what form chr is in:
    n.chr <- nchr(cross)
    ## Logical: include chr if TRUE.
    if (is.logical(chr)) {
      if (length(chr) != n.chr) 
        stop(paste("If logical, chr argument must have length", 
                   n.chr))
      chr <- (1:n.chr)[chr]
    }
    ## Numeric chr has index into chromosome name vector.
    if(is.numeric(chr)) {
      ## Negative numbers: chr to exclude.
      if(all(chr < 1)) 
        chr <- (1:n.chr)[chr]
      ## Check range.
      if(any(chr < 1 | chr > n.chr)) 
        stop("Chromosome numbers out of range.")
    }
    else {
      ## Character chr has names of chromosomes.
      if(any(is.na(match(chr, names(x$geno))))) 
        stop("Not all chromosome names found.")
      chr <- match(chr, names(cross$geno))
    }
    chr <- sort(unique(chr))

    ## Subset all MCMC elements based on chr.
    kept <- !is.na(match(seq(length(cross$geno)), chr))
    sub$mainloci <- sub$mainloci & kept[mainloci$chrom]
    sub$pairloci <- if(restrict.pair)
      sub$pairloci & kept[pairloci$chrom1] & kept[pairloci$chrom2]
    else
      sub$pairloci & (kept[pairloci$chrom1] | kept[pairloci$chrom2])
    sub$gbye <- sub$gbye & kept[gbye$chrom]
    ## Arrange to drop other chr via pull.grid function.
    for(i in seq(n.chr)[-chr])
      sub$region[i,] <- c(i,0,-1)
  }

  ## Incorporate new subset on top of any previous subsetting.
  x$subset$iterdiag <- x$subset$iterdiag[sub$iterdiag]
  x$subset$mainloci <- x$subset$mainloci[sub$mainloci]
  x$subset$gbye <- x$subset$gbye[sub$gbye]
  x$subset$pairloci$order <- x$subset$pairloci$order[sub$pairloci]
  x$subset$region <- sub$region
  x
}
##############################################################################
qb.save <- function(cross, qbObject,
                    dir = ".",
                    file = paste(".RData.", Name, sep = "."))
{
  crossName <- deparse(substitute(cross))
  qbName <- deparse(substitute(qbObject))
  Name <- substring(qbName, 3)
  file <- file.path(dir, file)
  out <- exists(qbName) & exists(crossName)
  if(out) {
    .tryResults <- try(save(cross, qbObject, file))
    out <- !(class(.tryResults) == "try-error")
  }
  invisible(out)
}
##############################################################################
qb.load <- function(cross, qbObject,
                    dir = system.file("external", package = "qtlbim"),
                    file = paste(Name, "RData", sep = "."))
{
  crossName <- deparse(substitute(cross))
  qbName <- deparse(substitute(qbObject))
  Name <- substring(qbName, 3)
  file <- file.path(dir, file)
  out <- exists(qbName) & exists(crossName)
  if (!out) {
    .tryResults <- try(load(file = file, parent.frame()))
    out <- !(class(.tryResults) == "try-error")
    ## Check of cross and qbObject now exist.
    if(out & (!exists(qbName) | !exists(crossName)))
      out <- FALSE
    if(out) {
      qbObject$output.dir <-
        file.path(dir, basename(qb.get(qbObject, "output.dir")))
      assign(qbName, qbObject, envir = parent.frame())
      step <- qb.get(qbObject, "step")
      if(is.null(step))
        step <- 2
      cross <- qb.genoprob(cross, step = 2)
      assign(crossName, cross, envir = parent.frame())
    }
  }
  invisible(out)
}
