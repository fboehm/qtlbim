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
qb.reorder <- function(qbObject, warn = FALSE,
                          pheno.col = qb.get(qbObject, "pheno.col", warn = warn))
{
  ## Need only do this once per phenotype (however, it is idempotent).
  cross <- qb.cross(qbObject, genoprob = FALSE, warn = warn)
  pheno.names <- names(qb.find.pheno(qbObject, cross, pheno.col, warn = warn))
  
  region <- qb.get(qbObject, "region", warn = warn)
  if(is.null(region)) {
    ## Subset regions for chromosomes.
    map <- pull.map(cross)
    region <- as.data.frame(lapply(map, range))
    region <- data.frame(chr = seq(ncol(region)),
                         start = unlist(region[1,]),
                         end = unlist(region[2,]))
    row.names(region) <- names(map)
  }
  
  if(is.legacy(qbObject)) {
    ## Assume if one phenotype has subset, they all do.
    if(length(qbObject$subset) >= 1) {
      if (!is.null(qbObject$subset[[pheno.names[1]]]))
        return(qbObject)
    }

    ## Sequence of iterations in iterdiag.
    iterdiag <- qb.get.legacy(qbObject, "iterdiag", sub = NULL, warn = warn, pheno.col = pheno.names, drop = FALSE)
    
    ## Need to reorder QTL in mainloci by niter, then chrom, then locus.
    mainloci <- qb.get.legacy(qbObject, "mainloci", sub = NULL, warn = warn, pheno.col = pheno.names, drop = FALSE)
    
    ## Subset for GxE interaction: order by niter, covar, chrom, locus.
    gbye <- qb.get.legacy(qbObject, "gbye", sub = NULL, warn = warn, pheno.col = pheno.names, drop = FALSE)
    
    ## Need to reorder QTL1, QTL2 in pairloci so that:
    ##      chrom1 < chrom2
    ##      locus1 < locus2 if chrom1 == chrom2
    pairloci <- qb.get.legacy(qbObject, "pairloci", sub = NULL, warn = warn, pheno.col = pheno.names, drop = FALSE)

    subsets <- list()
    for(pheno in names(iterdiag)) {
      subs <- list(iterdiag = seq(nrow(iterdiag[[pheno]])))
      if (!is.null(mainloci[[pheno]]))
        subs$mainloci <- order(mainloci[[pheno]]$niter, mainloci[[pheno]]$chrom, mainloci[[pheno]]$locus)
      if (!is.null(gbye[[pheno]]))
        subs$gbye <- order(gbye[[pheno]]$niter, gbye[[pheno]]$covar, gbye[[pheno]]$chrom, gbye[[pheno]]$locus)
    
      if (is.null(pairloci[[pheno]]))
        subs$pairloci <- list(order = 0, flip = 0, left = 0, right = 0)
      else {
        s <- list(flip = (pairloci[[pheno]]$chrom1 > pairloci[[pheno]]$chrom2) |
                  ((pairloci[[pheno]]$chrom1 == pairloci[[pheno]]$chrom2) &
                   (pairloci[[pheno]]$locus1 > pairloci[[pheno]]$locus2)))
        ## order chromosome/locus sets
        if (any(s$flip)) {
          s$left <- c("chrom1","locus1")
          s$right <- c("chrom2","locus2")
          if (!is.null(pairloci[[pheno]]$ad)) {
            ## f2
            s$left <- c(s$left,"ad","varad")
            s$right <- c(s$right,"da","varda")
          }
          n.p <- nrow(pairloci[[pheno]])
          p.chr <- n.p * (match(c("chrom1","chrom2"), names(pairloci[[pheno]])) - 1)
          p.pos <- n.p * (match(c("locus2","locus2"), names(pairloci[[pheno]])) - 1)
          s$order <- order(pairloci[[pheno]]$niter,
                           as.matrix(pairloci[[pheno]])[seq(n.p) + p.chr[1 + s$flip]],
                           as.matrix(pairloci[[pheno]])[seq(n.p) + p.pos[1 + s$flip]],
                           as.matrix(pairloci[[pheno]])[seq(n.p) + p.chr[2 - s$flip]],
                           as.matrix(pairloci[[pheno]])[seq(n.p) + p.pos[2 - s$flip]])
        }
        else {
          s$left <- s$right <- 0
          s$order <- seq(nrow(pairloci[[pheno]]))
        }
        subs$pairloci <- s
      }
      subs$region <- region
      subsets[[pheno]] <- subs
    }
    
    ## Now attach subset to qb object.
    qbObject$subset <- subsets
  }
  else
    qbObject$args$region <- region

  qbObject
}
##############################################################################
qb.split.chr <- function(qbObject, split = qb.mainmodes(qbObject, ...)$valleys,
                        ...)
{
  qbObject$args$split.chr <- split
  qbObject
}
##############################################################################
is.legacy <- function(qbObject) is.null(qbObject$args)
##############################################################################
qb.legacy <- function(qbObject, remove = FALSE, ...)
{
  ## Convert from several objects to one combined qb object.

  ## pull.grid used in scan.R, slice.R in summary, plot
  ##                   and for qb.inter, qb.scantwo.smooth.
  ## scanone, etc that use cross.name to pass to summary/plot: what is really needed?
  ##     mostly need pull.grid, mainloci, pairloci, iterdiag.
  ##     can be condensed with care (will take some time)
  ##     requires some rethinking of qb.intertwo, qb.inter, but not too much.
  ## Need to also check qb.hpdone, plot.qb. 
  ##
  ## Check what new$args can be dropped and what are missing from call.
  ## Make function qb.cross.class to get class(cross) (don't need args$cross).

  ## Check if already converted.
  if(!is.legacy(qbObject)) {
    if(is.null(qbObject$mcmc.samples)) {
      ## Upgrade recent qb objects, with only one trait.
      pheno.col <- qb.get(qbObject, "pheno.col")
      pheno.names <- qb.pheno.names(qbObject)[pheno.col]
      mcmc <- list()
      for(pheno in pheno.names) {
        mcmc[[pheno]] <- list()
        for(i in c("iterdiag","mainloci","pairloci","covariates","gbye"))
          mcmc[[pheno]][[i]] <- qb.get(qbObject, i, pheno.col = pheno)
      }
      if(length(pheno.names) > 1)
        mcmc$sigma <- qb.get(qbObject, "sigma", pheno.col = pheno.names[1])
      qbObject$mcmc.samples <- mcmc
      
      ## Remove tacked on MCMC runs.
      for(i in c("iterdiag","mainloci","pairloci","covariates","gbye"))
        qbObject[[i]] <- NULL
    }
    new <- qbObject
  }
  else {
    ## Old arguments now bundled as "args".
    new <- list(args = qbObject)
    class(new$args) <- "list"
    
    ## Only need region from old subset list.
    new$args$region <- new$args$subset$region
    new$args$subset <- NULL
    
    ## Drop objects not actually used, particularly large ones.
    for(i in c("yvalue","fixcoef","rancoef"))
      new$args[[i]] <- NULL
    
    ## Cross object. Ideally run through qb.genoprob already.
    cross <- qb.cross(qbObject, genoprob = FALSE, warn = FALSE)
    
    ## Reduce to phenotypes needed.
    ## Make sure names are preserved.
    phenos <- c(new$args$pheno.col, new$args$covar)
    for(i in c("sex","pgm")) {
      tmp <- match(i, tolower(names(cross$pheno)))
      if(!is.na(tmp) & is.na(match(tmp, phenos)))
        phenos <- c(phenos, tmp)
    }
    phenos <- names(cross$pheno)[phenos[phenos>0]]
    cross$pheno <- data.frame(cross$pheno[, phenos, drop = FALSE])
    names(cross$pheno) <- phenos
    
    ## Reset pheno.col and covar to reflect change.
    n.pheno <- length(new$args$pheno.col)
    new$args$pheno.col <- seq(n.pheno)
    new$args$covar <- if(new$args$nfixcov)
      n.pheno + seq(new$args$nfixcov)
    else
      0
    new$args$covar <- c(new$args$covar,
                        if(new$args$nrancov) {
                          n.pheno + new$args$nfixcov +
                            seq(new$args$nrancov)
                        }
                        else {
                          0
                        })
    
    ## Assign qb.genoprob attributes to args list if not there.
    defaults <- qb.genoprob.defaults(cross)
    for(i in names(defaults))
      if(is.null(new$args[[i]]))
        new$args[[i]] <- defaults[[i]]
    
    ## Attach cleaned cross object.
    new$cross.object <- clean(cross)
    
    ## External MCMC files now internal.
    ## Organized as list of lists for multiple trait extension.
    pheno.names <- qb.pheno.names(qbObject, cross)[new$args$pheno.col]

    ## Redo reorder if multiple traits, as probably done wrong.
    if(length(pheno.names) > 1)
      qbObject <- qb.reorder(qbObject)

    ## Get MCMC samples.
    mcmc <- list()
    for(pheno in pheno.names) {
      mcmc[[pheno]] <- list()
      for(i in c("iterdiag","mainloci","pairloci","covariates","gbye"))
        mcmc[[pheno]][[i]] <-
          qb.get.legacy(qbObject, i,, FALSE, pheno.col = pheno, ...)
    }
    if(length(pheno.names) > 1)
      mcmc$sigma <- qb.get.legacy(qbObject, "sigma",, FALSE, pheno.col = pheno.names[1], ...)
    new$mcmc.samples <- mcmc
    
    ## Set up number of iterations correctly.
    new$args$n.iter <- nrow(mcmc[[1]]$iterdiag)

    class(new) <- c("qb", "list")
  }
  if(remove)
    qb.remove(qbObject, external.only = TRUE)

  ## Set up chromosome split for multiple linked loci.
  ## qb.split.chr(new)
  
  new
}
##############################################################################
qb.genoprob.defaults <- function(cross)
{
  ## Assign calc.genoprob attributes to args list.
  defaults <- formals(calc.genoprob)
  defaults$cross <- NULL
  defaults$step <- 2
  defaults$stepwidth <- "variable"
  defaults$map.function <- "haldane"
  ## Other arguments are off.end, error.prob.

  ## Add extra argument to qb.genoprob.
  defaults$tolerance <- 1e-6

  tmp <- cross$geno[[1]]$prob
  if(!is.null(tmp)) {
    for(i in names(defaults)) {
      ## Use genoprob value if already present.
      ## Else use default.
      tmp2 <- attr(tmp, i)
      if(!is.null(tmp2))
        defaults[[i]] <- tmp2
    }
  }
  defaults
}
##############################################################################
qb.get <- function(qbObject, element,
                   pheno.col = qb.get(qbObject, "pheno.col", warn = warn)[1],
                   warn = TRUE, ...)
{
  ## Check if MCMC data element.
  is.mcmc.data <- element %in% c("iterdiag","mainloci","pairloci","covariates","gbye","sigma")
    
  if(!is.legacy(qbObject)) { ## new style qb object
    if(is.mcmc.data) {
      if(is.character(pheno.col)) {
        pheno.names <- names(qbObject$mcmc.samples)
        tmp <- match(pheno.col, pheno.names)
        if(is.na(tmp))
          stop(paste("cannot match pheno.col:", tmp))
        pheno.col <- tmp
      }

      x <- qbObject$mcmc.samples[[pheno.col]][[element]]
      if(is.null(x))
        x <- qbObject$mcmc.samples[[element]]
      if(is.null(x))
        x <- qbObject[[element]]
    }
    else { ## cross.object or element of args.
      if(element == "cross.object" | element == "args")
        x <- qbObject[[element]]
      else
        x <- qbObject$args[[element]]
    }
  }
  else { ## old style qb object.
    ## Should not get to this point after qb.legacy(),
    ## but kept for backward compatibility.
    if(warn)
      warning("Getting legacy objects: use qb.legacy() to update")
    x <- qb.get.legacy(qbObject, element, pheno.col = pheno.col,
                       warn = warn, ...)
  }
  x
}
##############################################################################
qb.find.pheno <- function(qbObject, cross = qb.cross(qbObject, genoprob = FALSE, warn = warn),
                          pheno.col, warn = FALSE)
{
  ## Match up pheno column as numeric or character.
  pheno.names <- qb.pheno.names(qbObject, cross)
  if(is.numeric(pheno.col))
    pheno.col <- pheno.names[pheno.col]
  pheno.names <- pheno.names[qb.get(qbObject, "pheno.col", warn = warn)]
  tmp <- match(pheno.col, pheno.names)
  if(any(is.na(tmp)))
    stop(paste("cannot match all pheno.col:", paste(pheno.col, collapse = ", ")))
  names(tmp) <- pheno.col
  attr(tmp, "pheno.names") <- pheno.names
  tmp
}
##############################################################################
qb.get.legacy <- function(qbObject, element,
                          sub = qbObject$subset,
                          warn = TRUE,
                          cross = qb.cross(qbObject, genoprob = FALSE, warn = warn),
                          pheno.col = qb.get(qbObject, "pheno.col", warn = warn),
                          drop = TRUE, legacy.dir = FALSE, ...)
{
  is.mcmc.data <- element %in% c("iterdiag","mainloci","pairloci","covariates","gbye","sigma")

  if(is.mcmc.data) { ## Get MCMC samples from flat files.

    ## Allow for multiple pheno.col; requires care on output.dir.
    pheno.col <- qb.find.pheno(qbObject, cross, pheno.col, warn = warn)
    pheno.names <- attr(pheno.col, "pheno.names")

    ## Get MCMC sample element. Return NULL if missing or empty.
    renamedElement<-paste(element,".dat",sep="")

    ## *** NEED TO FIX THIS FOR COMBINED FILES.
    if(legacy.dir)
      filename <- file.path(qbObject$output.dir[pheno.col],
                            renamedElement)
    else
      filename <- file.path(qbObject$output.dir[1],
                            renamedElement)
    if (!file.exists(filename))
      return(NULL)
    tmp <- scan(filename, n = 1, quiet = TRUE)
    x <- NULL
    if (length(tmp))
      x <- as.data.frame(read.table(filename))
    if (is.null(x)) 
      return(NULL)
    
    ## Specify variable names.
    is.bc <- (qbObject$cross == "bc")
    eff1 <- "add"
    eff2 <- "aa"
    var1 <- "varadd"
    var2 <- "varaa"
    if (!is.bc) {
      eff1 <- c(eff1,"dom")
      eff2 <- c(eff2,"ad","da","dd")
      var1 <- c(var1,"vardom")
      var2 <- c(var2,"varad","varda","vardd")
    }
    v <- var1
    if (qbObject$epistasis)
      v <- c(v,var2)
    if (qbObject$qtl_envi) {
      v <- c(v,"envadd")
      if (!is.bc) v <- c(v,"envdom")
    }
    if (qbObject$envi) v <- c(v,"varenv")
    v <- c(v,"var")

    xnames <- switch(element,
                     covariates = {
                       n.cov <- qb.get(qbObject, "nfixcov",warn=warn) + qb.get(qbObject,"nrancov",warn=warn)
                       if (n.cov)
                         paste("cov", seq(n.cov), sep = "")
                       else
                         c("cov1")
                     },
                     sigma = {
                       if (is.null(x))
                         c("cov1")
                       else
                         paste("cov", seq(ncol(x)), sep = "")
                     },
                     gbye = c("niter","n.gbye","covar","chrom","locus",eff1,var1),
                     mainloci = c("niter","nqtl","chrom","locus",eff1,var1),
                     pairloci = c("niter","n.epis","chrom1","locus1","chrom2","locus2",eff2,var2),
                     iterdiag = c("niter","nqtl","mean","envvar",v))

    ## assign the names to X
    if (is.null(x)) {
      if (!is.null(xnames)) {
        x <- data.frame(matrix(NA, 0, length(xnames)))
        names(x) <- xnames
      }
    }
    else { ## Element has MCMC samples.
      if(ncol(x) == length(xnames)) { ## Separate MCMC samples per phenotype.
        names(x) <- xnames
        x <- list(a = x)
        names(x) <- names(pheno.col)[1]
      }
      else {
        ## Second column (after niter) has phenotype number.
        ## Split MCMC sample element into a list.
        if(element == "covariates") {
          names(x) <- c("pheno.col", xnames)
          x <- split(x[, -1, drop = FALSE], x[, 1], drop = FALSE)
        }
        else {
          names(x) <- c(xnames[1], "pheno.col", xnames[-1])
          x <- split(x[, -2, drop = FALSE], x[, 2], drop = FALSE)
        }
        names(x) <- pheno.names
        x <- x[names(pheno.col), drop = FALSE]
      }

      ## Reorder for covariates is same as iterdiag.
      if (element == "covariates" | element == "sigma")
        element <- "iterdiag"

      if(!is.null(sub)) {
        if("iterdiag" %in% names(sub)) {
          if(length(pheno.names) > 1) {
            ## This should not happen.
            warning("making NULL reorder subset for qb object with multiple phenotypes")
            sub <- NULL
          }
          else {
            ## This could happen with old legacy stuff.
            sub <- list(a = sub)
            names(sub) <- pheno.names
          }
        }
        else
          names(sub) <- pheno.names

        ## Reorder and subset element by phenotype.
        for(i in names(x)) {
          mysub <- sub[[i]][[element]]
          if (element == "pairloci") {
            if (any(mysub$flip))
              x[[i]][mysub$flip, c(mysub$left,mysub$right)] <-
                x[[i]][mysub$flip, c(mysub$right,mysub$left)]
            x[[i]] <- x[[i]][mysub$order,, drop = FALSE]
          }
          else
            x[[i]] <- x[[i]][mysub,, drop = FALSE]
        }
      }
    }
    
    ## Assign genoprob attributes to args list.
    defaults <- qb.genoprob.defaults(cross)
    step <- qb.get(qbObject, "step", warn = warn)
    if(is.null(step))
      step <- defaults$step

    ## Pull grid (which automatically uses region subsets.
    grid <- pull.grid(qbObject, offset = TRUE, mask.region = FALSE,
                      cross = cross,
                      step = step,
                      off.end = defaults$off.end,
                      stepwidth = defaults$stepwidth,
                      warn = warn, ...)
    chrpos <- paste(grid$chr,
                    unlist(apply(as.matrix(table(grid$chr)), 1,
                                 function(x) {
                                   seq(0, length = x)
                                 })),
                    sep = ":")
    if (element == "mainloci" | element == "gbye") {
      for(i in names(x))
        x[[i]]$locus <- grid$pos[match(paste(x[[i]]$chrom, x[[i]]$locus, sep = ":"), chrpos)]
    }
    if (element == "pairloci") {
      for(i in names(x)) {
        x[[i]]$locus1 <- grid$pos[match(paste(x[[i]]$chrom1, x[[i]]$locus1, sep = ":"), chrpos)]
        x[[i]]$locus2 <- grid$pos[match(paste(x[[i]]$chrom2, x[[i]]$locus2, sep = ":"), chrpos)]
      }
    }
    ## Revert to data frame if only one pheno.col. This is consistent with old style.
    if(length(x) == 1 & drop)
      x <- x[[1]]
  }
  else { ## Old style qb elements.
    if(element == "region")
      x <- qbObject$subset$region
    else
      x <- qbObject[[element]]
  }
  x
}

##############################################################################
qb.remove <- function(qbObject, verbose = TRUE, external.only = FALSE)
{
  if(is.character(qbObject)) {
    qbName <- qbObject
    qbObject <- get(qbName)
  }
  else
    qbName <- deparse(substitute(qbObject))

  tmp <- qb.get(qbObject, "output.dir", warn = FALSE)
  if (any(dirname(tmp) != system.file("external", package = "qtlbim"))) {
    if (verbose)
      warning(paste("Removing external director", ifelse(length(tmp) == 1, "y ", "ies "),
                    paste(tmp, collapse = ", ")),
              call. = FALSE, immediate. = TRUE)
    for(i in tmp)
      unlink(i, recursive = TRUE)
  }
  if(!external.only) {
    if (verbose)
      warning(paste("Removing internal", qbName),
              call. = FALSE, immediate. = TRUE)
    remove(list = qbName, pos = 1)
  }
  invisible()
}
##############################################################################
qb.exists <- function(qbObject)
{
  if(!inherits(qbObject, "qb"))
    stop(paste("Object is not of class qb"))
  if(is.legacy(qbObject)) {
    warning("Object is legacy qb object; please use qb.legacy to upgrade object.")
    tmp <- qb.get(qbObject, "output.dir")
    if (any(!file.exists(tmp))) {
      cat("\n")
      stop(paste("Object contains no MCMC samples in", tmp), call. = FALSE)
    }
  }
  invisible()
}
##############################################################################
qb.recover <- function(cross, traitName,
                        ## Find existing directory name if it exists.
                        output.dir = system(paste("ls -d ./", traitName[1], "_*",
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
  if (length(output.dir) == 0)
    stop("No trait MCMC directories exist to be recovered")
  if (length(output.dir) > 1) {
    print(output.dir)
    stop("Multiple trait MCMC directories exist: select one as output.dir")
  }

  step <- attr(cross$geno[[1]]$prob, "step")
  if (is.null(step))
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

  ## Object is of class "qb".
  class(qbObject) <- "qb"

  ## Get the actual number of iterations.
  ## Need to first have qbObject with some settings before qb.get call.
  qbObject$n.iter <- n.iter <- nrow(qb.get(qbObject, "iterdiag", warn = FALSE))
  qbObject$n.burnin <- n.burnin

  ## Set up subset reordering.
  qbObject$subset <- NULL
  qb <- qb.reorder(qbObject)
  qb.legacy(qb)
}
##############################################################################
subset.qb <- function(x, nqtl = 1, pattern = NULL, exact = FALSE, chr,
                      region, offset = TRUE, restrict.pair = TRUE,
                      pheno.col = qb.get(x, "pheno.col"), ...)
{
  ## Checks.
  if (missing(nqtl) & missing(pattern) & missing(exact) & missing(chr) &
     missing(region))
    return(x)

  nqtl <- nqtl[1]
  if (!is.null(pattern)) {
    if (is.character(pattern))
      stop("pattern must be numeric, not character")
    nqtl <- max(nqtl, length(pattern))
  }
  
  ## Get elements (except covariate, which matches iterdiag)
  ## Note that previous subsetting is done here (see Incorporate below).
  pheno.col <- pheno.col[1]
  iterdiag <- qb.get(x, "iterdiag", pheno.col = pheno.col)
  mainloci <- qb.get(x, "mainloci", pheno.col = pheno.col)
  pairloci <- qb.get(x, "pairloci", pheno.col = pheno.col)
  gbye <- qb.get(x, "gbye", pheno.col = pheno.col)

  ## And set up subset list on nqtl.
  ## These are TRUE/FALSE indicators except for region.
  sub <- list()

  ## Subset on number of QTL (led by iterdiag).
  if (exact) {
    ## Only exactly nqtl QTL (see also pattern below).
    sub$iterdiag <- iterdiag$nqtl == nqtl
    if (!sum(sub$iterdiag))
      stop(paste("empty object: no iterations with number of QTL =", nqtl))
  }
  else {
    ## At least nqtl QTL (see also pattern below).
    sub$iterdiag <- iterdiag$nqtl >= nqtl
    if (!sum(sub$iterdiag))
      stop(paste("empty object: no iterations with number of QTL >=", nqtl))
  }
  ## Now adjust other MCMC elements based on iterdiag.
  iters <- iterdiag[sub$iterdiag,"niter"]
  sub$mainloci <- !is.na(match(mainloci$niter, iters))
  sub$pairloci <- !is.na(match(pairloci$niter, iters))
  sub$gbye <- !is.na(match(gbye$niter, iters))

  ## Subset on pattern of chromosomes.
  if (!is.null(pattern)) {
    mypat <- table(pattern)
    if (exact) {
      ## Create function to only exactly match pattern retained in subset.
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
    if (!sum(sub$iterdiag))
      stop(paste("empty object: no patterns like", pattern))
    iters <- iterdiag[sub$iterdiag, "niter"]
    sub$mainloci <- sub$mainloci & !is.na(match(mainloci$niter, iters))
    sub$gbye <- sub$gbye & !is.na(match(gbye$niter, iters))
    sub$pairloci <- sub$pairloci & !is.na(match(pairloci$niter, iters))
  }

  cross <- qb.cross(x, genoprob = FALSE)
  
  ## Subset of regions in chromosomes.
  sub$region <- qb.get(x, "region")
  if (!missing(region)) {
    region <- as.list(region)
    region$chr <- qb.find.chr(x, region$chr, cross = cross, sort.chr = FALSE)
    iters <- rep(TRUE, length(unique(mainloci$niter)))
    if (!offset)
      cross.map <- pull.map(cross)
    for(i in seq(length(region$chr))) {
      if (!offset) {
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
    if (!sum(sub$iterdiag))
      stop(paste("empty object: no regions like", as.data.frame(region)))

    iters <- iterdiag[sub$iterdiag, "niter"]
    sub$mainloci <- sub$mainloci & !is.na(match(mainloci$niter, iters))
    sub$gbye <- sub$gbye & !is.na(match(gbye$niter, iters))
    sub$pairloci <- sub$pairloci & !is.na(match(pairloci$niter, iters))

    if(restrict.pair) {
      ## Drop linked QTL outside of region.
      if (!is.null(mainloci)) {
        tmp <- rep(FALSE, nrow(mainloci))
        for(i in seq(length(region$chr))) {
          tmp <- tmp | (mainloci$chrom == region$chr[i] &
                        (mainloci$locus < region$start[i] |
                         mainloci$locus > region$end[i]))
        }
        sub$mainloci <- sub$mainloci & !tmp
      }
      if (!is.null(gbye)) {
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
          ## Drop if either chrom is chr and locus not in region.
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
    else { ## Still restrict pairloci to have at least one in region.
      if (!is.null(pairloci)) {
        tmp <- rep(TRUE, nrow(pairloci))
        for(i in seq(length(region$chr))) {
          ## Keep if neither chrom is chr or
          ## one of chrom is chr and its locus is in region.
          tmp <- tmp & ((pairloci$chrom1 != region$chr[i] &
                         pairloci$chrom2 != region$chr[i]) |
                        ((pairloci$chrom1 == region$chr[i] &
                          (pairloci$locus1 >= region$start[i] &
                           pairloci$locus1 <= region$end[i])) |
                         (pairloci$chrom2 == region$chr[i] &
                          (pairloci$locus2 >= region$start[i] &
                           pairloci$locus2 <= region$end[i]))))
        }
        sub$pairloci <- sub$pairloci & tmp
      }
    }
  }

  ## Subset of chromosomes.
  if (!missing(chr)) {
    n.chr <- length(cross$geno)
    chr <- qb.find.chr(x, chr, names(cross$geno))

    ## Subset all MCMC elements based on chr.
    kept <- !is.na(match(seq(n.chr), chr))
    sub$mainloci <- sub$mainloci & kept[mainloci$chrom]
    sub$pairloci <- if (restrict.pair)
      sub$pairloci & kept[pairloci$chrom1] & kept[pairloci$chrom2]
    else
      sub$pairloci & (kept[pairloci$chrom1] | kept[pairloci$chrom2])
    sub$gbye <- sub$gbye & kept[gbye$chrom]
    ## Arrange to drop other chr via pull.grid function.
    for(i in seq(n.chr)[-chr])
      sub$region[i,] <- c(i,0,-1)
  }

  if(is.legacy(x)) {
    ## Legacy: incorporate new subset on top of any previous subsetting.
    x$subset$iterdiag <- x$subset$iterdiag[sub$iterdiag]
    x$n.iter <- length(x$subset$iterdiag)
    x$subset$mainloci <- x$subset$mainloci[sub$mainloci]
    x$subset$gbye <- x$subset$gbye[sub$gbye]
    x$subset$pairloci$order <- x$subset$pairloci$order[sub$pairloci]
    x$subset$region <- sub$region
  }
  else {
    ## Subset MCMC samples.
    ## Now adjusted for multiple trait format.
    if(is.null(x$mcmc.samples))
      x <- qb.legacy(x)
    mcmc.pheno <- x$mcmc.samples[[pheno.col]]
    mcmc.pheno$iterdiag <- mcmc.pheno$iterdiag[sub$iterdiag, ]
    if(!is.null(mcmc.pheno$covariates))
      mcmc.pheno$covariates <-
        as.matrix(mcmc.pheno$covariates)[sub$iterdiag, ]
    mcmc.pheno$mainloci <- mcmc.pheno$mainloci[sub$mainloci, ]
    if(!is.null(mcmc.pheno$gbye))
      mcmc.pheno$gbye <- mcmc.pheno$gbye[sub$gbye, ]
    if(!is.null(mcmc.pheno$pairloci))
      mcmc.pheno$pairloci <- mcmc.pheno$pairloci[sub$pairloci, ]
    x$mcmc.samples[[pheno.col]] <- mcmc.pheno

    x$args$n.iter <- nrow(mcmc.pheno$iterdiag)
    x$args$region <- sub$region
  }
  x
}
##############################################################################
qb.save <- function(cross, qbObject, dir = ".", Name= substring(qbName, 3))
{
  .Deprecated("save")
  
  qbCross <- deparse(substitute(cross))
  qbName <- deparse(substitute(qbObject))

  my.mcmc <- paste(Name, "MCMC", sep = ".")
  tmp <- file.path(dir, my.mcmc)
  if (!file.exists(tmp))
    dir.create(tmp)
  for(i in c("iterdiag.dat", "mainloci.dat", "pairloci.dat", "gbye.dat",
             "covariates.dat"))
    file.copy(file.path(qbObject$output.dir, i), file.path(dir, my.mcmc, i),
              overwrite = TRUE)
  qbObject$output.dir <- my.mcmc
  assign(qbName, qbObject, envir = parent.frame())
  ## Clean cross to make it smaller (qb.load will rebuild).
  cross <- clean(cross)
  
  ## Save in right place.
  out <- exists(qbName) & exists(qbCross)
  if (out) {
    .tryResults <- try(save(list = c(qbCross,qbName),
                            file = file.path(dir, paste(Name, "RData", sep = "."))))
    out <- !(class(.tryResults) == "try-error")
  }
  invisible(out)
}
##############################################################################
qb.load <- function(cross, qbObject,
                    dir = system.file("external", package = "qtlbim"),
                    file = paste(Name, "RData", sep = "."))
{
  .Deprecated("load")
  
  crossName <- deparse(substitute(cross))
  qbName <- deparse(substitute(qbObject))
  Name <- substring(qbName, 3)
  file <- file.path(dir, file)
  out <- exists(qbName) & exists(crossName)
  if (!out) {
    .tryResults <- try(load(file = file, parent.frame()))
    out <- !(class(.tryResults) == "try-error")
    ## Check of cross and qbObject now exist.
    if (out & (!exists(qbName) | !exists(crossName)))
      out <- FALSE
    if (out) {
      qbObject$output.dir <-
        file.path(dir, basename(qb.get(qbObject, "output.dir", warn = FALSE)))

      ## Assign qbObject in legacy format.
      assign(qbName, qbObject, envir = parent.frame())

      ## Set step if not done already.
      step <- qb.get(qbObject, "step", warn = FALSE)
      if (is.null(step))
        step <- 2

      ## Assign cross object.
      cross <- qb.genoprob(cross, step = 2)
      assign(crossName, cross, envir = parent.frame())

      ## Re-assign qbObject after updating legacy object.
      assign(qbName, qb.legacy(qbObject), envir = parent.frame())
    }
  }
  invisible(out)
}
##############################################################################
qb.niter <- function(qbObject)
  qb.get(qbObject, "n.iter")
##############################################################################
qb.nqtl <- function(qbObject,
                    iterdiag = qb.get(qbObject, "iterdiag", ...),
                    mainloci = qb.get(qbObject, "mainloci", ...),
                    match.iter = TRUE, ...)
{
  iterdiag.nqtl <- rep(0, nrow(iterdiag))
  if(!is.null(mainloci)) if(nrow(mainloci)) {
    ## Fix nqtl in samples:
    ## iterdiag[, "nqtl"] may be wrong due to subsetting earlier.
    tmp <- table(mainloci[, "niter"])
    tmp2 <- match(names(tmp), as.character(iterdiag[, "niter"]), nomatch = 0)
    iterdiag.nqtl[tmp2] <- c(tmp[tmp2 > 0])
    tmp <- sum(tmp == 0)
    if(tmp)
      warning(paste(tmp, "mismatches of mainloci iterations and iterdiag"))
    
    if (!match.iter) {
      iterdiag.nqtl <- iterdiag.nqtl[iterdiag.nqtl > 0]
    }
  }
  iterdiag.nqtl
}
##############################################################################
## qb.cross 
##          This function is used to supply the default argument for
##  a cross in the functions
##                    "plot.qb"
##                    "qb.pattern"
##                    "qb.BayesFactor"
##                    "subset.qb"
##                    "qb.scanone"
##                    "qb.scantwo"
## The cross is extracted from the options stored in the qb argument.
##                    
## arguments
##         qbObject         An object of class qb.
##
## returns
##         An object of class "f2" (inheriting from class "cross").
##
## errors/exceptions:
##         If the name "cross" is not found in the options object returned by
## qb.get, then the "stop" function is called.
##

qb.cross <- function(qbObject, genoprob = TRUE, ...)
{
  ## Try to get cross object imbedded in qb object.
  cross <- qb.get(qbObject, "cross.object", ...)
  
  if(!is.null(cross)) {
    if(genoprob & is.null(cross$geno[[1]]$prob)) {
      ## Update genotype probabilities if needed.
      cross <- qb.genoprob(cross, step = qb.get(qbObject, "step"),
                           error.prob = qb.get(qbObject, "error.prob"),
                           off.end = qb.get(qbObject, "off.end"),
                           map.function = qb.get(qbObject, "map.function"),
                           stepwidth = qb.get(qbObject, "stepwidth"))
    }
    cross
  }
  else { ## Legacy qb object.
    rm(cross)
    cross.name <- qb.get(qbObject, "cross.name", ...)
    if (is.null(cross.name))
      stop("need to have cross.name as character string in qb object")
    get(cross.name)
  }
}

qb.cross.class <- function(qbObject)
{
  if(is.legacy(qbObject))
    qb.get(qbObject, "cross")
  else
    class(qbObject$cross.object)[1]
}
                           
##############################################################################
## qb.cex
##     This private function is used only once to set the default plotting
## parameter cex in the function "plot.qb.effects".
##
## arguments
##     qbObject    An object of class qb.
##
##     min.cex     A minimum (floating point) value by which symbols and
##                 text should be scaled relative to the default.
## returns
##     The return value is used in "plot.qb.effects", to set the cex
## parameter in the "plot" function.  This is a scale factor by which
## symbols and text will be scaled relative to the default. Note: the
## "par" function for setting plotting parameters also has an argument
## called cex which behaves differently.
## 

qb.cex <- function(qbObject, min.cex = 3.85)
{
  tmp <- qb.niter(qbObject)
  if (tmp)
    2 ^ (2 - min(min.cex, max(2, (log10(tmp)))))
  else
    1
}
##############################################################################
split.pattern <- function(patterns, epistasis = FALSE)
{
  ## WORK IN PROGRESS.
  chrs <- apply(as.matrix(patterns), 1,
                  function(x) table(strsplit(x, ",", fixed = TRUE)))
  if(epistasis) {
    epis <- lapply(chrs,
                   function(x)
                   {
                     tmp <- grep(":", names(x))
                     if(length(tmp)) {
                       out <- as.data.frame(strsplit(names(x[tmp]), ":", fixed = TRUE))
                       names(out) <- names(x[tmp])
                       out
                     }
                     else
                       numeric(0)
                   })
    names(chrs) <- names(epis) <- patterns
  }
  else
    epis <- NULL
  chrs <- lapply(chrs,
                 function(x)
                 {
                   tmp <- grep(":", names(x))
                   if(length(tmp))
                     x[-tmp]
                   else
                     x
                 })
  names(chrs) <- names(epis) <- patterns
  list(chr = chrs, epi = epis)
}
##############################################################################
qb.match.pattern <- function(qbObject, targets, exact = TRUE,
                             patterns = qb.makepattern(qbObject, ...),
                             ...)
{
  ## Change vector of character patterns (chr or epi separated by ",")
  ## into matrix with chr as row, pattern as column.
  ## Entry in matrix is number of copies of each chr.
  ## This handles split chr correctly.
  depat <- function(pattern) {
    strs <- apply(as.matrix(pattern), 1,
                  function(x) table(strsplit(x, ",", fixed = TRUE)))
    if(is.matrix(strs)) {
      strs <- as.data.frame(strs)
      levs <- row.names(strs)
      n.strs <- ncol(strs)
      levals <- rep(levs, n.strs)
      n.strs <- rep(seq(n.strs) - 1, rep(length(levs), n.strs))
    }
    else {
      levals <- unlist(lapply(strs, names))
      levs <- sort(unique(levals))
      tmp <- sapply(strs, length)
      n.strs <- rep(seq(length(strs)) - 1, tmp)
    }
    out <- matrix(0, length(levs), length(strs))
    mat <- unlist(lapply(levals, function(x, y) match(x, y), levs))
    out[mat + length(levs) * n.strs] <- unlist(strs)
    dimnames(out) <- list(levs, pattern)
    out
  }
  pat <- depat(patterns)
  tar <- depat(targets)

  ## Find matches of chr names.
  tmp <- match(dimnames(tar)[[1]], dimnames(pat)[[1]], nomatch = 0)

  res <- matrix(FALSE, length(patterns), length(targets))
  dimnames(res) <- list(patterns, targets)

  if(all(tmp == 0)) ## No matches at all. Unlikely.
    return(res)

  ## Contract targets to elements matching patterns.
  if(any(tmp == 0)) {
    nottar <- apply(tar[tmp == 0, ], 2, sum) == 0
    tar <- tar[tmp > 0,, drop = FALSE]
  }
  else
    nottar <- rep(TRUE, ncol(tar))

  ## Contract patterns to elements matching targets.
  if(length(tmp) < nrow(pat))
    notpat <- apply(pat[-tmp,, drop = FALSE], 2, sum) == 0
  else
    notpat <- rep(TRUE, ncol(pat))
  pat <- pat[tmp,, drop = FALSE]
  
  patfn <- if(exact)
    function(x, target) all(target == x)
  else
    function(x, target) all(target <= x)

  ## Inner product comparison. Must be a more elegant way, but I forget.
  for(i in targets) {
    res[, i] <- apply(pat, 2, patfn, tar[,i])
    ## Need to check others that are not in pat.
    if(exact)
      res[, i] <- res[, i] & notpat
  }
  ## Check if something in targets is not in patterns.
  for(i in patterns)
    res[i, ] <- res[i, ] & nottar
  
  res
}
##############################################################################
qb.split.names <- function(geno.names, locus, split.chr)
{
  if(length(split.chr)) { ## Set up names for split chr.
    ## Want a way to replace geno.names[mainloci$chrom] with
    ## the split.chr name. Need to use mainloci$locus to see what split on chr.
    ##
    ## Need to think about how this affects target for qb.close?
    ## Also need to think about priors for BayesFactor.
    ##
    ## Dumb loop on number of splits per chr.
    ## There has got to be a more clever way.
    for(i in names(split.chr)) {
      ii <- (i == geno.names)
      if(any(ii)) {
        suffix <- rep(1, sum(ii))
        for(j in seq(length(split.chr[[i]])))
          suffix[split.chr[[i]][j] < locus[ii]] <- j + 1
        geno.names[ii] <- paste(geno.names[ii], suffix, sep = ".")
      }
    }
  }
  geno.names
}
##############################################################################
qb.makepattern <- function(qbObject, epistasis = TRUE,
                           geno.names = qb.geno.names(qbObject),
                           iterdiag = qb.get(qbObject, "iterdiag", ...),
                           mainloci = qb.get(qbObject, "mainloci", ...),
                           pairloci = qb.get(qbObject, "pairloci", ...),
                           ...)
{
  ## Make pattern with chr separated by commas, epistatic pairs joined by colon.
  out <- rep("NULL", nrow(iterdiag))
  names(out) <- iterdiag$niter

  split.chr <- qb.get(qbObject, "split.chr")
  split.chr <- split.chr[names(split.chr) %in% geno.names]

  split.names <- qb.split.names(geno.names[mainloci$chrom],
                                mainloci$locus, split.chr)

  tmp <- unlist(tapply(split.names,
                       mainloci$niter,
                       paste, collapse = ",", sep = ""))
  out[names(tmp)] <- tmp

  if (epistasis) if(!is.null(pairloci)) {
    split.names <- qb.split.names(geno.names[pairloci$chrom1],
                                  pairloci$locus1, split.chr)
    split.names2 <- qb.split.names(geno.names[pairloci$chrom2],
                                   pairloci$locus2, split.chr)
    
    tmp <- unlist(tapply(paste(split.names, split.names2, sep = ":"),
                         pairloci$niter,
                         paste, collapse = ",", sep = ""))
    out[names(tmp)] <- paste(out[names(tmp)], tmp, sep = ",")
  }
  
  out
}
##############################################################################
qb.find.chr <- function(qbObject, chr = NULL,
                        geno.names = qb.geno.names(qbObject, cross),
                        cross = qb.cross(qbObject, genoprob = FALSE),
                        sort.chr = TRUE, drop.duplicate = TRUE, warn = FALSE)
{
  ## chr may be passed as logical, numeric or character.
  ## Internally we use numeric index to chromosomes
  ## as this coincides with how MCMC samples are stored.

  ## This is adapted from R/qtl's subset.cross.
  n.chr <- length(geno.names)
  if(is.null(chr))
    chr <- seq(n.chr)
  else {
    if (is.logical(chr)) {
      if (length(chr) != n.chr) 
        stop("If logical, chr argument must have length ", 
             n.chr)
      chr <- (1:n.chr)[chr]
      if(!length(chr))
        stop("No chr selected")
    }
    else {
      if (is.numeric(chr)) {
        if (all(chr < 1)) {
          if(all(chr) >= -n.chr)
            chr <- (1:n.chr)[chr]
        }
        else chr <- sort(chr)
        if (any(chr < 1 | chr > n.chr)) 
          stop("Chromosome numbers out of range.")
      }
      else { ## is.character
        if (any(!(chr %in% geno.names))) 
          stop("Not all chromosome names found.")
        chr <- match(chr, geno.names)
      }
      if (length(chr) != length(unique(chr)) & drop.duplicate) {
        chr <- unique(chr)
        if(warn)
          warning("Dropping duplicate chromosomes")
      }
    }
  }
  if(sort.chr)
    chr <- sort(chr)
  chr
}
##############################################################################
qb.pheno.names <- function(qbObject,
                           cross = qb.cross(qbObject, genoprob = FALSE))
  names(cross$pheno)
##############################################################################
qb.geno.names <- function(qbObject,
                          cross = qb.cross(qbObject, genoprob = FALSE))
  names(cross$geno)
##############################################################################
qb.nloci <- function(qbObject,
                     cross = qb.cross(qbObject, genoprob = FALSE),
                     step = qb.get(qbObject, "step"),
                     off.end = qb.get(qbObject, "off.end"),
                     stepwidth = qb.get(qbObject, "stepwidth"))
  length(unlist(pull.loci(cross, step, off.end, stepwidth)))
##############################################################################
pull.loci <- function(cross,
                      step = attr(cross$geno[[1]]$prob, "step"),
                      off.end = attr(cross$geno[[1]]$prob, "off.end"),
                      stepwidth = attr(cross$geno[[1]]$prob, "stepwidth"),
                      region = NULL)
{
  ## Need to check step, off.end, stepwidth.

  tmpfn <- function(x, step, off.end, stepwidth) {
    create.map(x$map, step, off.end, stepwidth)
  }
  loci <- lapply(cross$geno, tmpfn, step, off.end, stepwidth)
  if(!is.null(region)) {
    for(i in region$chr)
      loci[[i]] <- loci[[i]][loci[[i]] >= region$start[i] - 0.1 &
                             loci[[i]] <= region$end[i] + 0.1]
  }
  class(loci) <- "map"
  loci
}
##############################################################################
pull.grid <- function (qbObject, offset = FALSE, spacing = FALSE,
                       mask.region = TRUE,
                       cross = qb.cross(qbObject, genoprob = FALSE, ...),
                       step = qb.get(qbObject, "step"),
                       off.end = qb.get(qbObject, "off.end"),
                       stepwidth = qb.get(qbObject, "stepwidth"),
                       ...) 
{
  ## Pull grid map of loci.
  if(mask.region)
    region <- qb.get(qbObject, "region", ...)
  else
    region <- NULL

  grid.map <- pull.loci(cross, step, off.end, stepwidth, region)

  pos <- unlist(grid.map)
  len <- sapply(grid.map, length)

  ## Pull map.
  cross.map <- pull.map(cross)

  ## Construct map position with optional offset from 0 start.
  if (!offset) {
    m <- sapply(cross.map, function(x) x[1])
    pos <- pos - rep(m, len)
  }

  ## Construct grid object with chr as first column.
  grid <- data.frame(chr = rep(seq(grid.map), len),
                     row.names = paste("c", names(pos), sep = ""))

  if (spacing) {
    ## If spacing, add columns for map (=pos), eq.spacing, xchr.
    ## This is used only to create scantwo object in qb.scantwo, plot.qb.scantwo.
    grid$map <- pos
    nmap <- length(pos)
    grid$eq.spacing <- unlist(lapply(grid.map, function(x) {
      lx <- length(x)
      if (lx) {
        ## kludge to determine equal spacing
        d <- diff(x)
        tbl <- table(d)
        maxtbl <- max(tbl)
        dmode <- as.numeric(names(tbl)[tbl == maxtbl])
        if (maxtbl * 2 > lx)
          c(1, d == dmode)
        else
          rep(0, lx)
      }
      else
        integer()
    }))
    xclass <- sapply(cross.map, attr, "class")
    grid$xchr <- rep(xclass == "X", len)
  }
  else {
    ## Otherwise second column is pos.
    grid$pos <- pos
  }
  grid
}
