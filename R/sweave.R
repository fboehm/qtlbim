#####################################################################
##
## $Id: sweave.R,v 1.15.2.4 2006/09/07 21:15:39 byandell Exp $
##
##     Copyright (C) 2006 Brian S. Yandell
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
qb.sweave <- function(cross, pheno.col = 1, n.iter = 3000, n.draws = 64,
                      scan.type = "2logBF", hpd.level = 0.5,
                      upper.threshold = thr,
                       SweaveFile = system.file("doc", "hyperslide.Rnw",
                         package = "qtlbim"),
                       SweaveExtra = NULL,
                       PDFDir = paste(pheno.name, "PDF", sep = ""),
                       remove.qb = TRUE)
{
  .qb.name <- deparse(substitute(cross))
  assign(".qb.name", .qb.name, ".GlobalEnv")
  cross <- subset(clean(cross),
                  chr = ("X" != unlist(lapply(cross$geno, class))))

  assign(".qb.cross", cross, ".GlobalEnv")

  if(any(!is.numeric(pheno.col)))
    pheno.col <- match(pheno.col, names(cross$pheno))
  if(is.na(pheno.col))
    stop("invalid pheno.col")
  assign(".qb.pheno", pheno.col, ".GlobalEnv")
  pheno.name <- qb.pheno.names( , cross)[pheno.col]

  assign(".qb.niter", n.iter, ".GlobalEnv")
  assign(".qb.draws", n.draws, ".GlobalEnv")

  thr <- if(class(cross)[1] == "bc")
    2
  else ## f2
    4
  assign(".qb.threshold", c(upper=upper.threshold), ".GlobalEnv")

  assign(".qb.scan.type", scan.type, ".GlobalEnv")
  assign(".qb.hpd.level", hpd.level, ".GlobalEnv")

  ## PDF directory.
  assign(".qb.PDFDir", PDFDir, ".GlobalEnv")
  if(!file.exists(PDFDir))
      dir.create(PDFDir)

  assign(".qb.SweaveFile", SweaveFile, ".GlobalEnv")
  assign(".qb.SweaveExtra", SweaveExtra, ".GlobalEnv")

  assign(".qb.remove", remove.qb, ".GlobalEnv")

  require("tools")
  Sweave(SweaveFile)
}
#################################################3
## Set up genetic architecture as list.
##
##   arch$qtl          data frame with chr and pos.
##   arch$pair.by.qtl  data frame with epistatic pairs indexed to qtl.
##   arch$pair.by.chr  list with chr and pos for epistatic pairs.
##   arch$chr.by.set   list with epistasis-connected chr in sets.
#################################################3
qb.arch <- function(object, ...) UseMethod("qb.arch")
#################################################3
qb.arch.default <- function(object, chr, pos, tolerance = 10, ...)
{
  ## Merge main QTL and epistatic QTL into one list.
  ## Make any QTL within tolerance of each other the same.
  ## Take some care with triplets, etc.
  n.pair <- nrow(object)
  if(is.null(n.pair) | length(n.pair) == 0)
    n.pair <- 0
  type <- rep(0, length(chr))
  if(n.pair) {
    type <- c(type, rep(seq(n.pair), rep(2, n.pair)))
    chr <- c(chr, t(object[, c("chr1","chr2")]))
    pos <- c(pos, t(object[, c("u.pos1","u.pos2")]))
  }
  create.arch(chr, pos, type, n.pair, tolerance)
}
#################################################3
create.arch <- function(chr, pos, type, n.pair, tolerance = 10)
{
  o <- order(chr, pos)
  d <- abs(diff(pos[o])) <= tolerance & diff(unclass(factor(chr[o]))) == 0
  newpos <- unlist(tapply(pos[o], cumsum(1 - c(0, d)),
                          function(x) rep(round(mean(x), 2), length(x))))
  data <- data.frame(chr = chr[o], pos = newpos, type = type[o])
  data$unique <- !duplicated(interaction(chr[o], newpos))
  arch <- list(qtl=data.frame(chr = data$chr[data$unique],
                 pos = data$pos[data$unique]))
  
  if(n.pair) {
    data$qtl <- cumsum(data$unique)
    pairs <- tapply(data$qtl, data$type, function(x) x)
    pairs[[1]] <- NULL
    ## Drop pairs that are the same (i.e. closer than width).
    pairs <- pairs[unlist(lapply(pairs,function(x)x[1]!=x[2]))]
    pairs <- t(as.data.frame(pairs))
    if(length(pairs)) {
      dimnames(pairs) <- list(NULL, paste("qtl", c("a","b"), sep = ""))
      arch$pair.by.qtl <- as.data.frame(pairs)
    }
  }
  qb.archalt(arch)
}
#################################################3
qb.arch.qb.BestPattern <- function(object, ...)
{
  ## Merge main QTL and epistatic QTL into one list.
  ## Take some care with triplets, etc.
  patterns <- split.pattern(object$patterns, epistasis = TRUE)
  pairs <- as.matrix(patterns$epi[[1]])
  n.pair <- ncol(pairs)
  if(is.null(n.pair) | length(n.pair) == 0)
    n.pair <- 0

  best <- object$model[[1]] 
  chr <- best$chrom
  pos <- best$locus
  
  type <- rep(0, length(chr))
  if(n.pair) {
    type <- c(type, rep(seq(n.pair), rep(2, n.pair)))
    tmp <- match(t(pairs),chr)
    chr <- ordered(c(chr, chr[tmp]), levels(chr))
    pos <- c(pos, pos[tmp])
  }
  create.arch(chr, pos, type, n.pair)
}
#################################################3
qb.archalt <- function(arch)
  {

  ## Alternate organizations of epistatic pairs. 
  if(!is.null(arch$pair.by.qtl)) {
    arch$pair.by.chr <- qb.archpairs(arch)
    arch$chr.by.set <- qb.pairgroup(arch, arch$pair.by.chr)
  }
  class(arch) <- c("qb.arch", "list")
  arch
}
#################################################3
qb.arch.step.fitqtl <- function(object, main = numeric(),
                                epistasis = data.frame(), ...)
{
  ## The object must be result of call to step.fitqtl.
  arch <- object$arch
  if(missing(main) & missing(epistasis))
    return(arch)

  ## pick up main only qtl
  n.main <- length(main)
  mains <- rep(0, n.main)
  if(n.main) {
    for(i in seq(n.main)) {
      tmp <- seq(nrow(arch$qtl))[arch$qtl$chr == main[i]]
      wh <- which.min(object$fit$result.drop[tmp, "Pvalue(F)"])
      if(any(wh))
        mains[i] <- tmp[wh]
    }
  }
  ## pick up epistatic pairs
  epistatis <- as.matrix(epistasis)
  n.epi <- nrow(epistasis)
  if(n.epi) {
    epis <- rep(0, n.epi)
    for(i in seq(n.epi)) {
      pair <- sort(unlist(epistasis[i, ]))
      wh <- (arch$qtl[as.character(arch$pair.by.qtl$q1), "chr"] == pair[1] &
             arch$qtl[as.character(arch$pair.by.qtl$q2), "chr"] == pair[2])
      if(any(wh))
        epis[i] <- seq(nrow(arch$pair.by.qtl))[wh]
    }
    mains <- sort(unique(c(mains,
                           match(unlist(arch$pair.by.qtl[epis, ]),
                                 row.names(arch$qtl)))))
    arch$pair.by.qtl <- arch$pair.by.qtl[epis,]
  }
  else
    arch$pair.by.qtl <- NULL

  ## Include only QTL selected as main or epistatic.
  arch$qtl <- arch$qtl[mains,]

  qb.archalt(arch)
}
#################################################3
qb.archpairs <- function(arch)
{
  if(is.null(arch$pair.by.qtl))
    return(NULL)
  out <- list()

  ## Make sure chr is character.
  tmp <- lapply(arch$pair.by.qtl,
                function(x,y) y[as.character(x), "chr"],
                arch$qtl)
  out$chr <- cbind(as.character(tmp[[1]]),
                   as.character(tmp[[2]]))

  n.pair <- nrow(out$chr)
  dimnames(out$chr) <- list(if(n.pair) paste("pair", seq(nrow(out$chr)))
                            else NULL,
                            paste("chr", c("a", "b"), sep = ""))
  out$pos <-
    as.data.frame(lapply(arch$pair.by.qtl,
                         function(x,y) y[as.character(x), "pos"],
                         arch$qtl))
  names(out$pos) <- paste("pos", c("a", "b"), sep = "")
  row.names(out$pos) <- paste("pair", seq(nrow(out$pos)))
  out
}
##############################################################
qb.pairgroup <- function(arch, pairs = qb.archpairs(arch))
{
  ## Chromosomes are now ordered factors. Change to character.
  n.pair <- nrow(pairs$chr)
  pairchr <- as.matrix(pairs$chr)
  
  ## All chr in pairs.
  cliques <- list()
  if(!is.null(n.pair)) {
    cliques[[1]] <- sort(unique(unlist(pairchr[1,])))
    if(n.pair > 1) {
      for(i in 2:n.pair) {
        tmp <- sort(unique(unlist(pairchr[i,])))
        ci <- 0
        j <- 0
        while(j < length(cliques) & ci == 0) {
          j <- j + 1
          if(any(match(tmp, cliques[[j]], nomatch = 0))) {
            ci <- j
            cliques[[j]] <- sort(unique(c(cliques[[j]], tmp)))
          }
        }
        if(ci == 0) {
          j <- 1 + length(cliques)
          cliques[[j]] <- tmp
        }
      }
    }
  }
  if(!length(cliques))
    cliques <- NULL
  else
    names(cliques) <- seq(cliques)
  cliques
}
#################################################3
summary.qb.arch <- function(object, ...)
{
  cat("main QTL loci:\n")
  print(t(object$qtl))
  if(!is.null(object$pair.by.qtl)) {
    cat("\nEpistatic pairs by qtl, chr, pos:\n")
    print(cbind(object$pair.by.qtl,
                object$pair.by.chr$chr,
                object$pair.by.chr$pos))

    cat("Epistatic chromosomes by connected sets:\n")
    cat(paste(sapply(object$chr.by.set,paste, collapse = ","),
              collapse = "\n"),
        "\n")
  }
  else
    cat("\nepistatic pairs: none\n")
  invisible()
}
#################################################3
print.qb.arch <- function(x, ...) summary(x, ...)
#################################################3
step.fitqtl <- function(cross, qtl, pheno.col = 1, arch,
                        cutoff = 0.05, trace = 1, steps = 100)
{
  main <- as.numeric(row.names(arch$qtl))
  epistasis <- arch$pair.by.qtl
  
  n.main <- length(main)
  if(!n.main)
    return(NULL)
  
  if(!is.null(epistasis)) {
    q1 <- epistasis[,1]
    q2 <- epistasis[,2]
  }
  else {
    q1 <- q2 <- numeric()
  }

  ## Set up initial formula.
  my.main <- paste("Q", main, sep = "")
  if(length(q1))
    my.epis <- paste("Q", q1, ":Q", q2, sep = "")
  else
    my.epis <- NULL
  my.formula <-
    formula(paste("y ~", paste(c(my.main, my.epis), collapse = "+")))

  ## Fit model.
  ## R/qtl fitqtl changed with 1.08-43, but not released yet.
  old.fitqtl <- compareVersion(qtlversion(), "1.08-43") < 0
  if(old.fitqtl)
    myfitqtl <- function(cross, pheno.col, ...)
      fitqtl(cross$pheno[[pheno.col]], ...)
  else
    myfitqtl <- function(cross, pheno.col, ...)
      fitqtl(cross, pheno.col, ...)
  
  cross.fit <- myfitqtl(cross, pheno.col, qtl,
                      formula = my.formula)
  if(trace == 3 & !is.null(cross.fit$result.drop))
    print(signif(cross.fit$result.drop[, -6], 3))
  else if(trace == 4)
    print(summary(cross.fit))

  ## Set up check list for loop.
  if(length(q1)) {
    check.main <- seq(n.main)[is.na(match(main, c(q1, q2)))]
    check.epis <- n.main + seq(length(q1))
  }
  else {
    check.main <- seq(n.main)
    check.epis <- NULL
  }
  checks <- c(check.main, check.epis)
  if(length(checks) > 1)
    max.p <- max(cross.fit$result.drop[checks, "Pvalue(F)"])
  else
    max.p <- 0

  ## Drop least significant epistatic or main only terms sequentially.
  step <- 0
  if(trace)
    drops <- list(drop=character(), LOD=numeric(), p=numeric())
  while(max.p > cutoff & step < steps & length(checks) > 1) {
    step <- step + 1
    wh <- which.max(cross.fit$result.drop[checks, "Pvalue(F)"])
    if(trace) {
      drops$drop[step] <- row.names(cross.fit$result.drop)[checks[wh]]
      drops$LOD[step] <- cross.fit$result.drop[checks[wh], "LOD"]
      drops$p[step] <- cross.fit$result.drop[checks[wh], "Pvalue(F)"]
      if(trace > 1) {
        cat("Step", step, ":", drops$drop[step],
            "LOD =", signif(drops$LOD[step], 3),
            "p =", signif(drops$p[step], 3), "\n")
      }
    }

    ## Adjust formula.
    l.main <- length(check.main)
    if(wh <= l.main) {
      main <- main[-check.main[wh]]
      n.main <- n.main - 1
      my.main <- paste("Q", main, sep = "")
    }
    else {
      q1 <- q1[l.main - wh]
      q2 <- q2[l.main - wh]
      if(length(q1))
        my.epis <- paste("Q", q1, ":Q", q2, sep = "")
      else
        my.epis <- NULL
    }
    my.formula <-
      formula(paste("y ~", paste(c(my.main, my.epis), collapse = "+")))

    ## Refit model.
    cross.fit <- myfitqtl(cross, pheno.col, qtl,
                        formula = my.formula)

    ## Adjust check list for next loop.
    if(length(q1)) {
      check.epis <- n.main + seq(length(q1))
      check.main <- seq(n.main)[is.na(match(main, c(q1, q2)))]
    }
    else {
      check.epis <- NULL
      check.main <- seq(n.main)
    }
    checks <- c(check.main, check.epis)
    if(length(checks) > 1) {
      max.p <- max(cross.fit$result.drop[checks, "Pvalue(F)"])

      if(trace == 3)
        print(signif(cross.fit$result.drop[, -6], 3))
      else if(trace == 4)
        print(summary(cross.fit))
    }
  }

  ## Summary reporting if trace is on
  if(trace) {
    if(length(drops$LOD)) {
      drops <- as.data.frame(drops)
      drops$LOD <- signif(drops$LOD, 3)
      drops$p <- signif(drops$p, 3)
      print(drops, right = FALSE)
    }
    if(trace > 1) {
      tmp <- cross.fit
      tmp$result.full <- signif(tmp$result.full[, c(1,2,4,5,7)], 3)
      if(!is.null(tmp$result.drop))
        tmp$result.drop <- signif(tmp$result.drop[, -6], 3)
      print(summary(tmp))
    }
  }

  ## Save final genetic architecture.
  arch$qtl <- arch$qtl[ as.character(main), ]
  if(length(q1))
    arch$pair.by.qtl <- data.frame(q1 = q1, q2 = q2)
  else
    arch$pair.by.qtl <- NULL
  
  out <- list(fit = cross.fit, arch = qb.archalt(arch))
  class(out) <- c("step.fitqtl", "list")
  out
}
#################################################3
anova.step.fitqtl <- function(object, object2, ...)
{
  fit1 <- summary(object$fit)$result.full
  if(missing(object2))
    return(print(fit1, quote = FALSE, na.print = ""))
  
  fit2 <- summary(object2$fit)$result.full
  if(fit1[1] < fit2[1]) {
    tmp <- fit1
    fit1 <- fit2
    fit2 <- tmp
  }
  difit <- fit1 - fit2
  difit["Error", ] <- fit1["Error", ]
  difit["Model", "MS"] <- difit["Model", "SS"] / difit["Model", "df"]
  stat <- difit["Model", "MS"] / difit["Error", "MS"]
  df1 <- difit["Model", "df"]
  df2 <- difit["Error", "df"]
  difit["Model", "Pvalue(Chi2)"] <- signif(1 - pchisq(stat * df1, df1), 3)
  difit["Model", "Pvalue(F)"] <- signif(1 - pf(stat, df1, df2), 3)
  print(difit, quote = FALSE, na.print = "")
}
