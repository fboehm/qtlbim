#####################################################################
##
## $Id: slice.R,v 1.12.2.11 2006/12/01 19:59:09 byandell Exp $
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
qb.sliceone <- function(qbObject, slice, epistasis = TRUE,
                        scan = c("main", "GxE", "epistasis"),
                        type.scan = type.scans,
                        covar = if(nfixcov) seq(nfixcov) else 0,
                        adjust.covar = NA,
                        chr = NULL,
                        sum.scan = "yes",
                        min.iter = 1,
                        aggregate = TRUE,
                        smooth = 3,
                        weight = c("sqrt","count","none","atten","ratten"),
                        split.chr = attr(qbObject, "split.chr"),
                        center.type = c("mode","mean","scan"),
                        verbose = FALSE, ...)
{
  type.scans <- c("heritability","LPD","LR","deviance","detection",
             "variance","estimate","cellmean","count","log10",
             "posterior","logposterior","2logBF","BF","nqtl",
             "npar","rss")
  nfixcov <- qb.get(qbObject, "nfixcov")

  out <- qb.commonone(qbObject, "sliceone", slice, epistasis, scan, type.scan, covar,
                      adjust.covar, chr, sum.scan, min.iter, aggregate,
                      smooth, weight, split.chr, center.type,, verbose, ...)
  ## Somehow we lose the slice column in this.
  class(out) <- c("qb.sliceone", class(out))
  out
}
###################################################################
summary.qb.sliceone <- function(object, chr = attr(object, "chr"), ...)
{
#  class(object) <- class(object)[-1]
#  summary(object, chr = chr, ...)
  NextMethod()
}
###################################################################
print.qb.sliceone <- function(x, ...) print(summary(x, ...))
###################################################################
plot.qb.sliceone <- function(x, ..., scan, auto.par = TRUE)
{
  if(attr(x, "type.scan") != "cellmean") {
## Need to fix this, but first fix slice col in qb.sliceone.
    if(missing(scan)) {
      scan <- names(x)[-(1:2)]
      scan <- scan[scan != "slice"]
    }
    if(any(scan == "slice")) {
      scan.name <- attr(x, "type.scan")
      attr(x, "type.scan") <- "slice"
      attr(x, "reference") <- mean(x$slice)
      return(invisible(plot.qb.scanone(x, ..., scan.name = scan.name,
                                       scan = scan)))
    }
    else
      return(invisible(plot.qb.scanone(x, ..., scan = scan)))
  }
  if(missing(scan)) {
    cellmean.plot(x, ..., scan = scan, auto.par = auto.par)
  }
  else {
    if(length(scan) == 1 & scan[1] == "slice")
      attr(x, "type.scan") <- paste("chr", attr(x, "slice")["chr"], "(cM)")
    plot.qb.scanone(x, scan = scan, ...)
  }
  invisible()
}
######################################################
cellmean.plot <- function(x, col = cols, ..., scan, auto.par = TRUE)
{
  is.bc <- (attr(x, "cross.class") == "bc")
  if(is.bc) {
    if(auto.par) {
      tmpar <- par(mfcol=c(2,1), mar = c(4.1,4.1,3.1,0.1))
      on.exit(par(tmpar))
    }
    cols <- c("blue","green")
    names(col) <- scans <- c("AA","HA")
    plot.qb.scanone(x, ..., scan = scans, col = col)
    names(col) <- scans <- c("AH","HH")
    plot.qb.scanone(x, ..., scan = scans, col = col)
  }
  else {
    if(auto.par) {
      tmpar <- par(mfcol=c(3,1), mar = c(4.1,4.1,3.1,0.1))
      on.exit(par(tmpar))
    }
    cols <- c("blue","green","red")
    names(col) <- scans <- c("AA","HA","BA")
    plot.qb.scanone(x, ..., scan = scans, col = col)
    names(col) <- scans <- c("AH","HH","BH")
    plot.qb.scanone(x, ..., scan = scans, col = col)
    names(col) <- scans <- c("AB","HB","BB")
    plot.qb.scanone(x, ..., scan = scans, col = col)
  }
}
######################################################
qb.slicetwo <- function(qbObject, chr, pos, type.scan = "2logBF", width = 10, ...)
{
  qb.exists(qbObject)
  if(is.null(qb.get(qbObject, "pairloci", ...)))
    stop("slicetwo not possible without epistasis")
  
  qb.name <- deparse(substitute(qbObject))
  chr <- qb.find.chr(qbObject, unlist(chr), sort.chr = FALSE)
  pos <- unlist(pos)
  names(chr) <- names(pos) <- NULL
  if(length(chr) == 1)
    chr <- rep(chr, 2)
  if(length(chr) != 2 | length(pos) != 2)
    stop("chr and pos must be of length 2")

  ## Below will break, but it is in direction desired.
  slice <- list()
  for(i in 1:2) {
    ## Slice of objective function.
    ii <- paste("slice", i, sep = "")
    slice[[ii]] <- qb.sliceone(qbObject, chr=chr[i],
                               slice=c(chr=chr[3-i],
                                 start = pos[3-i] - width,
                                 end = pos[3-i] + width),
                               sum.scan = "no",
                               type.scan = type.scan, scan = "epistasis", ...)
    
    ## Slice of Cockerham parameter estimates.
    tmp <- qb.sliceone(qbObject, chr=chr[i],
                       slice=c(chr=chr[3-i],
                         start = pos[3-i] - width,
                         end = pos[3-i] + width),
                       type.scan = "estimate", scan = "epistasis",
                       aggregate = FALSE, sum.scan = "no", ...)
    
    ## Slice of Cell means.
    tmp2 <-
      qb.sliceone(qbObject, chr=chr[i],
               slice=c(chr=chr[3-i],
                 start = pos[3-i] - width,
                 end = pos[3-i] + width),
               type.scan = "cellmean", sum.scan = "no", ...)

    vars <- names(tmp)[-(1:2)]
    vars <- vars[-length(vars)]
    cells <- names(tmp2)[-(1:2)]
    attrs <- attributes(slice[[ii]])
    slice[[ii]]$slice <- NULL
    names(slice[[ii]])[3] <- attr(slice[[ii]], "type.scan")
    for(i in vars)
      slice[[ii]][[i]] <- tmp[[i]]
    for(i in cells)
      slice[[ii]][[i]] <- tmp2[[i]]

    ## Set up metadata needed for summary, plot.
    attr(slice[[ii]], "type.scans") <- c(obj = attr(slice[[ii]], "type.scan"),
                                      est = "effects",
                                      mean = attr(tmp2, "type.scan"))
    attr(slice[[ii]], "scans") <- list(obj = attr(slice[[ii]], "scan"),
                                       est = vars, 
                                       mean = attr(tmp2, "scan"))
    attr(slice[[ii]], "scan") <- c(type.scan, vars, attr(tmp2, "scan"), "slice")
    attr(slice[[ii]], "refs") <- c(obj = attr(slice[[ii]], "reference"),
                                   est = attr(tmp, "reference"),
                                   mean = attr(tmp2, "reference"))
  }
  slice$cross.object <- qb.cross(qbObject, genoprob = FALSE)
  
  class(slice) <- c("qb.slicetwo", "list")
  attr(slice, "pheno.col") <- qb.get(qbObject, "pheno.col")
  attr(slice, "type.scan") <- type.scan
  attr(slice, "chr") <- chr
  attr(slice, "pos") <- pos
  attr(slice, "step") <- qb.get(qbObject, "step")
  slice
}
######################################################
summary.qb.slicetwo <- function(object, ...)
{
  rbind(summary(object$slice1),
        summary(object$slice2))
}
######################################################
print.qb.slicetwo <- function(x, ...) print(summary(x, ...))
######################################################
plot.qb.slicetwo <- function(x, byrow = TRUE,
                             figs = fig.options,
                             auto.par = TRUE,
                             col = cols, lty = 1,
                             ...)
{
  ## Just fixed this problem with summary.qb.scanone
  ## in terms of catching chr properly for slice.
  ## However, negative effects not showing properly.
  
  chr <- attr(x, "chr")
  pos <- attr(x, "pos")

  cross <- x$cross.object
  markers <- find.marker(cross, names(cross$geno)[chr], pos)

  is.bc <- class(cross)[1] == "bc"
  if(is.bc) {
    scans <- c("AA","AH","HA","HH")
    cols <- c("blue","purple","green","red")
    col <- array(col, 4)
    names(col) <- scans
  }
  else {
    scans <- list(A = c("AA","HA","BA"),
                  H = c("AH","HH","BH"),
                  B = c("AB","HB","BB"))
    cols <- c("blue","purple","red")
    col <- array(col, 3)
    names(cols) <- names(scans)
  }
  fig.options <- c("profile", "effects", "cellmean", "effectplot")
  figs <- match.arg(figs, fig.options, several.ok = TRUE)
#  if(!epistasis) {
#    tmp <- fig.options[3:4][match(figs, fig.options[3:4], nomatch = 0)]
#    if(length(tmp) < length(figs) & missing(epistasis))
#      warning("profile and effects not shown: no epistasis")
#    figs <- tmp
#  }
  epistasis <- TRUE
  is.profile <- any(pmatch(tolower(figs), "profile", nomatch = 0))
  is.effects <- any(pmatch(tolower(figs), "effects", nomatch = 0))
  is.cellmean <- any(pmatch(tolower(figs), "cellmean", nomatch = 0))
  is.effectplot <- any(pmatch(tolower(figs), "effectplot", nomatch = 0))
  if(is.effectplot & is.null(cross$geno[[1]]$draws))
    cross <- sim.geno(cross, step = attr(x, "step"))
  num.plots <- is.profile + is.effects + is.effectplot +
    (3 - 2 * (is.bc | !epistasis)) * is.cellmean

  if(auto.par) {
    if(byrow)
      tmpar <- par(mfrow=c(2, num.plots), mar=c(4.1,4.1,3.1,0.1))
    else
      tmpar <- par(mfcol=c(num.plots, 2), mar=c(4.1,4.1,3.1,0.1))
    on.exit(par(tmpar))
  }

  for(i in 1:2) {
    chrbychr <- paste("\nchr ", chr[i], "@", round(pos[i]),
                      " by ", chr[3-i], "@", round(pos[3-i]), sep = "")
    ii <- paste("slice", i, sep = "")
    slice.object <- x[[ii]]
    if(is.profile) {
      ## Slice of objective function.
      type.scan <- attr(slice.object, "type.scan") <- attr(x[[ii]], "type.scans")[1]
      attr(slice.object, "scan") <- attr(x[[ii]], "scans")[[1]]
      attr(slice.object, "reference") <- attr(x[[ii]], "refs")[1]
      plot(slice.object, scan = type.scan, main = paste("epistasis", type.scan, chrbychr))
      abline(v = pos[i], lty = 2, col = "red")
    }
    ## Add back in objects from objn.
    if(is.effects) {
      ## Slice of Cockerham parameter estimates.
      type.scan <- attr(slice.object, "type.scan") <- attr(x[[ii]], "type.scans")[2]
      attr(slice.object, "scan") <- attr(x[[ii]], "scans")[[2]]
      attr(slice.object, "reference") <- attr(x[[ii]], "refs")[2]
      plot(slice.object, scan = attr(slice.object, "scan"),
           main = paste("epistasis", type.scan, chrbychr))

      abline(v = pos[i], lty = 2, col = "red")
    }
    if(is.cellmean) {
      ## Slice of Cell means.
      type.scan <- attr(slice.object, "type.scan") <- attr(x[[ii]], "type.scans")[3]
      attr(slice.object, "scan") <- attr(x[[ii]], "scans")[[3]]
      attr(slice.object, "reference") <- attr(x[[ii]], "refs")[3]
      if(is.bc) {
        index <- if(epistasis)
          c(1,i+1,4-i,4)
        else
          c(1,4)
        ltys <- array(lty, 4)
        plot(slice.object, scan = scans[index],
             col = col[index], lty = ltys[index],
             main = paste(type.scan, chrbychr))
        abline(v = pos[i], lty = 2, col = "red")
      }
      else {
        tmp <- c("A","H","B")
        if(!epistasis)
          tmp <- "A"
        for(j in tmp) {
          if(epistasis)
            names(col) <- paste(tmp, j, sep = "")
          plot(slice.object, scan = scans[[j]], col = col,
               main = paste(type.scan,
                 ifelse(epistasis, paste("slice =", j), ""),
                 chrbychr))
          abline(v = pos[i], lty = 2, col = "red")
        }
      }
    }
    if(is.effectplot) {
      cols <- col
      if(is.bc)
        cols <- col[c(1, length(col))]
      effectplot(cross, attr(x, "pheno.col"), col = cols,
                 mname1 = markers[i], mname2 = markers[3-i],
                 main = paste("interaction plot", chrbychr))
    }
  }
  invisible()
}
