
#*******************************************************************************
# fitting main effects
fit.main <- function(object, x.main, prior.scale=2.5, prior.df=1, 
                     select=TRUE, vars.keep=NULL, p=0.05, one.time=50, ...)
{ 
  if(!"glm"%in%class(object)&!"polr"%in%class(object))
    stop("model not allowed yet") 
  start.time = Sys.time()
  f = object
  vars = model.matrix(f, na.action=na.pass)
  if(nrow(vars)!=nrow(x.main)) stop("missing data not allowed") 
  if(colnames(vars)[1]=="(Intercept)") vars = vars[, -1, drop=F]
  intercept = (names(coef(f))[1]=="(Intercept)")
  n.fit = if(intercept) length(coef(f))-1 else length(coef(f))
  if(is.null(vars.keep)) vars.keep = n.fit 
  else if(vars.keep > n.fit) vars.keep = n.fit 
  pre.main = f$main
  if(select){
    n.main = ncol(x.main)                                       
    n.time = as.integer(n.main/one.time)+1
    prior.scale = min(prior.scale)
    prior.df = min(prior.df)
    for(j in 1:n.time){
      start = if((j-1)*one.time + 1 > 0) (j-1)*one.time + 1 else 1
      end = if(j*one.time <= n.main) j*one.time else n.main
      if(start > end) break
      vars = data.frame(vars, x.main[,start:end,drop=F], check.names=F)
      vars = vars[,unique(colnames(vars)),drop=F]
      f = update(f, data=data.frame(vars), prior.scale=prior.scale, prior.df=prior.df)
      vars = GetEffects(f, vars.keep=vars.keep, p=p)
      if(j >= n.time) break
    }
    if(!is.null(vars)&ncol(vars)>0)
      f = update(f, data=data.frame(vars), prior.scale=2.5, prior.df=1)
    else stop("no variables detected")  
  }
  else {
    vars = data.frame(vars, x.main, check.names=F)
    vars = vars[,unique(colnames(vars)),drop=F]
    f = update(f, data=data.frame(vars), prior.scale=prior.scale, prior.df=prior.df)
  }
  if(vars.keep+1 <= ncol(vars))
    main = as.matrix(vars[, (vars.keep+1):ncol(vars), drop=F])
  else {
    main = NULL
    warning("no main effects detected")
  }
  f$main = cbind(pre.main, main)
  f$main = f$main[,unique(colnames(f$main)),drop=F] 
  stop.time = Sys.time()
  f$minutes = round(difftime(stop.time, start.time, units="min"),3)
  f
}

#*******************************************************************************
# fitting GxG or GxE      
fit.interaction <- function(object, x.main, inter.type = c("GxG","GxE"), 
                            x.inter=NULL, all.gg=FALSE, intcovar,
                            prior.scale=2.5, prior.df=1, 
                            select=TRUE, vars.keep=NULL, p=0.01, one.time=100, ...)
{
  if(!"glm"%in%class(object)&!"polr"%in%class(object))
    stop("model not allowed yet") 
  start.time = Sys.time()
  f = object
  vars = model.matrix(f, na.action=na.pass)
  if(nrow(vars)!=nrow(x.main)) stop("missing data not allowed") 
  if(colnames(vars)[1]=="(Intercept)") vars = vars[, -1, drop=F]
  intercept = (names(coef(f))[1]=="(Intercept)")
  n.fit = if(intercept) length(coef(f))-1 else length(coef(f))
  if(is.null(vars.keep)) vars.keep = n.fit 
  else if(vars.keep > n.fit) vars.keep = n.fit
  inter.type = inter.type[1]
  if(inter.type=="GxG"&is.null(x.inter)){
    if(all.gg) x.inter = make.inter(x.main, x.main, back=1)
    else x.inter = make.inter(f$main, x.main, back=1)
  }
  if(inter.type=="GxE"&is.null(x.inter)) 
    x.inter = make.inter(x.main, intcovar, back=0)
  n.inter = ncol(x.inter)
  n.main = ncol(x.main)
  pre.main = f$main
  if(select){
    n.time = as.integer(n.inter/one.time)+1
    prior.scale = min(prior.scale)
    prior.df = min(prior.df)
    s.inter = if(n.main >= n.inter) prior.scale else prior.scale*n.main/n.inter
    for(j in 1:n.time){
      start = if((j-1)*one.time + 1 > 0) (j-1)*one.time + 1 else 1
      end = if(j*one.time <= n.inter) j*one.time else n.inter
      if(start > end) break
      if(n.fit > 0) p.s = c(rep(prior.scale, n.fit),rep(s.inter, end-start+1))
      else p.s = rep(s.inter, end-start+1) 
      vars = data.frame(vars, x.inter[, start:end, drop=F], check.names=F)
      vars = vars[,unique(colnames(vars)),drop=F]
      f = update(f, data=data.frame(vars), prior.scale=p.s, prior.df=prior.df)
      vars = GetEffects(f, vars.keep=vars.keep, p=p)
      n.fit = if(!is.null(vars)) ncol(vars) else 0
      if(j >= n.time) break
    }
    if(!is.null(vars)&ncol(vars)>0)
      f = update(f, data=data.frame(vars), prior.scale=2.5, prior.df=1)
    else stop("no variables detected") 
  }
  else {
    if(length(prior.scale)==1){
      s.inter = if(n.main >= n.inter) prior.scale else prior.scale*n.main/n.inter
      prior.scale = c(rep(prior.scale, n.fit), rep(s.inter, n.inter))
    }
    vars = data.frame(vars, x.inter, check.names=F)
    vars = vars[,unique(colnames(vars)),drop=F]
    f = update(f, data=data.frame(vars), prior.scale=prior.scale, prior.df=prior.df)
  }
  f$main = pre.main
  stop.time = Sys.time()
  f$minutes = round(difftime(stop.time, start.time, units="min"),3)
  f
}
#*******************************************************************************
# get significant effects from a fitted model
# vars.keep is the first vars.keep variables (don't include intercept) to be kept
GetEffects <- function(object, vars.keep=0, p=0.05, sig=TRUE, na.action=na.pass)
{
  if(!"glm"%in%class(object)&!"polr"%in%class(object))
    stop("object not allowed yet")
  
  if(any(class(object)=="glm")) {
    x = model.matrix(object, na.action=na.action)
    colnames(x) = names(object$coefficients)
    pvalue = summary(object, dispersion = object$dispersion)$coefficients[,4]
    if(colnames(x)[[1]]=="(Intercept)"){
      x = x[,-1,drop=F]
      pvalue = pvalue[-1]
    }
  }

  if(any(class(object)=="polr")) {
    x = model.matrix(object)[,-1,drop=F]
    tvalue = summary(object)$coefficients[,3]
    df.r = object$df.residual
    if(df.r > 0) pvalue = 2 * pt(-abs(tvalue), df.r)
    else pvalue = 2 * pnorm(-abs(tvalue))
    pvalue = pvalue[1:ncol(x)]
  }
                                              
  if((vars.keep+1)<=ncol(x)&any(pvalue[(vars.keep+1):ncol(x)]>p)) {
    vars.delete = vars.keep + which(pvalue[(vars.keep+1):ncol(x)]>p)
    if(sig) x.new = x[,-vars.delete,drop=F]
    else x.new = x[,vars.delete,drop=F]
  }
  else {
    if(sig) x.new = x
    else x.new = NULL
  }
  x.new
}
#*******************************************************************************
# construct main-effect variable matrix
make.main <- function (cross, loci.names = c("marker","position"), 
          model = c("Cockerham","codominant","additive","dominant","recessive","overdominant","full"),
          imprint = TRUE, multipoint = TRUE, cc.col = NULL, fill.missing = TRUE, geno.order = TRUE, 
          effect.order = TRUE, aggregate = TRUE, pos.digits = 1, ...)
{
  if(class(cross)[2]!="cross") stop("not an object of class 'cross'")
  loci.names = loci.names[1]
  model = model[1]
  if(!any(loci.names == c("marker","position"))) {
    loci.names = c("marker")
    warning("You give wrong name. use marker names.") 
  }
  if(!any(model == c("Cockerham","codominant","additive","dominant","recessive","overdominant","full"))){
    model == c("Cockerham")
    warning("You give wrong model. use Cockerham.") 
  }
  
  if(!multipoint) { #remove markers with less than 3 (for F2) or 2 (for BC) genotypes
    for(i in 1:nchr(cross)) {
      geno.m = apply(cross$geno[[i]]$data,2,max,na.rm=T)
      if(class(cross)[1] == "f2") m = 3
      else m = 2 
      if(any(geno.m < m)) {
        geno.0 = geno.m[geno.m < m]
        cross = drop.markers(cross,names(geno.0))
        warning("markers with incomplete genotypes removed")
      }
    }
  }
  if(!multipoint&geno.order) { # order genotypes: AA common, AB, BB rare
    func = function(marker) {
      x = table(marker)
      marker = (x[1]<x[3])*(4-marker) + (x[1]>=x[3])*marker
    }
    for(i in 1:nchr(cross))
      cross$geno[[i]]$data = apply(cross$geno[[i]]$data,2,func)
  } 
  
  n.chr = nchr(cross)
  n.ind = nind(cross)
  n.mar = nmar(cross)
  x = vector(mode="list",length=n.chr)
  names(x) = names(cross$geno)

  if(multipoint){  
    if( is.null(cross$geno[[1]]$prob) ) {
      warning("first run qb.genoprob(cross, step=0)") 
      cross = qb.genoprob(cross, step=0) 
    }
    else if(attr(cross$geno[[1]]$prob,"step")!=0) cross = qb.genoprob(cross, step=0)  
    grid = pull.loci(cross) # in order to deal with loci in marker intervals
  }
  
  if(!multipoint) {
    grid = cross$geno
    for(i in 1:n.chr) {
      grid[[i]] = cross$geno[[i]]$map
      n.gen = max(cross$geno[[i]]$data[,1],na.rm=T) # number of genotypes
      cross$geno[[i]]$prob = array(1, c(n.ind,n.mar[i],n.gen))
      for(k in 1:n.gen) cross$geno[[i]]$prob[,,k] = (cross$geno[[i]]$data==k)*1
      colnames(cross$geno[[i]]$prob) = names(cross$geno[[i]]$map)
      if(n.gen==3) dimnames(cross$geno[[i]]$prob)[[3]] = c("AA","AB","BB")
      else {
        if(class(cross$geno[[i]])=="A") dimnames(cross$geno[[i]]$prob)[[3]] = c("AA","AB")
        else dimnames(cross$geno[[i]]$prob)[[3]] = c("g1","g2")
      }
      attr(cross$geno[[i]]$prob,"map") = cross$geno[[i]]$map
      attr(cross$geno[[i]]$prob,"error.prob") = 1e-04
      attr(cross$geno[[i]]$prob,"step") = 0
      attr(cross$geno[[i]]$prob,"off.end") = 0
      attr(cross$geno[[i]]$prob,"map.function") = "haldane"
      attr(cross$geno[[i]]$prob,"stepwidth") = "variable"
      attr(cross$geno[[i]]$prob,"tolerance") = 1e-06
      if(fill.missing) {
        if(is.null(cc.col)|length(unique(cross$pheno[,cc.col]))!=2) {
          w = apply(cross$geno[[i]]$data,2,table)
          wsum = apply(w,2,sum)
          for(k in 1:n.gen) {
            z = rep(w[k,]/wsum,each=n.ind)
            na.index = which(is.na(cross$geno[[i]]$prob[,,k]))
            cross$geno[[i]]$prob[,,k][na.index] = z[na.index] 
          }
        }
        if(!is.null(cc.col)&length(unique(cross$pheno[,cc.col]))==2) {
          cc = cross$pheno[,cc.col]
          cc.type = unique(cc)
          cross1 = subset(cross, ind=(cc==cc.type[1]))
          w1 = apply(cross1$geno[[i]]$data,2,table)
          w1sum = apply(w1,2,sum)
          cross2 = subset(cross, ind=(cc==cc.type[2]))
          w2 = apply(cross2$geno[[i]]$data,2,table)
          w2sum = apply(w2,2,sum)
          for(k in 1:n.gen) {
            z1 = rep(w1[k,]/w1sum,each=n.ind)
            z2 = rep(w2[k,]/w2sum,each=n.ind)
            ccrep = rep(cc,n.mar[i])
            na.index1 = which(is.na(cross$geno[[i]]$prob[,,k])&ccrep==cc.type[1])
            na.index2 = which(is.na(cross$geno[[i]]$prob[,,k])&ccrep==cc.type[2])
            cross$geno[[i]]$prob[,,k][na.index1] = z1[na.index1] 
            cross$geno[[i]]$prob[,,k][na.index2] = z2[na.index2] 
          }
        }
      }
    }
  }
   
  for(i in 1:n.chr)
  {
    if(class(grid[[i]])=="matrix") grid[[i]] = grid[[i]][1,]
    grid[[i]] = round(grid[[i]],digits=pos.digits)
    if(loci.names == "marker")   
      names(grid[[i]]) = paste(names(grid[[i]]))  
    if(loci.names == "position"){ 
      names(grid[[i]]) = paste("c",names(grid)[i],".",grid[[i]],sep="")
    #  if(!dataframe&old.name) 
    #    names(grid[[i]]) = paste("(",names(grid)[i],"@",grid[[i]],")",sep="") 
    }           
    if(dim(cross$geno[[i]]$prob)[3]==2|class(cross$geno[[i]])=="X"){
      x[[i]] = as.matrix(0.5*(cross$geno[[i]]$prob[,,2]-cross$geno[[i]]$prob[,,1]))
      colnames(x[[i]]) = names(grid[[i]])
    }
      
    if(dim(cross$geno[[i]]$prob)[3]==3&class(cross$geno[[i]])!="X"){
      if (model == "codominant") {
        x1 = as.matrix(cross$geno[[i]]$prob[,,3])
        x2 = as.matrix(cross$geno[[i]]$prob[,,2])
        colnames(x1) = paste(names(grid[[i]]),"a",sep="")
        colnames(x2) = paste(names(grid[[i]]),"d",sep="")
        x[[i]] = cbind(x1,x2)
      }
      if (model == "Cockerham") {
        x1 = as.matrix(cross$geno[[i]]$prob[,,3]-cross$geno[[i]]$prob[,,1]); x1 = x1/sqrt(2)
        x2 = as.matrix(0.5*(cross$geno[[i]]$prob[,,2]-cross$geno[[i]]$prob[,,1]-cross$geno[[i]]$prob[,,3]))
        colnames(x1) = paste(names(grid[[i]]),"a",sep="")
        colnames(x2) = paste(names(grid[[i]]),"d",sep="")
        x[[i]] = cbind(x1,x2)
      }
      if (model == "additive") {
        x[[i]] = as.matrix(cross$geno[[i]]$prob[,,3]-cross$geno[[i]]$prob[,,1])
        colnames(x[[i]]) = names(grid[[i]])
      }
      if (model == "dominant") {
        x[[i]] = as.matrix(cross$geno[[i]]$prob[,,1]+cross$geno[[i]]$prob[,,2])
        colnames(x[[i]]) = names(grid[[i]])
      }
      if (model == "recessive") {
        x[[i]] = as.matrix(cross$geno[[i]]$prob[,,3]+cross$geno[[i]]$prob[,,2])
        colnames(x[[i]]) = names(grid[[i]])
      }
      if (model == "overdominant") {
        x[[i]] = as.matrix(cross$geno[[i]]$prob[,,2])
        colnames(x[[i]]) = names(grid[[i]])  
      }
      if (model == "full") {
        x1 = as.matrix(cross$geno[[i]]$prob[,,1])
        x2 = as.matrix(cross$geno[[i]]$prob[,,2])
        x3 = as.matrix(cross$geno[[i]]$prob[,,3])
        colnames(x1) = paste(names(grid[[i]]),"1",sep="")
        colnames(x2) = paste(names(grid[[i]]),"2",sep="")
        colnames(x3) = paste(names(grid[[i]]),"3",sep="")
        x[[i]] = cbind(x1,x2,x3)
      }
    }
    
    if(dim(cross$geno[[i]]$prob)[3]==4&class(cross$geno[[i]])!="X"){
      if (model != "full") {
        x1 = as.matrix(cross$geno[[i]]$prob[,,4]-cross$geno[[i]]$prob[,,1])
        x2 = as.matrix(0.5*(cross$geno[[i]]$prob[,,2]+cross$geno[[i]]$prob[,,3]-cross$geno[[i]]$prob[,,1]-cross$geno[[i]]$prob[,,4]))
        x3 = as.matrix(0.5*(cross$geno[[i]]$prob[,,3]-cross$geno[[i]]$prob[,,2]))
        colnames(x1) = paste(names(grid[[i]]),"a",sep="")
        colnames(x2) = paste(names(grid[[i]]),"d",sep="")
        colnames(x3) = paste(names(grid[[i]]),"i",sep="")
        x[[i]] = cbind(x1,x2)
        if (imprint) x[[i]] = cbind(x[[i]],x3)
      }
      if (model == "full") {
        x1 = as.matrix(cross$geno[[i]]$prob[,,1])
        x2 = as.matrix(cross$geno[[i]]$prob[,,2])
        x3 = as.matrix(cross$geno[[i]]$prob[,,3])
        x4 = as.matrix(cross$geno[[i]]$prob[,,4])
        colnames(x1) = paste(names(grid[[i]]),"1",sep="")
        colnames(x2) = paste(names(grid[[i]]),"2",sep="")
        colnames(x3) = paste(names(grid[[i]]),"3",sep="")
        colnames(x4) = paste(names(grid[[i]]),"4",sep="")
        x[[i]] = cbind(x1,x2,x3,x4)
      }
    }
    
    if (ncol(x[[i]])>n.mar[i]&effect.order) {
      o = order(rep(1:n.mar[i],as.integer(ncol(x[[i]])/n.mar[i])))
      x[[i]] = x[[i]][,o]
    }
    x[[i]] = data.frame(x[[i]][,unique(colnames(x[[i]])),drop=F]) 
  }
  if(aggregate){
    x.main = x[[1]]
    if(nchr(cross)>1) 
      for(j in 2:nchr(cross)) x.main = data.frame(x.main, x[[j]])
    x.main = x.main[,unique(colnames(x.main)),drop=F] 
  }
  else x.main = x
  cross$x.main = x.main
  cross$model = model
  cross
}
#*******************************************************************************
# construct interaction matrix
make.inter <- function(x1, x2, back = 1, na.action = na.pass)  
{
#back = 0, 1
  if(class(x1)=="factor"){
    f = ~ x1
    mf = model.frame(f, na.action=na.action)
    x1 = model.matrix(f, mf)[,-1,drop=F]
  }
  if(class(x2)=="factor"){
    f = ~ x2
    mf = model.frame(f, na.action=na.action)
    x2 = model.matrix(f, mf)[,-1,drop=F]
  }
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  if( is.null(colnames(x1))|is.null(colnames(x2)) ) 
    stop("no colnames for the inputs")
  f = ~ x1:x2 - 1
  mf = model.frame(f, na.action=na.action)
  z = model.matrix(f, mf)
  if(ncol(x1)==1|ncol(x2)==1) colnames(z) = paste( paste("x1",colnames(x1),sep=""),":",
                                                   paste("x2",colnames(x2),sep=""), sep="" )
  znames = strsplit(colnames(z),split=":",fixed=T)
  znames0 = lapply(znames, function(x) substr(x,3,nchar(x)-back))
  index = unlist(lapply(znames0, function(x) (x[1]!=x[2])))
  z1 = z[,index,drop=F]
  z1names = strsplit(colnames(z1),split=":",fixed=T)
  z1names0 = lapply(z1names, function(x) substr(x,3,nchar(x))) 
  z1names1 = lapply(z1names0, function(x) sort(c(x[1],x[2])))
  z1names2 = unique(z1names1)
  z2 = z1[,match(z1names2,z1names1),drop=F] 
  colnames(z2) = lapply(z1names2, function(x) paste(x[1],":",x[2],sep=""))
  z2
}
#*******************************************************************************
# clean data for matched case-control
clean.cc <- function(cross, cc.col, stra.col)
{
  stra = cross$pheno[,stra.col]
  cc = cross$pheno[,cc.col]
  z = table(stra,cc)
  w = z[,1]*z[,2]
  stra.0 = names(w[w!=0])
  ind = which(stra%in%stra.0)
  cross.new = subset(cross, ind=ind)
  cross.new
}
#*******************************************************************************
# select 1:1 or 1:m with the fewest alleles
select.cc <- function(cross, cc.col, stra.col, design = c("1:m","1:1"))
{
  design = design[1]
  if(!design%in%c("1:m","1:1")) stop("design not allowed")
  stra = cross$pheno[,stra.col]
  cc = cross$pheno[,cc.col]
  vars = cross$geno[[1]]$data
  if(nchr(cross)>1)
    for(j in 2:nchr(cross)) vars = cbind(vars, cross$geno[[j]]$data)
  data = as.data.frame(cbind(c(1:nind(cross)), cc, vars))
  by.stra = split(data, stra)

  pairs = function(x) {
    cc = x[,2]
    control = x[cc==0,]
    case = x[cc==1,]
    n.control = nrow(control)
    n.case = nrow(case)
    id = c(control[,1],case[,1])
    if(n.control>0&n.case>0){
      control.var = control[,3:ncol(x)]
      case.var = case[,3:ncol(x)]
      if(design=="1:1"){
        dif = 0; j0 = 1; h0 = 1
        if(n.control>1|n.case>1){
          for(j in 1:n.control)
            for(h in 1:n.case) {
              c0 = as.numeric(control.var[j,])
              c1 = as.numeric(case.var[h,])
              dif0 = sum(abs(c1-c0),na.rm = T)
              if(dif0 > dif) {
                j0  = j; h0 = h
                dif = dif0
              }
          }
        }
        id = c(control[j0,1],case[h0,1])
      }
      if(design=="1:m"){
        dif = 0; h0 = 1
        if(n.case>1){
          for(h in 1:n.case) {
            c0 = as.matrix(control.var)
            c1 = as.numeric(case.var[h,])
            dif0 = 0
            for(j in 1:n.control) dif0 = dif0 + sum(abs(c1-c0[j,]),na.rm = T)
            if(dif0 > dif) {
              h0 = h
              dif = dif0
            }
          }
        }
        id = c(control[,1],case[h0,1])
      }
    }
    id
  }

  id = unlist(lapply(by.stra, pairs))
  cross = subset(cross, ind=id)
  cross
}
#*******************************************************************************
# geno is a matrix of all genotypes 
geno.freq <- function(geno, freq=TRUE)
{ 
  counts = t(apply(geno, 2, table))
  counts[,c(1,3)] = t(apply(counts[,c(1,3)], 1, sort, decreasing=T))
  colnames(counts) = c("c", "h", "r")
  p.HWE = apply(counts, 1, function(x){
                    pa = (x[1] + x[2]/2)/sum(x) 
                    chisq.test(x, p = c(pa^2, 2*pa*(1-pa), (1-pa)^2))$p.value }
               )
  if(freq) counts = t(apply(counts, 1, function(x) x/sum(x)))
  missing = apply(geno, 2, function(x) length(x[is.na(x)]))
  results = list()
  results$genotypes = c("c: common homozygote", "h: heterzygote", "r: rare homozygote")
  results$freq = cbind(missing, counts, p.HWE)
  results
}
#*******************************************************************************
ROC <- function(dn, pre, m = 100, add = FALSE, lwd = 2, lty = 1, col = "black")               
{
  # draw roc curve
  n <- length(dn)                   # sample size 
  d <- length(dn[dn != 0])          # number of disease (dn=1)
  sr <- (max(pre) - min(pre)) / m   # size of reclassifying interval
  tpf <- rep(NA, m)                 # true positive fraction
  fpf <- rep(NA, m)                 # false positive fraction
  for (i in 1:(m + 1)){
    tpf[i] <- sum(pre >= max(pre) - sr*(i - 1) & dn == 1) / d
    fpf[i] <- sum(pre >= max(pre) - sr*(i - 1) & dn == 0) / (n-d)
  }  
  if(!add) plot(c(0, 1), c(0, 1), xlim=c(0, 1), ylim=c(0, 1), xlab = "False Positive Rate",
                ylab = "True Positive Rate", type = "l",lty=1, col = "gray")    # draw a diagonal line
  lines(fpf, tpf, type = "l", lwd = lwd, lty = lty, col = col)

  # calculate area under curve based on non-parametric statistic (Mann-Whitney U)    
  tp2 <- rep(0, m)                  # number of true positive for each interval of cut-off values
  fp2 <- rep(0, m)                  # number of false positive for each interval of cut-off values
  for (i in 1:m){
    tp2[i] <- sum((pre >= min(pre) + sr*(i - 1) & pre < min(pre) + sr * (i)) & dn == 1)
    fp2[i] <- sum((pre >= min(pre) + sr*(i - 1) & pre < min(pre) + sr * (i)) & dn == 0)      
  }

  dtp2 <- rep(0, m)                 # degression of tp2
  afp2 <- rep(0, m)                 # acumulation of fp2
  dtp2[1] <- (d - tp2[1])
  for (i in 1:(m - 1)){
    dtp2[i + 1] <- dtp2[i] - tp2[i + 1]
    afp2[i + 1] <- afp2[i] + fp2[i]
  }

  a <- rep(0, m)
  for (i in 1:m) {a[i] <- fp2[i] * dtp2[i] + (fp2[i] * tp2[i]) / 2}
  auc <- sum(a) / (d * (n - d))     # area under curve                 

  # test hypothesis of auc=0.5 (Mann-Whitney U test)
  q1 <- rep(0, m)                   # degression (leijian) of tp2
  q2 <- rep(0, m)                   # acumulation (leijia) of fp2
  for (i in 1:m){
    q1[i] <- fp2[i] * ((dtp2[i])^2 + tp2[i] * dtp2[i] + ((tp2[i])^2) / 3)
    q2[i] <- tp2[i] * ((afp2[i])^2 + fp2[i] * afp2[i] + ((fp2[i])^2) / 3)
  }
  sq1 <- sum(q1) / (d^2 * (n - d))
  sq2 <- sum(q2) / (d * (n - d)^2)
  se <- sqrt((auc * (1 - auc) + (d - 1) * (sq1 - auc^2) + ((n - d) - 1) * (sq2 - auc^2)) / (d * (n - d)))
  z <- (auc - 0.5) / se             # z score
  p <- 2 * pnorm(-abs(z))           # p-value for two-tailed test
  lci <- auc - 1.96 * se            # left boundary of 95% confidence interval
  rci <- auc + 1.96 * se            # right boundary of 95% confidence interval
#  w <- matrix( c(auc, z, p, lci, rci), 5, 1, dimnames = list(c("Area under curve (AUC)",
#             "Mann-Whitney U statistic (Z)", "p-value", "95% confidence interval", ""), c("")) )     
  w = list(TPF=tpf,FPF=fpf,AUC=auc,Test.statistics=z,P.value=p,CI=c(lci,rci))
  w         
}
#*******************************************************************************
# Average predictive comparison
geno.effects <- function(object, x.covar = NULL, x.main) 
{
  x.fit = model.matrix(object)
  beta = object$coef
  eta = as.numeric(x.fit%*%beta)  # linear predictor
  geno.mean = mean(exp(eta)/(1+exp(eta))) # population mean
  
  if(!is.null(x.covar)){
    x.covar = as.matrix(x.covar)
    if(nrow(x.covar)>nrow(x.fit)) x.covar = x.covar[as.numeric(rownames(x.fit)),]
    covar = colnames(x.covar)  # all covariates
    beta.covar = beta[covar]
    beta.covar[is.na(beta.covar)] = 0
    names(beta.covar) = covar
  }
  x.main = as.matrix(x.main)
  if(nrow(x.main)>nrow(x.fit)) x.main = x.main[as.numeric(rownames(x.fit)),]
  main = colnames(x.main)  # all main effects
  beta.main = beta[main]
  beta.main[is.na(beta.main)] = 0
  names(beta.main) = main
  x.fiti = x.fit[,!colnames(x.fit)%in%c(covar,main,"(Intercept)"),drop=F]
  beta.i = NULL  # interactions
  if(dim(x.fiti)[2]!=0) beta.i = beta[colnames(x.fiti)]
  
  if(max(x.main[,1])==1/2^0.5) geno = rbind(c(-1/2^0.5,-0.5),c(0,0.5),c(1/2^0.5,-0.5))
  else geno = rbind(c(0,0),c(0,1),c(1,0)) # only deal with Cockerham and Codominant
  dimnames(geno) = list(c("c","h","r"),c("a","d"))
  marker = unique(substr(main,1,nchar(main)-1))
  n.mar = length(marker) 
  
  prob1 = array(NA,c(n.mar,3)) # average genotypic effects for each marker
  dimnames(prob1) = list(marker,c("c","h","r"))
  for(j in 1:n.mar){
    m.col = c(2*j-1,2*j)
    mOld = as.numeric(x.main[,m.col]%*%beta.main[m.col])
    mNames = main[m.col] 
    for(k in 1:3) {
      x.mNew = geno[k,,drop=F]
      mNew = as.numeric(x.mNew%*%beta.main[m.col])
      eta.new = eta - mOld + mNew
      if(!is.null(beta.i)){
        for(h in 1:length(beta.i)){
          inter = unlist(strsplit(names(beta.i[h]),split=".",fixed=T))
          if(length(inter)>2) stop("wrong variable names: with .")
          x = cbind(x.covar,x.main)[,inter,drop=F] 
          if(mNames[1]%in%inter){
            iOld = x.fiti[,h]*beta.i[h]
            z = x[,mNames[1]!=inter]*geno[k,1]
            iNew = z*beta.i[h]
            eta.new = eta.new - iOld + iNew 
          } 
          if(mNames[2]%in%inter){
            iOld = x.fiti[,h]*beta.i[h]
            z = x[,mNames[2]!=inter]*geno[k,2]
            iNew = z*beta.i[h]
            eta.new = eta.new - iOld + iNew 
          }
        }
      }
      prob1[j,k] = mean(exp(eta.new)/(1+exp(eta.new)))
    }
  }
  results = list()
  results$geno.mean = geno.mean
  results$genotypes = c("c: common homozygote", "h: heterzygote", "r: rare homozygote")
  results$one = prob1
  
  if(!is.null(beta.i)){
    inter = list()
    for(i in 1:length(beta.i)) {
      inter[[i]] = sort(unlist(strsplit(names(beta.i[i]),split=".",fixed=T)))
      for(k in 1:2)
        if(!inter[[i]][k]%in%covar) 
          inter[[i]][k] = substr(inter[[i]][k],1,nchar(inter[[i]][k])-1) 
    }
    inter = unique(inter)
    prob2 = list()
    for(i in 1:length(inter)) {
      m.col1 = c(paste(inter[[i]][1],c("a","d"),sep=""))
      if(inter[[i]][1]%in%covar) m.col1 = inter[[i]][1]
      m.col2 = c(paste(inter[[i]][2],c("a","d"),sep=""))
      if(inter[[i]][2]%in%covar) m.col2 = inter[[i]][2]
      mOld1 = as.numeric(cbind(x.covar,x.main)[,m.col1,drop=F]%*%c(beta.covar,beta.main)[m.col1])
      mOld2 = as.numeric(cbind(x.covar,x.main)[,m.col2,drop=F]%*%c(beta.covar,beta.main)[m.col2])
      m1 =3; n1 = c("c","h","r")
      if(inter[[i]][1]%in%covar) {
        Cov = sort(unique(x.covar[,inter[[i]][1]]))
        m1 = length(Cov)
        n1 = as.character(Cov)
      }
      m2 =3; n2 = c("c","h","r") 
      if(inter[[i]][2]%in%covar) {
        Cov = sort(unique(x.covar[,inter[[i]][2]]))
        m2 = length(Cov)
        n2 = as.character(Cov)
      } 
      prob2[[i]] = array(NA,c(m1,m2))
      names(prob2)[i] = paste(inter[[i]][1],inter[[i]][2],sep=" ")
      dimnames(prob2[[i]]) = list(n1,n2) 
      for(k1 in 1:m1) 
        for(k2 in 1:m2){
          x.mNew1 = if(inter[[i]][1]%in%covar) Cov[k1] else geno[k1,,drop=F] 
          mNew1 = as.numeric(x.mNew1%*%c(beta.covar,beta.main)[m.col1])
          x.mNew2 = if(inter[[i]][2]%in%covar) Cov[k2] else geno[k2,,drop=F]
          mNew2 = as.numeric(x.mNew2%*%c(beta.covar,beta.main)[m.col2])
          eta.new = eta - mOld1 + mNew1 - mOld2 + mNew2
          for(h in 1:length(beta.i)){
            iOld = x.fiti[,h]*beta.i[h]
            inter0 = sort(unlist(strsplit(names(beta.i[h]),split=".",fixed=T)))
            x = cbind(x.covar,x.main)[,inter0,drop=F]
            z.old = x[,!inter0%in%c(m.col1,m.col2),drop=F]
            z.new = cbind(geno[k1,1],geno[k1,2],geno[k2,1],geno[k2,2])
            if(inter[[i]][1]%in%covar) z.new = cbind(Cov[k1],geno[k2,1],geno[k2,2])
            if(inter[[i]][2]%in%covar) z.new = cbind(geno[k1,1],geno[k1,2],Cov[k2])  
            colnames(z.new) = c(m.col1,m.col2) 
            if(ncol(z.old)==2) z = z.old
            if(ncol(z.old)==0) z = z.new[,inter0,drop=F]
            if(ncol(z.old)==1){
              z.new = z.new[,inter0[!inter0%in%colnames(z.old)]]
              z.new = rep(z.new,nrow(z.old))
              z = cbind(z.old,z.new)
            }
            iNew = z[,1]*z[,2]*beta.i[h]
            eta.new = eta.new - iOld + iNew 
          }
          prob2[[i]][k1,k2] = mean(exp(eta.new)/(1+exp(eta.new)))
        }
    }
    results$two = prob2 
  }
  
  class(results) = "Average genotypic effects"
  return(results)
}
#*******************************************************************************
bdrop <- function(object, vars.keep=0, p=0.001, prior.scale=2.5, prior.df=1, trace=TRUE)
{
  for(j in 1:100){
    if(!is.null(GetEffects(object, vars.keep=vars.keep, p=p, sig=F))){
      x = GetEffects(object, vars.keep=vars.keep, p=p)
      object = update(object, data=data.frame(x), prior.scale=prior.scale, prior.df=prior.df)
      if(trace) {
        cat("times =", j, "\n")
        print(summary(object, dispersion=object$dispersion)$coef, digits=3)
        cat("\n")
      }
     }
     else break
  }  
  object
}
#*******************************************************************************
summary.bglm <- function (object, ...) 
{
    out <- summary(object, dispersion = object$dispersion)
    coef.table <- out$coefficients
    herit <- coef.table[,1]^2 * apply(model.matrix(object), 2, var)/
             var(object$linear.pred + object$resid)
    herit <- as.matrix(herit)
    colnames(herit) <- "Heritability"
    OR <- exp(coef.table[,1])
    lower.95 <- exp(coef.table[,1] - 2*coef.table[,2])
    upper.95 <- exp(coef.table[,1] + 2*coef.table[,2])
    OR <- cbind(OR, lower.95, upper.95) 
    colnames(OR) <- c("OR", "lower .95", "upper .95") 
    out$herit <- herit
    out$OR <- OR
    class(out) <- c("summary.glm", "summary.bglm")
    return(out)
}
#*******************************************************************************
summary.bpolr <- function (object, ...) 
{
  out <- summary(object)
  tvalue <- out$coefficients[,3]
  df.r <- object$df.residual
  if(df.r > 0) pvalue <- 2*pt(-abs(tvalue), df.r)
  else pvalue <- 2*pnorm(-abs(tvalue))
  pvalue <- as.matrix(pvalue)
  colnames(pvalue) <- "Pr(>|t|)"
  out$coefficients <- cbind(out$coefficients, pvalue)
  class(out) <- c("summary.polr", "summary.bpolr")
  return(out)
}


























