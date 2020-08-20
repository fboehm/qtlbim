qb.outbredF2 <- function(dir="",phefile,genfile,mapfile,filestem="outbredF2") {
	#convert outbred F2 to inbred F2 allowing for four alleles per locus
	#read in data using QTL Express format

	if (missing(mapfile)) stop("The mapfile argument is missing.\n")
	if (missing(genfile)) stop("The genfile argument is missing.\n")
	if (missing(phefile)) stop("The phefile argument is missing.\n")

	#load and process mapfile
	if (!missing(dir) && dir!="") setwd(dir)
	file <- mapfile
	n.chr <- as.integer(scan(file,nlines=1,what="integer",quiet=TRUE))
	cat(n.chr,"\n",file="map.txt",append=TRUE)
	chr.name <- vector(mode="character",length=n.chr)
	marker.map <- vector(mode="list",length=n.chr)
	position <- vector(mode="list",length=n.chr)
	for (n in 1:n.chr) {
		tmp <- scan(file,nlines=1,skip=2*n,what="character",quiet=TRUE)
		n.mark.chr <- as.integer(tmp[2])
		chr.name[n] <- tmp[1]
		tmp2 <- scan(file,nlines=1,skip=(2*n)+1,what="character",quiet=TRUE)
		marker.index <- seq(from=1,to=length(tmp2),by=2)
		marker.map[[n]] <- tmp2[marker.index]
		position.index <- seq(from=2,to=length(tmp2),by=2)
		intervals <- as.integer(tmp2[position.index])
		#convert cM distances between markers into recombination fractions using the Kosambi mapping function
		position[[n]] <- 0.5*(exp(4*intervals/100)-1)/(exp(4*intervals/100)+1)
		cat(n.mark.chr,position[[n]],"\n",file="map.txt",sep=" ",append=TRUE)
		write(marker.map[[n]],file="map.txt",append=TRUE,ncolumns=1)
	}

	#load and process genfile
	file <- genfile
	n.markers <- scan(file,nlines=1,what="integer",quiet=TRUE)
	stopifnot(length(unlist(marker.map))==as.integer(n.markers))
	marker.gen <- scan(file,nlines=1,skip=1,what="character",quiet=TRUE)
	pop.code <- scan(file,nlines=1,skip=2,what="character",quiet=TRUE)
	na.strings <- c(scan(file,nlines=1,skip=4,what="character",quiet=TRUE),"-","NA")
	gen <- read.table(file,skip=5,header=FALSE,colClasses="character",na.strings=na.strings)
	dat1 <- gen[(gen$V5==pop.code[1] | gen$V5==pop.code[2]),-(2:4)]
	dat2 <- gen[gen$V5==pop.code[3],-(2:5)]
	dat3 <- gen[gen$V5==pop.code[4],-(2:5)]
	dat4 <- gen[,1:3]
	F0inds <- nrow(dat1)
	F1inds <- nrow(dat2)
	F2inds <- nrow(dat3)
	print("Recoding...",quote=FALSE)
	genotype <- .C("outbredF2_4way",as.character(file),as.integer(F0inds),as.integer(F1inds),as.integer(F2inds),PACKAGE="qtlbim")
	genotype <- read.table("newgenfile.txt",header=FALSE)
	print("The recoded genotypes are:",quote=FALSE)
	cat("AA is ",pop.code[1],"/",pop.code[1],"\n",sep="")
	cat("AB is ",pop.code[1],"/",pop.code[2],"\n",sep="")
	cat("BB is ",pop.code[2],"/",pop.code[2],"\n",sep="")
	cat("not BB is ",pop.code[1],"/?\n",sep="")
	cat("not AA is ",pop.code[2],"/?\n",sep="")
	cat("\n")
	#make the order of markers in the genotype matrix match the order of markers in the map file
	r <- unlist(marker.map)
	q <- order(r)
	#order the markers in the genotype matrix
	genotype <- genotype[,order(marker.gen)]
	sorted_genotype <- matrix(data=NA,nrow=nrow(genotype),ncol=ncol(genotype))
	for (i in 1:ncol(genotype)) sorted_genotype[,q[i]] <- genotype[,i]
	#order the individual IDs in the genotype file
	sorted_genotype <- cbind(dat3[,1],sorted_genotype)
	o <- order(sorted_genotype[,1])
	cleangeno <- sorted_genotype[o,]
	
	#load and process phefile
	file <- phefile
	pheno.names <- scan(file,nlines=1,skip=1,what="character",quiet=TRUE)
	na.strings <- c(scan(file,nlines=1,skip=2,what="character",quiet=TRUE),"-","NA")
	pheno <- read.table(file,skip=3,header=FALSE,na.strings=na.strings)
	#order the individuals in the phenotype file
	p <- order(pheno[,1])
	#determine which individuals are genotyped but not phenotyped
	qqq <- setdiff(as.character(sorted_genotype[o,1]),as.character(pheno[p,1]))
	if (length(qqq)!=0) {
		qqq1 <- sapply(1:length(qqq),function(a) {which(sorted_genotype[o,1]==qqq[a],arr.ind=TRUE)})
		cleangeno <- cleangeno[-qqq1,]
	}
	write.table(cleangeno[,-1],"gen.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,na="0")
	#determine which individuals are phenotyped but not genotyped
	rrr <- setdiff(as.character(pheno[p,1]),as.character(sorted_genotype[o,1]))
	cleanpheno <- pheno[p,]
	if (length(rrr)!=0) {
		rrr1 <- sapply(1:length(rrr),function(a) {which(pheno[p,1]==rrr[a],arr.ind=TRUE)})
		cleanpheno <- cleanpheno[-rrr1,]
	}
	colnames(cleanpheno) <- c("ID",pheno.names)
	write.table(cleanpheno,"phe.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,na="-")

	#overwritten to include map duplication since the map for a 4-way cross has to be a matrix of 2 x n.markers rather than a vector of n.markers
	read.cross.karl <- function (dir, genfile, mapfile, phefile) {
		if (missing(genfile)) 
			genfile <- "gen.txt"
		if (missing(mapfile)) 
			mapfile <- "map.txt"
		if (missing(phefile)) 
			phefile <- "phe.txt"
		if (!missing(dir) && dir != "") {
			genfile <- file.path(dir, genfile)
			mapfile <- file.path(dir, mapfile)
			phefile <- file.path(dir, phefile)
		}
		geno <- as.matrix(read.table(genfile, na.strings = "0"))
		pheno <- as.matrix(read.table(phefile, na.strings = "-", header = TRUE))
		tempmap <- scan(mapfile, what = character(), quiet = TRUE)
		n.chr <- as.numeric(tempmap[1])
		n.mar <- 1:n.chr
		g <- map <- geno.data <- vector("list", n.chr)
		cur <- 2
		min.mar <- 1
		names(g) <- as.character(1:n.chr)
		for (i in 1:n.chr) {
			n.mar[i] <- as.numeric(tempmap[cur])
			cur <- cur + 1
			geno.data[[i]] <- geno[, min.mar:(min.mar + n.mar[i] - 1)]
			min.mar <- min.mar + n.mar[i]
			r <- as.numeric(tempmap[cur:(cur + n.mar[i] - 2)])
			d <- 0.25 * log((1 + 2 * r)/(1 - 2 * r)) * 100
			map[[i]] <- round(c(0, cumsum(d)), 2)
			cur <- cur + n.mar[i] - 1
			names(map[[i]]) <- tempmap[cur:(cur + n.mar[i] - 1)]
			dimnames(geno.data[[i]]) <- list(NULL, names(map[[i]]))
			cur <- cur + n.mar[i]
			g[[i]] <- list(data = geno.data[[i]], map = map[[i]])
			mar.names <- names(map[[i]])
			twodig <- grep("[Dd][1-9][0-9][Mm]", mar.names)
			onedig <- grep("[Dd][1-9][Mm]", mar.names)
			xchr <- grep("[Dd][Xx][Mm]", mar.names)
			chr.num <- NULL
			if (length(twodig) > 0) 
				chr.num <- c(chr.num, substr(mar.names[twodig], 2, 3))
			if (length(onedig) > 0) 
				chr.num <- c(chr.num, substr(mar.names[onedig], 2, 2))
			if (length(xchr) > 0) 
				chr.num <- c(chr.num, rep("X", length(xchr)))
			if (is.null(chr.num)) {
				chr.num <- length(mar.names)
				names(chr.num) <- "1"
			}
			else {
				chr.num <- table(chr.num)
			}
			m <- max(chr.num)
			if (m > sum(chr.num)/2 && m > 1) 
				names(g)[i] <- names(chr.num)[chr.num == m][1]
			if (names(g)[i] == "X" || names(g)[i] == "x") 
				class(g[[i]]) <- "X"
			else class(g[[i]]) <- "A"
		}
		n.mar1 <- sapply(g, function(a) ncol(a$data))
		n.mar2 <- sapply(g, function(a) length(a$map))
		n.phe <- ncol(pheno)
		n.ind1 <- nrow(pheno)
		n.ind2 <- sapply(g, function(a) nrow(a$data))
		if (any(n.ind1 != n.ind2)) {
			print(c(n.ind1, n.ind2))
			stop("Number of individuals in genotypes and phenotypes do not match.")
		}
		if (any(n.mar1 != n.mar2)) {
			print(c(n.mar, n.mar2))
			stop("Numbers of markers in genotypes and marker names files do not match.")
		}
		cat(" --Read the following data:\n")
		cat("\t", n.ind1, " individuals\n")
		cat("\t", sum(n.mar1), " markers\n")
		cat("\t", n.phe, " phenotypes\n")
		if (is.null(colnames(pheno))) 
			dimnames(pheno) <- list(NULL, paste("phenotype", 1:n.phe, sep = ""))
		if (max(geno[!is.na(geno)]) <= 2) 
			type <- "bc"
		else if (max(geno[!is.na(geno)]) <= 5) 
			type <- "f2"
		else type <- "4way"
		cross <- list(geno = g, pheno = pheno)
		class(cross) <- c(type, "cross")
		cross.type <- class(cross)[1]
		if (cross.type == "f2") 
			max.gen <- 5
		else if (cross.type == "bc") 
			max.gen <- 2
		else max.gen <- 14
		u <- unique(geno)
		if (any(!is.na(u) & (u > max.gen | u < 1))) 
			stop("There are strange values in the genotype data : ", paste(u, collapse = ":"), ".")
		if (type == "4way")
			for (i in 1:n.chr) cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,cross$geno[[i]]$map) 
		cross$pheno <- as.data.frame(cross$pheno)
		list(cross, FALSE)
	}

	cross <- read.cross(format="karl",genfile="gen.txt",phefile="phe.txt",mapfile="map.txt",genotypes=NULL,alleles=c("A","B","C","D"))
	names(cross$geno) <- chr.name
	#write.cross does not work for 4way crosses
	#write.cross(cross,format="csv",filestem=filestem)
	file.remove("gen.txt")
	file.remove("phe.txt")
	file.remove("map.txt")
	file.remove("newgenfile.txt")
	return(cross)
}
