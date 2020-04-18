library(parallel)

# Simulate a random orthogonal matrix of dimensions m*R using the Steifel Manifold 
rustiefel <- function (m, R=m) 
{
  X <- matrix(rnorm(m * R), m, R)
  tmp <- eigen(t(X) %*% X)
  X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
}

subsample.dataset<- function( D, n=1000 ) {

        N = length(D$y)
	if ( n > N ) {
	   return(D)
	} else {
	  idx = sample( N, size=n, replace=FALSE)
      	  D.subset = list( y=D$y[idx], geno=D$geno[idx,], cov=D$cov[idx,], map=D$map )
	  return(D.subset)
	}
}
	     
# Generate a set of orthogonal matrices of max dimension <=bloocksize suitable for encrypting genotype dosage and phenotype datasets of the same dimension as the one in D. The dimension of the last matrix is reduced in order to ensure the sum of the dimensions exactly equals the number of individuals in D. If blocksize<=0 then a single block to encrypt the dataset is generated

make.encrypter <- function( D, blocksize=0 ) {
	     N = nrow(D$geno)
	     if ( blocksize <= 0  ) blocksize = N;

	     start = 1
	     end = min( blocksize, N)
	     encrypter = list()
	     block = list()
	     b = 1
	     df = NULL
	     while ( end <= N ) {
		    if ( end > N-100 ) end = N
	            bsize = end-start+1
		    block[[b]] = rustiefel( bsize, bsize )
		    df = rbind( df, c(start, end, bsize ))
		    start = end+1
		    if ( end == N ) {
		       end = N+1
		    } else {
		       end = min( start + blocksize-1, N )
		    }
		    b = b+1
		 }           
	     df = data.frame(df)
	     names(df) = c( "start", "end", "size" )
	     encrypter$blocks = df
	     encrypter$block = block
	     return ( encrypter )
}

# Encypt a dataset D using the encrypter

encrypt.D <- function( D, encrypter ) {

	    blocks = encrypter$blocks
	    block = encrypter$block
	    encrypted = list( y = rep( 0, length=length(D$y)), geno = matrix( 0, nrow=nrow(D$geno), ncol=ncol(D$geno)), cov=matrix(0, nrow=nrow(D$cov), ncol=ncol(D$cov)))
	    colnames(encrypted$geno) = colnames(D$geno)
	    colnames(encrypted$cov) = colnames(encrypted$cov)
	    for( i in 1:nrow(blocks) ) {
	    	 P = encrypter$block[[i]]
		 idx = blocks$start[i]:blocks$end[i]
	    	 encrypted$y[idx] = P %*% D$y[idx]
		 encrypted$geno[idx,] = P %*% D$geno[idx,]
		 encrypted$cov[idx,] = P %*% D$cov[idx,]
		 encrypted$map = D$map
	    }
	    return( encrypted ) 
}
 
# create a dataset comprising a genotype dosage matrix, phenotype vector and covariate matrix

build.D <- function( y, dosages, cov=NULL, map=NULL ) {

  y = y[!is.na(y)]
	ids = intersect( names(y) , rownames(dosages))
	y = y[match(ids, names(y), nomatch=0)]
        dosages = dosages[match(ids, rownames(dosages), nomatch=0),]
	N = length(y)

	if ( !is.numeric(y)) {
	   warning( "y is not numeric")
	   return(NULL)
	}	   

	if ( is.null(cov)) {
	    cov = matrix(1, nrow=N, ncol=1)
	    colnames(cov) = c("intercept")
	}

	if ( sum(is.na(y)) + sum(is.na(dosages)) + sum(is.na(cov)) > 0 ) {
	   warning( "data contain missing values" )
	   return(NULL)
	}
	
	if ( N != nrow(dosages ) | ( !is.null(cov) & N != nrow(cov) )) {
	   warning( "dimensions incompatible")
	   return(NULL)
	}

	y = scale(y)
	geno = scale(dosages)
	cov = scale(cov)

	D = list( y=y, geno=geno, cov=cov, map=map )
	return(D)
}

# Perform a simple quantitative trait GWAS without a mixed model, regressing phenotype on genotype dosage.

basic.gwas <- function( D, mc.cores=10 ) {
	   
           to = floor((1:mc.cores)*ncol(D$geno)/mc.cores)
	   to[mc.cores] = ncol(D$geno)
	   start = c(1, to[1:(mc.cores-1)]+1)
	   lp = mclapply( 1:mc.cores, function( segment ) {
#	   cat("segment", segment, start[segment],to[segment],"\n")
	   r = as.numeric(cor( D$y, D$geno[,start[segment]:to[segment]]))
	   r2 = r*r
	   n = length(D$y)
	   t = r * sqrt( (n-2)/(1-r2) )
	   logP = -pt( abs(t), df=n-2, lower=FALSE, log=TRUE)/log(10) 
	   return(logP) })
	   logP = unlist(lp)
	   return( data.frame(cbind(D$map, logP)))
}


# simulate a sequence of orthogonal matrices for Figure 1.

# n = dimension of matrices (ie number of individuals)
# m = number of simulations per lambda point

simulate.P <- function( n=100, m=1000, pdf.file="SimFigure.pdf") {
	
	P = rustiefel(n,n)
	e = eigen(P)
	inv = solve(e$vectors) + 0+0i
	vec = e$vectors + 0+0i
	lambda = seq(0,10,0.01)
	rn = matrix(0, nrow=length(lambda),ncol=m)
	ru = matrix(0, nrow=length(lambda),ncol=m)
	i=1
	for( x in lambda ) {
		PP = Re(vec %*% (diag(e$values^x) %*% inv ))
		for(j in 1:m) {
			y = rnorm(n)
			rn[i,j] = cor(y, PP %*% y)
			y = runif(n)<0.5
			ru[i,j] = cor(y, PP %*% y)
		}
		i = i+1
	}
	mn = apply(rn,1, mean)
	mu = apply(ru,1, mean)
	plot(lambda, mn, t="l", xlab=expression(lambda), ylab="mean correlation")
	lines(lambda,mu, col="red")
	abline(h=0, col="grey")
	abline(v=1:9)
	return(list(rn, ru))
}

correlation.P <- function( D, n=1000, m=1000, pdf.file="SimFigureReal.pdf" ) {

	Ds = subsample.dataset(D,n)      	
	P = rustiefel(n,n)
	e = eigen(P)
	inv = solve(e$vectors) + 0+0i
	vec = e$vectors + 0+0i
	lambda = seq(0,10,0.01)
#	lambda = seq(0,0.02,0.01)
	w.snps = sample(nrow(Ds$geno),size=m,replace=FALSE)
	snps = Ds$geno[,w.snps]
	L = length(lambda)
	rrn = mclapply( 1:L, function(k, snps, lambda, inv, e, m ) {
		cat("k ",k,lambda[k],"\n")
		PP = Re(vec %*% (diag(e$values^lambda[k]) %*% inv ))
		snps.enc = PP %*% snps
		
		r = sapply( 1:m, function(i) { cor( snps[,i], snps.enc[,i] )})
	},  snps, lambda, inv, e, m, mc.cores=20)
	browser();
	rn = do.call("rbind", rrn )
	mn = apply(rn,1, mean)
	sd.mn = apply(rn,1,sd)
	pdf(pdf.file)
	par(mfrow=c(3,1))
	plot(lambda, mn, t="l", xlab=expression(lambda), ylab="mean correlation")
	lines(lambda,mn+sd.mn, col="red")
	lines(lambda,mn-sd.mn, col="red")
	abline(h=0, col="grey")
	abline(v=1:9)
	dev.off()
	return(rn)
}

make.converge.genos <- function( dir="/Net/sparse/data/sanja/WTCHG/1--SVproject/Analyses/CONVERGE/4--QTLAnalysis/0--taggingSNPsCONVERGE/dosage/byChr/", infile="D.RData", outfile="converge.dosages.RData") {

files = list.files( path=dir, pattern="snpDataCONVERGEchr", full=TRUE)

genos = NULL
map = NULL
ids = NULL

for( chr in 1:23 ) {
    f = paste0( dir, "/snpDataCONVERGEchr", chr, ".RData")
    load(f) 
    cat(f, "\n")
    genos = cbind( genos, chr.data$genotype)
    map = rbind( map, chr.data$alleles[,c(1,2)] )
}

load(infile)
names(map) = c("chr", "bp")
ids = intersect( rownames(D$y), rownames(genos))
D$y = D$y[match(ids,rownames(D$y),nomatch=0),]
D$geno=genos[match(ids, rownames(genos),nomatch=0),]
D$map= map

save(D,file=outfile)
}

check.logistic.regression <- function( D, encoder, gwas.results, N=10000) {
			   w = which( gwas.results$logP>3)
			   idx = sample( 1:nrow(gwas.results), size=N, replace=FALSE)
			   idx = unique(sort(c(w,idx)))
			   yy = ifelse(D$y > 0 , 1, 0) 
			   D.subset = list(y=yy,map=D$map,cov=D$cov,geno=D$geno[,idx])
			   
			   D.enc.subset = list(y= encoder %*% yy,map=D$map,cov=NULL,geno=encoder %*% D$geno[,idx])
			   unit = rep(1,length=length(D$y))
			   unit.enc = encoder %*% unit

			   D$f0 = glm(D$y ~ 1, family="binomial")
			   D.enc.subset$f0.enc = glm(D.enc.subset$y ~ unit.enc-1, family="binomial")
			   logP = sapply( 1:length(idx), function( i, D.subset, D.enc.subset ) {
			   	browser()
			   	   f = glm( D$y ~ D$geno[,i], family="binomial")
				   a = anova( D$f0, f, test="Chisq")
				   pval = a[1,5]
				   f.enc = glm( D.enc$y ~ unit.enc+D.enc$geno[,i]-1, family="binomial")
				   a = anova( D.enc$f0, f.enc, test="Chisq")
				   pval.enc = a[1,5]
				   return(-log10(c(pval, pval.enc))) }, D.subset, D.enc.subset )
		           results = data.frame( D$map[idx,], logP=logP[,1], logP.enc=logP[,2])
			   return(results)
}

fit.logistic <- function( y, X, beta0=c(0,0)) {
  results = optim( par=beta0, fn=logistic.log.likelihood, 
                   gr=gr.logistic.log.likelihood, lower=c( -30, -30), upper=c( 30, 30), hessian=TRUE, method="L-BFGS-B", X=X, y=y )
}

gr.logistic.log.likelihood <- function( beta,X,y ) {
  lp = X %*% beta
  e = exp(lp)
  u = e/(1+e)
  return(- t(y-u) %*% X)
}

logistic.log.likelihood <- function( beta,X,y ) {
  lp = X %*% beta
  e = exp(lp)
  L = -t(y) %*% lp + sum( log(1+e))
  #  cat(beta, L, "\n")
  return(L)
}

test.logistic<-function( N=500 , beta0=0, beta1=2, p=0.2) {
  
  x = rbinom(N,size=2,p=p)
  X = cbind(rep(1,length=N), x)
  beta = c(beta0,beta1)
  u = X %*% beta
  q = exp(u)/(1+exp(u))
  y = runif(N) < q
  g = glm( y ~ x, family="binomial")
  f = fit.logistic( y, X)
  X1 = as.matrix(X[,1],ncol=1)
  f0 = fit.logistic( y, X1, beta=0)
  chi.s = 2*(logistic.log.likelihood(f$par,X,y)-logistic.log.likelihood(f0$par,X1,y))
  cat("glm", g$coef, "\n")
  cat("optim", f$par, "\n")
  cat("chis", chi.s, "\n")
  browser()
  P = rustiefel(N)
  y.P = P %*% y
  X.P = P %*% X
  X1.P = P %*% X1
  f.P = fit.logistic( y.P, X.P)
  cat("optim.P", f.P$par, "\n")
  f0.P = fit.logistic( y.P, X1.P, beta=0)
  chi.s.P = 2*(logistic.log.likelihood(f.P$par,X.P,y.P)-logistic.log.likelihood(f0.P$par,X1.P,y.P))
  cat("chis.P", chi.s.P, "\n")
  seq1 = seq(-2,2,0.001)
  seq2 = seq(1,4,0.001)
  mat = matrix(0, ncol=length(seq2), nrow=length(seq1))
  mat.P = matrix(0, ncol=length(seq2), nrow=length(seq1))
  i = 0
  for( b0 in seq1) {
    i <- i+1
    j = 0
    for(b1 in seq2) {
      j = j+1
      beta = c(b0,b1)
      mat[i,j] = logistic.log.likelihood(beta,X,y)
      mat.P[i,j] = logistic.log.likelihood(beta,X.P,y.P)
    }
  }
  par(mfrow=c(2,1))
  image(mat)
  image(mat.P)
}

logistic.log.likelihood <- function( y, x, b0, b1) {
  u = b0+b1*x
  uu = exp(u)
  luu = log(1+uu)
  sluu = sum(luu)
  uu1 = 1/uu
  L = sum(y*u) - sluu
  N = length(y)
  u1 = 
  DL0 = sum(y) - sum(sluu)
  DL1 = sum(y*x) - sum(x*uu*uu1)
}

# This function CFW.example() will download public data from Nicod et al 2016 Nature Genetics
# The paper analyses 200 outbred CFW mice for many traits, see http://outbredmice.org. 
# The function will by default download genotype dosages just for chromosome 11 
# and the phenotypes for Platelets (after correction for covariates). 
# These data are chosen because there is a QTL for Platelets on chr11 
# (see http://mtweb.cs.ucl.ac.uk/qtls/Haem.PLT_chr.11.97045910.pdf)
# It will save these data in a Rdata object D that is organised appropriately for
# encryption and analysis by other functions in this file.
# It then simulates an approprately sized encryyption matrix and 
# performs association mapping using both original and encrypted data.
# It also computes so,me diagnostic plots and statistics exploring the
# degree of randomness in the encrypted data (evaluated as the correlation 
# between the orioginal and encrypted dosages), and the concordance betgween the 
# logP values of the association statistics

CFW.example <- function(dosage.url="http://mtweb.cs.ucl.ac.uk/dosages/chr11.prunedgen.final.maf001.0.98.RData", 
                              phenotype.url="http://mtweb.cs.ucl.ac.uk/phenotypes/CFW_residuals.txt",
                              phenotype="Haem.PLT",
                              Rdata="chr11.PLT.Rdata",
                              pdf.file="chr11.PLT.pdf",
                              seed=123456789) {
  set.seed(seed)
  if ( !file.exists(Rdata)) {
    load(url(dosage.url))
    u = unlist(nameList)
    u1=sub("recal.reheadered.bam", "", u)
    u1=sub("_recal.reheadered.bam", "", u)
    u2=sub("___","/",u1)
    colnames(pruned_dosages) = u2
  
    p = read.delim(url(phenotype.url))
    y = p[[phenotype]]
    names(y) = p$Animal_ID
    y = y[!is.na(y)]
  
    D = build.D( y=y, dosages=t(pruned_dosages), map=pruned_pos )
    save(D,file=Rdata)
  }
  else {
    load(Rdata)
  }
browser()
  af = apply(D$geno,2,function(x) { mean(x>0)/2})
  medf= ( af <0.95 & af > 0.05 )
  
 D$geno = D$geno[,medf]
 D$map = D$map[medf,]
  g = basic.gwas(D)
  e = make.encrypter(D)
  De = encrypt.D(D,e)
  ge = basic.gwas(De)
  d = mean(abs(ge$logP-g$logP))
  dmax= max(abs(ge$logP-g$logP))
  cat(phenotype, "gwas mean diff", d, "max diff" , dmax, "\n" )
  print(summary(abs(ge$logP-g$logP)))
#  r = sapply(1:ncol(D$geno),function(i,d,de) { cor(d[,i],de[,i])}, D$geno, De$geno)
  r = mapply(cor,as.data.frame(D$geno),as.data.frame(De$geno))
  pdf(pdf.file)
  par(mfrow=c(2,2))
  hist(r,main="(A) Cor of original vs encrypted dosages")
  n = nrow(De$geno)
  z = r * sqrt( (n-2)/(1-r*r) )
  hist(z,freq=FALSE,main="(B) Z of original vs encrypted dosages")
  x=seq(-3,3,0.01)
  lines(x,dnorm(x),col="red")
  qqnorm(z)
  abline(a=0,b=1,col="red",main = "(C) QQplot of Z")
  plot(ge$logP,g$logP,pch=20, cex=0.1, xlab="logP encrypted", ylab="logP original", main="(D) PLT GWAS ")
  abline(a=0,b=1,col="red")
  dev.off()
  cat("correlation of genotypes\n")
  print (summary(abs(r)))
  cat("correlation of phenotypes\n")
  print(cor(D$y,De$y))
  
}
	