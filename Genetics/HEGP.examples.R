library(parallel)
library(abind)
library(zoo)

# HEGP.example.R
# Author: Richard Mott, UCL Genetics Institute.
# (C) 2019 Richard Mott

# The function CFW.example() will download public data from Nicod et al 2016 Nature Genetics
# That paper analyses outbred CFW mice for many traits, see http://outbredmice.org. 
# The function will by default download genotype dosages just for chromosome 11 
# and the phenotypes for Platelets (after correction for covariates). 
# These data are chosen because there is a QTL for Platelets on chr11 
# (see http://mtweb.cs.ucl.ac.uk/qtls/Haem.PLT_chr.11.97045910.pdf)
# It will save these data in a Rdata object D that is organised appropriately for
# encryption and analysis by other functions in this file.
# It then simulates an approprately sized encryyption matrix and 
# performs association mapping using both original and encrypted data.
# It also computes diagnostic plots and statistics exploring the
# degree of randomness in the encrypted data (evaluated as the correlation 
# between the orioginal and encrypted dosages), and the concordance betgween the 
# logP values of the association statistics

# A pdf file of diagnostic plots called "chr11.PLT.pdf"

# Text output should be:
  
# > source("HEGP.examples.R")
# > CFW.example(mixed.model=TRUE)
# estimated heritability 0.02534315
# Built kinship dim 1329 1329
# estimated heritability 0.025049
# gwas logP mean diff 0.003140902 max diff 0.02634605 cor 0.9999849
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 7.400e-08 1.103e-03 2.158e-03 3.141e-03 3.865e-03 2.635e-02
# correlation of plaintext and ciphertext genotypes
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 1.330e-06 7.942e-03 1.676e-02 2.061e-02 2.989e-02 8.331e-02
# correlation of plaintext and ciphertext phenotypes
# [,1]
# [1,] 0.01833442
# [1] 1

CFW.example <- function( Rdata="chr11.PLT.Rdata", 
                        mixed.model = FALSE,
                        pdf.file=ifelse( mixed.model, "chr11.PLT.mm.pdf", "chr11.PLT.pdf" ),
                        seed=928349) {
  set.seed(seed)
  pdf(pdf.file)
  load(Rdata)

  if (  mixed.model & is.null(D$kinship) ) D$kinship = make.kinship(D$geno)
  
  # GWAS before encryption
  
  if ( mixed.model == FALSE ) {
    g = basic.gwas(D)
  } else {
    g = basic.mm.gwas(D)
  }
  
  # Encrypt 
  e = make.encrypter(D)
  De = encrypt.D(D,e, kinship=TRUE)
  
  # GWAS after encryption
  
  if ( mixed.model == FALSE) {
    ge = basic.gwas(De)
  } else {
    ge = basic.mm.gwas(De)
  }
  
  d = mean(abs(ge$logP-g$logP))
  dmax= max(abs(ge$logP-g$logP))
  cat("gwas plaintext vs ciphertext logP mean diff", d, "max diff" , dmax, "correlation", cor(ge$logP,g$logP), "\n" )
  print(summary(abs(ge$logP-g$logP)))
  #  r = sapply(1:ncol(D$geno),function(i,d,de) { cor(d[,i],de[,i])}, D$geno, De$geno)
  r = mapply(cor,as.data.frame(D$geno),as.data.frame(De$geno))
 
  par(mfrow=c(2,2))
  hist(r,main="(A) Cor of plaintext and ciphertext dosages")
  n = nrow(De$geno)
  z = r * sqrt( (n-2)/(1-r*r) )
  hist(z,freq=FALSE,main="(B) Z of plaintext and ciphertext dosages")
  x=seq(-3,3,0.01)
  lines(x,dnorm(x),col="red")
  qqnorm(z)
  abline(a=0,b=1,col="red",main = "(C) QQplot of Z")
  plot(ge$logP,g$logP,pch=20, cex=0.1, xlab="logP encrypted", ylab="logP original", main="(D) PLT GWAS ")
  abline(a=0,b=1,col="red")
  dev.off()
  cat("correlation of plaintext and ciphertext genotypes\n")
  print (summary(abs(r)))
  cat("correlation of plaintext and ciphertext phenotypes\n")
  print(cor(D$y,De$y))
  
  return(1)
}


# Generate a set of orthogonal matrices of max dimension <=bloocksize suitable for encrypting genotype dosage and phenotype datasets of the same dimension as the one in D. The dimension of the last matrix is reduced in order to ensure the sum of the dimensions exactly equals the number of individuals in D. If blocksize<=0 then a single block to encrypt the dataset is generated

make.encrypter <- function( D, blocksize=0 ) {
  N = nrow(D$geno)
  if ( blocksize <= 0  ) blocksize = N;
  
  start = 1
  end = min( blocksize, N)
  encrypter = list(N=N)
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

encrypt.D <- function( D, encrypter, invert=FALSE, kinship=FALSE ) {
  
  blocks = encrypter$blocks
  block = encrypter$block
  encrypted = list( y=rep( 0, length=length(D$y)), geno = matrix( 0, nrow=nrow(D$geno), ncol=ncol(D$geno)), cov=matrix(0, nrow=nrow(D$cov), ncol=ncol(D$cov)))
  colnames(encrypted$geno) = colnames(D$geno)
  colnames(encrypted$cov) = colnames(encrypted$cov)
  for( i in 1:nrow(blocks) ) {
    P = encrypter$block[[i]]
    if ( invert ) P = t(P)
    idx = blocks$start[i]:blocks$end[i]
    encrypted$y[idx] = P %*% D$y[idx]

    encrypted$geno[idx,] = P %*% D$geno[idx,]
    encrypted$cov[idx,] = P %*% D$cov[idx,]
    encrypted$map = D$map
  }
  names(encrypted$y) = names(D$y)
  rownames(encrypted$geno) = rownames(D$geno)
  rownames(encrypted$cov) = rownames(D$cov)
  if ( kinship ) {
    encrypted$kinship = make.kinship(encrypted$geno)
    colnames(encrypted$kinship) = rownames(encrypted$geno)
    rownames(encrypted$kinship) = rownames(encrypted$geno)
  }
  return( encrypted ) 
}



# create a dataset comprising a genotype dosage matrix, phenotype vector and covariate matrix

build.D <- function( y, dosages, cov=NULL, map=NULL, kinship=FALSE ) {
  
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
  
  y.s = scale(y)
  names(y.s) = names(y)
  af = apply(dosages, 2, mean)
  af = ifelse( af < 0.5, af, 1-af)
  geno = safe.scale(dosages)
  cov = safe.scale(cov)
  
  D = list( y=y.s, geno=geno, cov=cov, map=map, maf=af )
  if ( kinship ) {
    K = make.kinship(D$geno)
    D$kinship = K
  }
  return(D)
}

safe.scale <- function( mat) {
  # scale the columns of a matrix dealing safely with situation where variance of column is 0
  # if a column is not numeric then it is converted to a column of 0's
  # row and column names are preserved.
  
  m = apply(mat,2, function(x) { 
    if ( is.numeric(x)) {
      s = sd(x)
      if ( s > 0.0) {
        x = (x-mean(x))/s
      } else {
        x = rep(0,len=length(x))
      }
    }
    else {
      x = rep(0,len=length(x))
    }
    })
  colnames(m) = colnames(mat)
  rownames(m) = rownames(mat)
  return(m)
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

basic.mm.gwas <- function( D, mc.cores=1) {
  y = D$y
  if( is.null(names(y))) names(y) = as.character(1:length(y))
  genotypes = D$geno
  rownames(genotypes) = names(y)
  mm = mixed.model.gwas( y, genotypes, kinship=D$kinship, nperm=0 ) 
  return( data.frame(cbind(D$map, logP=-log10(mm$pval[1,]))))
}

# Simulate a random orthogonal matrix of dimensions m*R using the Steifel Manifold (function adapted from R package rsteifel)
rustiefel <- function (m, R=m) 
{
  X <- matrix(rnorm(m * R), m, R)
  tmp <- eigen(t(X) %*% X)
  X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
}


library(RColorBrewer)
figure1 <- function( D, N=20, m=30, pdf.file="fig1.pdf", seed=1234, idx=1501) {
  
  set.seed(seed)
  pdf(pdf.file)
  par(mar=c(3,4.1,4,2.1))
  par(mfrow=c(3,2))
  snps = sample( 1:ncol(D$geno), size=m)
  r = rustiefel(nrow(D$geno))
  geno.e = r %*% D$geno[,snps]
  y.e = r %*% D$y
  filler = matrix(0,ncol=4,nrow=N)
  g = D$geno[1:N,snps]
  
  gg=apply(g,2,function(x) { r=range(x); round(2*(x-r[1])/(r[2]-r[1]+1.0e-6))})
  yy = D$y[1:N,1]
  rr=range(yy)
  yy = 4*(yy-rr[1])/(rr[2]-rr[1]+1.0e-6) -2
  mat = cbind(yy,filler,gg )
  mat.e = cbind(y.e[1:N,1],filler,geno.e[1:N,] )
  
  
  
  image(t(mat),axes=F, col=brewer.pal(11,'RdGy'))
  title("A           Unencrypted Dosage Matrix", adj=0)
  axis(1, at=c(0,0.5), labels=c("y", "G"),font=4, tick=FALSE)
  
  image(t(mat.e),axes=F, col=brewer.pal(11,'RdGy'))
  title("B           Encrypted Dosage Matrix", adj=0)
  axis(1, at=c(0,0.5), labels=c("Py", "PG"), font=4, tick=FALSE)
  
  
  x = D$geno[,idx]
  rx=range(x); 
  gg.idx=round(2*(x-rx[1])/(rx[2]-rx[1]+1.0e-6))
  hist(gg.idx,breaks=50,col="black",main=NA,axes=FALSE,xlab=NULL,ylab=NULL,freq=TRUE)
  axis(1, at = c(0,1,2), labels = TRUE, tick = TRUE)
  title("C           Genotype Dosages", adj=0)
  box()
 
  gg.idx = r %*% D$geno[,idx]
  hist(gg.idx,main=NA,axes=FALSE,xlab=NULL,ylab=NULL,freq=FALSE)
  title("D           Encrypted Dosages", adj=0)
  nx = seq(-3,3,0.01)
  nd = dnorm(nx, mean=mean(gg.idx), sd=sd(gg.idx))
  lines(nx, nd, col="red")
  box()
  
  plot(0,type='n',axes=FALSE,ann=FALSE)

 
  qqnorm(gg.idx,main=NULL, pch=20, cex=0.5)
  abline(a=0,b=1,col="red")
  title("E           QQ Plot of Encrypted Dosages", adj=0)
  
  dev.off()
}
# mixed model code


estimate.mixed.model <- function( y, kinship, make.positive=TRUE ) {
  y = y[!is.na(y)]
  if ( length(y) < 1 ) return(NULL)
  
  use.ids = intersect( names(y), colnames(kinship))
  match.kinship = match( use.ids, colnames(kinship), nomatch=0)
  K = kinship[match.kinship,match.kinship]
  K.eigen.trunc = K.eigen = eigen(K,symmetric=TRUE)
  if ( make.positive ) {
    w.eigen = which( K.eigen$values/K.eigen$values[1] > 1.0e-8 )
    neigen = length(w.eigen) # number of non-trivial principal components
    K.eigen.trunc = K.eigen
    K.eigen.trunc$values = K.eigen.trunc$values[w.eigen]
    K.eigen.trunc$vectors = K.eigen.trunc$vectors[,w.eigen]
  }
  match.phen = match( use.ids, names(y), nomatch=0)
  y = y[match.phen]
  y = scale(y)
  
  z = t(K.eigen.trunc$vectors) %*% y # transformed phenotype
  zz = z*z
  lambda = K.eigen.trunc$values
  
  opt.theta1 = optimise( interval=c( 0,1 ),
                         f =function(theta, zz, lambda ) { # log likelihood for mvn 
                           u = theta[1]*lambda+1-theta[1]
                           u = ifelse ( abs(u)<1.0e-10, 1.0e-10 , u )
                           return(sum(zz/u) + sum( log(u)))
                         }, zz=zz, lambda=lambda)
  
  vg = opt.theta1$minimum[1]
  ve = 1-vg
  cat("estimated heritability", vg,"\n")
  E = K.eigen.trunc$vectors
  Lam = vg*K.eigen.trunc$values +ve
  V = Lam^(-0.5)    
  inv.sqrt.sigma = ( E %*% diag(V) ) %*% t(E)
  
  mixed.model.data = list( y=y, K=K, vg=vg, ve=ve, inv.sqrt.sigma=inv.sqrt.sigma, eigen=K.eigen.trunc )
  return(mixed.model.data)
}


mixed.model.gwas <- function( y, genotypes, kinship=NULL, nperm=0 ) {
  # y is an  n-vector of phenotypes. 
  # genotypes is a matrix of genotype dosages, n rows and p columns
  
  # Each element of y should be named (ie names(y) returns meaningful IDs as should rownames(genotypes)) 
  # so that the code can make sure the phenotypes and genotypes are correctly aligned and ordered. 
  # K is the optional kinship matrix. If omitted then it is computed from the genotypes
  # nperm is the number of permutations to perform to determine genomewise significance.
  
  genos = apply( genotypes, 2, function( g ) {
    s = sd(g, na.rm=TRUE)
    mu = mean(g, na.rm=TRUE)
    if ( s > 0 ) {
      g = ifelse ( is.na(g), 0.0, (g-mu)/s )
    } else {
      g = rep(0, length(g));
    }
    return(g)
  })
  
  
  if ( is.null(kinship)) {
    kinship = make.kinship( genos )
  }
  
  mmd = estimate.mixed.model( y, kinship )
  
  use.ids = rownames(mmd$y)
  genos = genotypes[match( use.ids, rownames(genotypes), nomatch=0),]
  
  if ( nrow(genos) != length(use.ids)) {
    cat( "ERROR sample ids in genotypes do not match phenotypes\n")
    return(NULL);
  }
  
  mm.transformation = mmd$inv.sqrt.sigma    
  mm.transformed.y = mm.transformation %*% mmd$y
  mm.transformed.g = mm.transformation %*% genos
  
  mm.gwas.cor =  cor( mm.transformed.y, mm.transformed.g )
  n = length(mm.transformed.y)
  df = n-2
  t.stat = mm.gwas.cor*sqrt((df-2)/(1.0e-10+sqrt(1-mm.gwas.cor*mm.gwas.cor)))
  pval = 2*pt( abs(t.stat), df, low=FALSE )
  
  if ( nperm > 0 ) {
    perms = matrix( NA, ncol=nperm, nrow=n)
    for( i in 1:nperm ) {
      perms[,i] = sample(mm.transformed.y,replace=FALSE)
    }
    
    mm.gwas.cor.perm = cor( perms, mm.transformed.g )
    t.stat.perm = mm.gwas.cor.perm*sqrt((df-2)/(1.0e-10+sqrt(1-mm.gwas.cor.perm*mm.gwas.cor.perm)))
    pval.perm = 2*pt( abs(t.stat.perm), df, low=FALSE)
    pval.perm.empirical = sapply( 1:ncol(pval.perm), function( i, pval.perm, pval ) { mean(pval[i] > pval.perm[,i]) }, pval.perm, pval ) 
    results = list( pval=pval, pval.perm=pval.perm, pval.perm.empirical=pval.perm.empirical )
  } else {
    results = list( pval=pval)
  }
  
  return(results)
}

make.kinship <- function( genotypes, scale.genos=FALSE ) {
  
  if ( scale.genos ) {
    genos = apply( genotypes, 2, function( g ) {
      s = sd(g, na.rm=TRUE)
      mu = mean(g, na.rm=TRUE)
      if ( s > 0 ) {
        g = ifelse ( is.na(g), 0.0, (g-mu)/s )
      } else {
        g = rep(0, length(g));
      }
      return(g)
    })
    rownames(genos) = rownames(genotypes)
    genotypes = genos
  }
  
  kinship = genotypes %*% t(genotypes)/(ncol(genotypes)-1)
  colnames(kinship) = rownames(genotypes)
  rownames(kinship) = rownames(genotypes)
  
  cat("Built kinship dim", dim(kinship), "\n")
  
  return(kinship)
}


