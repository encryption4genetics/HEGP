cayley.transform <- function ( A ) {
	n= nrow(A)
	I = diag(n)
	return((I-A) %*% solve( I+A))
}

skew.symmetric <- function( n=10) {
	r = rnorm(n*(n-1)/2)
	A = matrix(0,nrow=n, ncol=n)
	A[upper.tri(A,diag=FALSE)] = r
	A = A-t(A)
	return(A)
}

library(rstiefel)

rustiefel <- function (m, R) 
{
  X <- matrix(rnorm(m * R), m, R)
  tmp <- eigen(t(X) %*% X)
  X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
}

simulate.P <- function( n=100, m=1000) {
	
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

require(RColorBrewer)

diverge.color <- function(data,pal_choice="RdGy",centeredOn=0){
  nHalf=50
  Min <- min(data,na.rm=TRUE)
  Max <- max(data,na.rm=TRUE)
  Thresh <- centeredOn
  pal<-brewer.pal(n=11,pal_choice)
  rc1<-colorRampPalette(colors=c(pal[1],pal[2]),space="Lab")(10)
  for(i in 2:10){
    tmp<-colorRampPalette(colors=c(pal[i],pal[i+1]),space="Lab")(10)
    rc1<-c(rc1,tmp)
  }
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  return(rc1)
  
}
simulate.image<- function( N=20, M=50, pdf.file="image.pdf") {
	
	rc1 <-diverge.color(seq(-3,3,0.1))
	if (!is.null(pdf.file)) pdf(pdf.file)
	dose = matrix( rbinom(size=2, n=N*M, p=0.5),nrow=N,ncol=M)
	y = matrix( rnorm(N), nrow=1)
	y.dose = matrix(0,nrow=N,ncol=M+10)
	y.dose[,1] = y
	y.dose[,11:ncol(y.dose)] = dose
	P = rustiefel(m=N,R=N)
	P.dose = P %*% dose
	P.y.dose= P %*% y.dose
	P.y.dose[,2:10]=0
	y.dose[,2:10] = 0
	par(mfrow=c(2,2))
	image( t(dose), axes=F,asp=1)
	image( t(y.dose), axes=F, asp=1, col=rc1)
	image( t(P.dose), axes=F, asp=1)
	image( t(P.y.dose), axes=F, asp=1, col=rc1)
	if (!is.null(pdf.file)) dev.off()
}