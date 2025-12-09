library(fracdiff)
library(LongMemoryTS)
T<-500
m<-floor(1+T^0.8)
d=0.4
series<-fracdiff.sim(n=T, d=d)$series
gph(X=series,m=m)

### AR=0, MA=0

data <- load("BayesARFIMA-049.RData")
dhat <- c()
N <- ncol(y_mat)
m<-floor(1+N^0.8)

for(i in 1:20)
{
  y <- y_mat[i,]
  dhat[i] <- gph(X=y,m=m,l=1)
}
dhat
c(quantile(dhat,0.025),mean(dhat),quantile(dhat,0.975))
rm(list=ls())

### AR=0.8, MA=0.1

data <- load("BayesARFIMA-049_ARhighMAlow.RData")
dhat <- c()
N <- ncol(y_mat)
m<-floor(1+N^0.8)

for(i in 1:50)
{
  y <- y_mat[i,]
  dhat[i] <- gph(X=y,m=m,l=1)
}
dhat
c(quantile(dhat,0.025,na.rm=T),mean(dhat,na.rm=T),quantile(dhat,0.975,na.rm=T))
rm(list=ls())

### AR=0.5, MA=0.1

data <- load("BayesARFIMA-049_ARmediumMAlow.RData")
dhat <- c()
N <- ncol(y_mat)
m<-floor(1+N^0.8)

for(i in 1:50)
{
  y <- y_mat[i,]
  dhat[i] <- gph(X=y,m=m,l=1)
}
dhat
c(quantile(dhat,0.025,na.rm=T),mean(dhat,na.rm=T),quantile(dhat,0.975,na.rm=T))
rm(list=ls())

### AR=0.1, MA=0.5

data <- load("BayesARFIMA-049_ARlowMAmedium.RData")
dhat <- c()
N <- ncol(y_mat)
m<-floor(1+N^0.8)

for(i in 1:20)
{
  y <- y_mat[i,]
  dhat[i] <- gph(X=y,m=m,l=1)
}
dhat
c(quantile(dhat,0.025,na.rm=T),mean(dhat,na.rm=T),quantile(dhat,0.975,na.rm=T))
rm(list=ls())

### AR=0.1, MA=0.8

data <- load("BayesARFIMA-049_ARlowMAhigh.RData")
dhat <- c()
N <- ncol(y_mat)
m<-floor(1+N^0.8)

for(i in 1:50)
{
  y <- y_mat[i,]
  dhat[i] <- gph(X=y,m=m,l=1)
}
dhat
c(quantile(dhat,0.025,na.rm=T),mean(dhat,na.rm=T),quantile(dhat,0.975,na.rm=T))
rm(list=ls())
