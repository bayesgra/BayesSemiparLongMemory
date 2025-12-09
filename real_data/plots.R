load("Nile/Nile.RData")

n = length(NileMin)
FF = abs(fft(NileMin)/sqrt(n))^2
P = (4/n)*FF[1:((n+1)/2)] # Only need the first (n/2)+1 values of the FFT result.
f = (0:((n+1)/2-1))/n # this creates harmonic frequencies from 0 to .5 in steps of 1/128.
par(mfrow=c(1,1))
plot(f, P, type="l",ylab="I",xlab=expression(lambda)) # This plots the periodogram; type = “l” creates a line plot. Note: l is lowercase L, not number 1.

y = log(P)
x = -2*log(f)

par(mfrow=c(1,1))
plot(x,y,xlab=expression(log(lambda)),ylab="log(I)")
abline(lm(y[2:length(y)]~x[2:length(x)]),col="red",lwd=2)

plot(x[2:length(x)],lm(y[2:length(y)]~x[2:length(x)])$residuals,ylab="Residuals",
     xlab=expression(log(lambda)))
abline(h=0,col="red")

plot(NhemiTemp,col="cadetblue")

data(NileMin)

data <- na.omit(data)
dim(data)
head(NileMin)
length(NileMin)

par(mfrow=c(1,1))
par(mar = c(3, 4.5, 2, 2) + 0.1)
par(oma = c(0, 0, 0, 0))
plot(NileMin,col="cadetblue",type="l",ylab="Annual Water Levels",xlab="")

acf(NileMin,main="",xlab="")

#####

load("S&P500/SP500.RData")

data <- read.csv("SP500.csv")

n = length(log_returns)
FF = abs(fft(log_returns)/sqrt(n))^2
P = (4/n)*FF[1:((n+1)/2)] # Only need the first (n/2)+1 values of the FFT result.
f = (0:((n+1)/2-1))/n # this creates harmonic frequencies from 0 to .5 in steps of 1/128.

par(mfrow=c(1,1))
plot(f, P, type="l",ylab="I",xlab=expression(lambda)) # This plots the periodogram; type = “l” creates a line plot. Note: l is lowercase L, not number 1.

y = log(P)
x = -2*log(f)

par(mfrow=c(1,1))
plot(x,y,xlab=expression(log(lambda)),ylab="log(I)")
abline(lm(y[2:length(y)]~x[2:length(x)]),col="red",lwd=2)

plot(x[2:length(x)],lm(y[2:length(y)]~x[2:length(x)])$residuals,ylab="Residuals",
     xlab=expression(log(lambda)))
abline(h=0,col="red")

par(mfrow=c(1,1))
par(mar = c(3, 4.5, 2, 2) + 0.1)
par(oma = c(0, 0, 0, 0))
plot(as.Date(data$DATE),data$SP500,type="l",
     xlab="",ylab="Prices",col="darkgoldenrod2")
log_returns <- diff(log(as.numeric(data$SP500)), lag=1)
plot(as.Date(data$DATE)[-1],log_returns,type="l",
     xlab="",ylab="Log-returns",col="darkgoldenrod2")

acf(na.omit(log_returns),xlab="",main="")

################

load("Temperatures/global_temp.RData")

library(zoo)
temp <- read.csv("HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv")

par(mfrow=c(1,1))
par(mar = c(3, 4.5, 2, 2) + 0.1)
par(oma = c(0, 0, 0, 0))
plot(as.yearmon(temp$Time), temp$Anomaly..deg.C.,type="l",ylab="Temperatures C Degrees",col="dodgerblue",xlab="")

str(temp)


tempdiff = diff(temp$Anomaly..deg.C.,differences=1)
plot(as.yearmon(temp$Time)[-1],tempdiff,type="l",col="dodgerblue",ylab="Differences",xlab="")

acf(tempdiff,xlab="",main="")

n = length(tempdiff)
FF = abs(fft(tempdiff)/sqrt(n))^2
P = (4/n)*FF[1:((n+1)/2)] # Only need the first (n/2)+1 values of the FFT result.
f = (0:((n+1)/2-1))/n # this creates harmonic frequencies from 0 to .5 in steps of 1/128.

par(mfrow=c(1,1))
plot(f, P, type="l",ylab="I",xlab=expression(lambda)) # This plots the periodogram; type = “l” creates a line plot. Note: l is lowercase L, not number 1.

y = log(P)
x = -2*log(f)

par(mfrow=c(1,1))
plot(x,y,xlab=expression(log(lambda)),ylab="log(I)")
abline(lm(y[2:length(y)]~x[2:length(x)]),col="red",lwd=2)

plot(x[2:length(x)],lm(y[2:length(y)]~x[2:length(x)])$residuals,ylab="Residuals",
     xlab=expression(log(lambda)))
abline(h=0,col="red")
