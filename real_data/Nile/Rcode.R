# https://medium.com/@karakulath.suraj/quantifying-climate-change-f304eff918cf

library(longmemo)

data(NileMin)

data <- na.omit(data)
dim(data)
head(NileMin)
length(NileMin)

par(mfrow=c(1,1))
par(mar = c(4, 4.5, 2, 1) + 0.1)
par(oma = c(0, 0, 0, 0) + 0.1)
plot(NileMin,col="cadetblue",type="l",ylab="Annual Water Levels",xlab="")

acf(NileMin,main="",xlab="")

n = length(NileMin)
FF = abs(fft(NileMin)/sqrt(n))^2
P = (4/n)*FF[1:((n+1)/2)] # Only need the first (n/2)+1 values of the FFT result.
f = (0:((n+1)/2-1))/n # this creates harmonic frequencies from 0 to .5 in steps of 1/128.

par(mfrow=c(1,2))
par(mar = c(5, 4.5, 2, 2) + 0.1)
plot(f, P, type="l",ylab="I",xlab=expression(lambda)) # This plots the periodogram; type = “l” creates a line plot. Note: l is lowercase L, not number 1.

y = log(P)
x = -2*log(f)

#par(mfrow=c(1,1))
plot(x,y,xlab=expression(log(lambda)),ylab="log(I)")
abline(lm(y[2:length(y)]~x[2:length(x)]),col="red",lwd=2)
par(mfrow=c(2,2))
plot(lm(y[2:length(y)]~x[2:length(x)]))

par(mfrow=c(1,2))
plot(x,y,xlab=expression(log(lambda)),ylab="log(I)")
abline(lm(y[2:length(y)]~x[2:length(x)]),col="red",lwd=2)
plot(x[2:length(x)],lm(y[2:length(y)]~x[2:length(x)])$residuals,ylab="Residuals",
     xlab=expression(log(lambda)))
abline(h=0,col="red")

##### Frequentist estimate

library(fracdiff)
library(LongMemoryTS)

m<-floor(1+n^0.8)
dhat <- c()
for(i in 2:m)
{
  dhat[i] <- gph(X=NileMin,m=i,l=1)
}
mean(dhat[-1])
quantile(dhat[-1],probs=c(0.025,0.975))
#[1] 0.3857184
#[1] 0.4274209
#  0.3647384 0.5897619 

dhat2 <- fdGPH(NileMin)
dhat2
#[1] 0.5038294
# sdas 0.1570167
# sdred 0.1420154

library(tseries)
adf.test(NileMin)
#Dickey-Fuller = -5.0255, Lag order = 8, p-value = 0.01
#alternative hypothesis: stationary

save.image("Nile.RData")

########## Semiparametrc Bayes

load("Nile.RData")
source("longmemory_DP.R")  

library(arfima)
#library(LongMemoryTS)
library(longmemo)
library(truncnorm)
library(fracdiff)
library(LongMemoryTS)

N <- length(NileMin)
m<-floor(1+N^0.8)

lambda <- 2 * 0:N * pi / N
lambda <- lambda[1:(N/2)]
n<-floor(N/2)
lambdaj<-2*pi/N*1:n
m<-floor(1+N^0.8)
l<-1
xj <- -2*log(lambdaj[l:m])

range01<-function(x){(x-min(x))/(max(x)-min(x))}

I  <- abs(fft(NileMin)/sqrt(N))^2
P  <- (4/N) * I[1:(N/2)]

logIj <- log(P[-1])            # remove zero frequency
Yj    <- logIj[l:m]
y     <- Yj

n <- length(y)

x  <- -2 * log(lambdaj[l:m])
x  <- cbind(1, x)      # intercept + slope
p  <- ncol(x)

# Rescale t = scaled x_j (or use omega directly)
t <- (x[,2] - min(x[,2])) / (max(x[,2]) - min(x[,2]))

ynew <- y[1:150]
omeganew <- lambdaj[l:(l+149)]
xnew <- x[1:150, ]
tnew <- t[1:150]

obj <- longmemory_dp_sampler(
  Y       = ynew,
  omega   = omeganew,
  H       = 20,            
  runs    = 30000,
  burn    = 0.2 * 30000,
  thin    = 1,
  a_V     = 1,
  b_V     = 2,
  a_sigma = 0.1,
  b_sigma = 0.1,
  a_xi    = 0.1,
  b_xi    = 0.1,
  d_init  = 0,
  tau_d   = 0.05,
  display = 1000
)

par(mfrow=c(1,1))

# posterior of the slope (memory parameter)
plot(density(obj$d), main="Posterior density of d")
abline(v = mean(obj$d), col="red", lwd=2)

mean(obj$d)
quantile(obj$d, probs=c(0.025, 0.975))

# trace plot
plot(obj$d, type="l", main="Trace plot of d")

save.image("Nile.RData")

################ Parametric Bayesian ARFIMA

if (!require(nsarfima))install.packages('nsarfima',repos="https://cran.csiro.au/")
if (!require(arfima))install.packages('arfima',repos="https://cran.csiro.au/")
if (!require(MASS))install.packages('MASS',repos="https://cran.csiro.au/")
if (!require(DataCombine))install.packages('DataCombine',repos="https://cran.csiro.au/")
if (!require(FastGP))install.packages('FastGP',repos="https://cran.csiro.au/")
if (!require(spam))install.packages('spam',repos="https://cran.csiro.au/")
if (!require(mvtnorm))install.packages('mvtnorm',repos="https://cran.csiro.au/")
if (!require(tmvtnorm))install.packages('tmvtnorm',repos="https://cran.csiro.au/")
if (!require(fracdiff))install.packages('fracdiff',repos="https://cran.csiro.au/")
if (!require(sarima))install.packages('sarima',repos="https://cran.csiro.au/")
if (!require(gsl))install.packages('gsl',repos="https://cran.csiro.au/")
if (!require(mcmcsae))install.packages('mcmcsae',repos="https://cran.csiro.au/")

library(fracdiff)
library(arfima)
library(nsarfima)
library(MASS)
library(FastGP)
library(mvtnorm)
library(spam)
library(truncnorm)
#library(gsl)
library(DataCombine)
library(tmvtnorm)
#library(sarima)
library(coda)
library(sarima)

start_time <- Sys.time()

args = commandArgs(trailingOnly=TRUE)
d.raw <- as.double(args[1])
d <- 0.3#(d.raw - 1)/100
phi.raw <- as.double(args[2])
phi <- 0.5#(phi.raw - 1)/100
theta.raw <- as.double(args[3])
theta <- 0.5#(theta.raw - 1)/100
max_sims <- 30000#as.integer(args[4])

n_chain_start <- 1#as.integer(args[6])
n_chain_end <- 20#as.integer(args[5])

true_pars <- list(d=d,phi=phi,theta=theta,mu=0,sig2=1)
n <- 1000

gen_plots <- function(mcmc, par_lims=c()){
  n_pars <- ncol(mcmc)
  par(mfrow=c(n_pars,2))
  
  for (i in 1:n_pars){
    cat(i)
    if (length(par_lims) == 0){
      plot(as.vector(mcmc[,i]), type='l', main=paste('Trace',colnames(mcmc)[i]))
      plot(density(mcmc[,i]), main=paste('Density',colnames(mcmc)[i]))
      abline(v=mean(mcmc[,i]))
    }
    else{
      cat(par_lims[i,],'\n')
      plot(as.vector(mcmc[,i]), type='l', main=paste('Trace',colnames(mcmc)[i]),
           ylim=par_lims[i,])
      plot(density(as.vector(mcmc[,i])), main=paste('Density',colnames(mcmc)[i]),
           xlim=par_lims[i,])
      abline(v=mean(mcmc[,i]))
    }
  }
  par(mfrow=c(1,1))
}

#source('K:/PhD/2024/Honours Revisited/arfima_library.R')
gamma_ratio_est1 <- function(x, a, b=0){
  x^a
  #(x+(a+b-1)/2)^(a-b)
}
gamma_ratio_est2 <- function(x, a, b=0){
  #x^a
  (x+(a+b-1)/2)^(a-b)
}
indic_vec <- function(n, val=1){
  matrix(rep(val, n), ncol=1, nrow=n)
}
gauss_hypergeom_2nd <- function(a,b,c,x){
  1 + x * a*b/c + x^2 * a*(a+1)*b*(b+1)/(c*(c+1)*2)
}

acf_arma.1.1 <- function(n,phi,theta){
  acf <- rep(phi,n)
  acf[1] <- (1+2*phi*theta+theta^2)/(1-phi*phi)
  acf[2] <- (phi+theta)*(1+phi*theta)/(1-phi*phi)
  acf[2:n] <- cumprod(acf[2:n])
  return(acf)
}
acf_arfima.0.d.0 <- function(n, d){
  h <- seq(0,n-1)
  ret_arr <- numeric(length(h))
  ret_arr[1] <- gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * gamma(h[1] + d) / gamma(1+h[1] - d)
  ret_arr[2:n] = 
    gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * 1/gamma_ratio_est2(h[-1]+d, 1-2*d)
  return(ret_arr)
}
acf_arfima.0.d.0.true <- function(n, d){
  h <- seq(0,n-1)
  ret_arr <- numeric(length(h))
  ret_arr[1] <- gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * gamma(h[1] + d) / gamma(1+h[1] - d)
  ret_arr[2:n] = 
    gamma(1-2*d)/ (gamma(1-d) * gamma(d)) * 
    exp(lgamma(h[-1]+d) - lgamma(1+h[-1]-d))
  
  return(ret_arr)
}

acf_c_func <- function(lambda, d,h,rho,p=1){
  lambda[abs(h)+1] * 
    (rho^(2*p)*gauss_hypergeom_2nd(d+h,1,1-d+h,rho) + gauss_hypergeom_2nd(d-h,1,1-d-h,rho)-1)
}
acf_arfima.1.d.1 <- function(n, d, phi, theta){
  acf0.d.0_arr <- acf_arfima.0.d.0(n,d)
  h <- seq(0,n-1)
  c1 <- acf_c_func(acf0.d.0_arr, d,-h,-phi)
  c2 <- acf_c_func(acf0.d.0_arr, d,1-h,-phi)
  c3 <- acf_c_func(acf0.d.0_arr, d,2-h,-phi)
  return(
    (theta*c1+(1+theta^2)*c2+theta*c3)/(phi*(phi^2-1))
  )
}

post_psi0_arfima0d0_cond_y <- function(Y, d, log=T){
  n <- length(Y)
  acf_vec <- acf_arfima.0.d.0.true(n,d)
  
  cov_mat <- toeplitz(acf_vec)
  Sig_inv <- tinv(cov_mat)
  log_det_covmat <- determinant(cov_mat, logarithm=T)$modulus
  
  one_vec <- indic_vec(n)
  
  onev_siginv <- t(one_vec)%*%Sig_inv
  onev_siginv_onev <- onev_siginv %*% one_vec
  
  E <- t(Y)%*%Sig_inv%*%Y - t(Y)%*%Sig_inv%*%one_vec * onev_siginv_onev^(-1) * onev_siginv%*%Y
  
  if (E < 0){
    return(-Inf)
  }
  return(
    -1/2 * log_det_covmat - 
      1/2*log(abs(onev_siginv_onev)) -
      (n-1)/2 * log(E)
  )
}
post_psi0_arma_cond_y <- function(Y, phi, theta, log=T){
  n <- length(Y)
  acf_vec <- acf_arma.1.1(n, phi,theta) #ARMAacf(ar=phi,ma=theta,lag.max=length(Y)-1)
  
  cov_mat <- toeplitz(acf_vec)
  #note: R struggles with verrrry small numbers, which occur here (e-150+)
  #round these to zero beyond roughly 8b precision
  cov_mat[cov_mat < 1e-10] <- 0
  Sig_inv <- tinv(cov_mat)
  log_det_covmat <- determinant(cov_mat, logarithm=T)$modulus
  
  one_vec <- indic_vec(n)
  
  onev_siginv <- t(one_vec)%*%Sig_inv
  onev_siginv_onev <- onev_siginv %*% one_vec
  
  E <- t(Y)%*%Sig_inv%*%Y - t(Y)%*%Sig_inv%*%one_vec * onev_siginv_onev^(-1) * onev_siginv%*%Y
  
  if (E < 0){
    return(-Inf)
  }
  return(
    -1/2 * log_det_covmat - 
      1/2*log(abs(onev_siginv_onev)) -
      (n-1)/2 * log(E)
  )
}
post_psi0_cond_y <- function(Y, d, phi, theta, log=T){
  n <- length(Y)
  #acf_vec <- acf_arfima.1.d.1.true(n,d,phi,theta)
  acf_vec <- tacvfARFIMA(phi=phi,theta=-theta,dfrac=d,maxlag=n-1)
  
  cov_mat <- toeplitz(acf_vec)
  Sig_inv <- tinv(cov_mat)
  log_det_covmat <- determinant(cov_mat, logarithm=T)$modulus
  #cov_mat <- as.Toeplitz(acf_vec)
  #Sig_inv <- cov_mat$solve(method='pcg')
  #log_det_covmat <- cov_mat$log_det()
  
  one_vec <- indic_vec(n)
  
  onev_siginv <- t(one_vec)%*%Sig_inv
  onev_siginv_onev <- onev_siginv %*% one_vec
  
  E <- t(Y)%*%Sig_inv%*%Y - t(Y)%*%Sig_inv%*%one_vec * onev_siginv_onev^(-1) * onev_siginv%*%Y
  
  if (E < 0){
    return(-Inf)
  }
  return(
    -1/2 * log_det_covmat - 
      1/2*log(abs(onev_siginv_onev)) -
      (n-1)/2 * log(E)
  )
}

  Y <- NileMin
  #fix_Y_idx <- 1
  #fix_Y_idx <- 3
  #Y <- YDF$Y[(1+(fix_Y_idx-1)*n):(fix_Y_idx*n)]
  
  par_sims <- matrix(rep(0,(length(true_pars)+5)*max_sims), ncol=length(true_pars)+5, nrow=max_sims,
                     dimnames=list(rep(0,max_sims),c( 
                       "d","phi","theta",
                       "prop_d","prop_phi","prop_theta",
                       "posterior_d","alpha_d","posterior_arma","alpha_arma")))
  
  par_count <- 0
  while(T){
    flag = F
    par0 = tryCatch(
      {
        par0 <- mle.arfima(Y, p=1, q=1, d.range=c(0,0.5), incl.mean=F)
        post_tmp <- post_psi0_cond_y(Y, par0$pars['d'], par0$pars['ar.1'],par0$pars['ma.1'])
        #cat(post_tmp, ',', par0$pars)
        if (is.na(post_tmp) || is.infinite(post_tmp)){
          stop()
        }
        par0
      }, 
      error=function(e){print(e);flag=T},
      finally={}
    )
    if (!is.logical(par0)){break}
    par_count <- par_count + 1
  }
  
  mle_arfima1d1 <- as.list(par0$pars[c('d','ar.1','ma.1')])
  names(mle_arfima1d1) <- c('d','phi','theta')
  
  #two chain analysis - start with the arma chain
  U <- diffseries(Y, d=mle_arfima1d1$d)[-1] #diffseries just includes the first point
  arma_fit <- arima(U, order=c(1,0,1))
  arma_mle <- as.list(arma_fit$coef[c('ar1','ma1')])
  names(arma_mle) <- c('phi','theta')
  
  arma_sigm <- arma_fit$var.coef[1:2,1:2]
  
  par_sims[1,c('phi','theta')] <- rtmvnorm(1, mean=unlist(arma_mle), sigma=arma_sigm,
                                           lower=c(-1,-1),upper=c(1,1))
  par_sims[1, 'posterior_arma'] <- post_psi0_arma_cond_y(Y=U, 
                                                         phi=par_sims[1,'phi'], 
                                                         theta=par_sims[1,'theta'])
  prev_accepted_post_arma <- par_sims[1, 'posterior_arma']
  
  num_accept_arma <- 1
  
  #setup the d chain
  Z <- xarmaFilter(x=Y, model=list(ar=mle_arfima1d1$phi, ma=mle_arfima1d1$theta), whiten=T)
  d_mle <- list(d=mle_arfima1d1$d)
  
  #run 1 had 0.03
  #d_sig2 <- 0.03 #tweak parameter?
  d_sig2 <- 0.015
  
  par_sims[1,'d'] <- rtruncnorm(1, mean=d_mle, sd=d_sig2, a=0, b=0.5)
  par_sims[1, 'posterior_d'] <- post_psi0_arfima0d0_cond_y(Z, par_sims[1,'d'])
  prev_accepted_post_d <- par_sims[1, 'posterior_d']
  
  num_accept_d <- 1
  prev_num_accept_d <- 0
  
  curr_sim <- 2
  while(curr_sim <= max_sims){
    cat('Sim arma chain: ',curr_sim, '  Num-Accept - ARMA:', num_accept_arma,' d: ',num_accept_d,'           \r')
    
    prop_arma_mu <- par_sims[curr_sim-1,c('phi','theta')]
    prop_arma <- rtmvnorm(1, mean=unlist(prop_arma_mu), sigma=arma_sigm,
                          lower=c(-1,-1),upper=c(1,1))
    prop_d_mu <- par_sims[curr_sim-1,'d']
    prop_d <- rtruncnorm(1, mean=prop_d_mu, sd=d_sig2, a=0,b=0.5)
    
    #generate z~arfima(0,d,0) from the last accepted (phi,theta)
    Z <- xarmaFilter(x=Y, model=list(ar=par_sims[curr_sim-1,'phi'],ma=par_sims[curr_sim-1,'theta']), whiten=T)
    
    #generate u~arma(1,1) from the last accepted d
    U <- diffseries(Y, d=par_sims[curr_sim-1,'d'])[-1] #diffseries just includes the first point
    
    post_arma_prop <- post_psi0_arma_cond_y(U, prop_arma[1],prop_arma[2])
    post_arma_prop2 <- post_psi0_arma_cond_y(U, par_sims[curr_sim-1,'phi'],par_sims[curr_sim-1,'theta'])
    
    post_d_prop <- post_psi0_arfima0d0_cond_y(Z, prop_d)
    post_d_prop2 <- post_psi0_arfima0d0_cond_y(Z, par_sims[curr_sim-1,'d'])
    
    #determine the acceptance of the ARMA parameters (phi,theta)
    alpha_arma <- 0
    if (!is.na(post_arma_prop) && !is.infinite(post_arma_prop)){
      dens_prop1 <- dtmvnorm(prop_arma, mean=prop_arma_mu, sigma=arma_sigm,
                             lower=c(-1,-1),upper=c(1,1))
      
      dens_prop2 <- dtmvnorm(prop_arma_mu, mean=as.vector(prop_arma), sigma=arma_sigm,
                             lower=c(-1,-1),upper=c(1,1))
      alpha_s <- post_arma_prop + dens_prop1 - post_arma_prop2 - dens_prop2
      alpha_arma <- min(1, exp(alpha_s))
    }
    w_arma <- rbinom(1,1,alpha_arma)
    if (w_arma == 1){
      par_sims[curr_sim, c('phi','theta')] <- prop_arma
      prev_accepted_post_arma <- post_arma_prop
      num_accept_arma <- num_accept_arma + 1
    } else {
      par_sims[curr_sim, c('phi','theta')] <- par_sims[curr_sim-1, c('phi','theta')]
    }
    
    #determine the acceptance of the long-memory parameter d
    alpha_d <- 0
    if (!is.na(post_d_prop) && !is.infinite(post_d_prop)){
      dens_prop1 <- dtruncnorm(prop_d, mean=prop_d_mu, sd=d_sig2, a=0,b=0.5)
      dens_prop2 <- dtruncnorm(prop_d_mu, mean=prop_d, sd=d_sig2, a=0,b=0.5)
      
      alpha_s <- post_d_prop + dens_prop1 - post_d_prop2 - dens_prop2
      alpha_d <- min(1, exp(alpha_s))
    }
    
    w_d <- rbinom(1,1,alpha_d)
    if (w_d == 1){
      par_sims[curr_sim, 'd'] <- prop_d
      prev_accepted_post <- post_d_prop
      num_accept_d <- num_accept_d + 1
    } else {
      par_sims[curr_sim, 'd'] <- par_sims[curr_sim-1, 'd']
    }
    
    par_sims[curr_sim, 'prop_d'] <- prop_d
    par_sims[curr_sim, 'posterior_d'] <- post_d_prop
    par_sims[curr_sim, 'alpha_d'] <- alpha_d
    
    par_sims[curr_sim, c("prop_phi","prop_theta")] <- prop_arma
    par_sims[curr_sim, 'posterior_arma'] <- post_arma_prop
    par_sims[curr_sim, 'alpha_arma'] <- alpha_arma
    
    curr_sim <- curr_sim + 1
  }
#  print(chain)
#  d_mat[chain,] <- c(quantile(par_sims[(0.2*max_sims):max_sims,'prop_d'],prob=0.025),
#                     mean(par_sims[(0.2*max_sims):max_sims,'prop_d']),
#                     quantile(par_sims[(0.2*max_sims):max_sims,'prop_d'],prob=0.975))
  #  save_details(location)
#  ESS[chain] <- effectiveSize(par_sims[,'prop_d'])
  
save.image("Nile.RData")
  
end_time <- Sys.time()
time=(end_time-start_time)

plot(par_sims[,"prop_d"],type="l")
mean(par_sims[6001:30000,"prop_d"])
quantile(par_sims[6001:30000,"prop_d"],probs=c(0.025,0.975))
save.image("Nile.RData")

load("Nile.RData")

plot(NileMin,col="cadetblue",type="l",ylab="Annual Water Levels")
acf(NileMin,main="")

plot(obj$beta[,2])
plot(density(obj$beta[17000:24000,2]),col="cadetblue",xlab="d",main="",ylab="Posterior for d",lwd=2)

### qqplot

load("Nile.RData")
str(obj)
plot(obj$beta[,2])
qqnorm(obj$beta[,2],main="",cex=0.5)
qqline(obj$beta[,2],col="cadetblue",lwd=2)
