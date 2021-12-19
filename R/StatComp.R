#' @title EM_mix
#' @description EM Algorithm for two-component Gaussian Mixture
#' @param mu1 the mean of Y1
#' @param sigma1 the variance of Y1
#' @param mu2 the mean of Y2
#' @param sigma2 the variance of Y2
#' @param p the probability of the mixture distribution
#' @param N the number of samples
#' @param rv the random number from mixed normal distribution
#' @import stats
#' @import MASS
#' @import bootstrap
#' @import boot
#' @import energy
#' @import Ball 
#' @import RANN
#' @import microbenchmark
#' @return the maximum likelihood estimation of the mean, variance and the mixed probablity
#' @examples
#' \dontrun{
#' rv_mix<-function(N,p){
#' u1<-rnorm(N,0,1)
#' u2<-rnorm(N,2,1)
#' j<-runif(N)
#' k<-as.integer(j>p)
#' u<-k*u1+(1-k)u2
#' return (u)
#' }
#' 
#' p<-0.6
#' N<-1000
#' set.seed(2335)
#' rv<-rv_mix(N,p)
#' p_0<-0.5
#' mu1<-0.5
#' mu2<-1.5
#' sigma1<-sigma2<-sd(rv)
#' EM_mix(mu1,sigma1,mu2,sigma2,p_0,N,rv)
#' }
#' @export
EM_mix<-function(mu1,sigma1,mu2,sigma2,p,N,rv){
  phi<-rep(0,N)
  for (i in 1:1000) {
    #E-step
    for (j in 1:N) {
      f1<-dnorm(rv[j],mu1,sigma1)
      f2<-dnorm(rv[j],mu2,sigma2)
      phi[j]<-p*f2/((1-p)*f1+p*f2)
    }
    #M-step
    mu1<-sum((1-phi)*rv)/sum(1-phi)
    s1<-sum((1-phi)*(rv-mu1)^2)/sum(1-phi)
    mu2<-sum(phi*rv)/sum(phi)
    s2<-sum(phi*(rv-mu2)^2)/sum(phi)
    sigma1<-sqrt(s1)
    sigma2<-sqrt(s2)
    p<-sum(phi)/N
  }
  return(c(mu1,sigma1,mu2,sigma2,p))
}
rv_mix<-function(N,p){
  u1<-rnorm(N,0,1)
  u2<-rnorm(N,2,1)
  j<-runif(N)
  k<-as.integer(j>p)
  u<-k*u1+(1-k)*u2
  return (u)
}
#' @title Gibbs_mix
#' @description Gibbs sampling for two-component Gaussian Mixture
#' @param mu1 the mean of Y1
#' @param sigma1 the EM-estimated variance of Y1
#' @param mu2 the mean of Y2
#' @param sigma2 the EM-estimated variance of Y2
#' @param p the probability of the mixture distribution
#' @param N the number of samples
#' @param rv the random number from mixed normal distribution
#' @import stats
#' @return the maximum likelihood estimation of the mean, variance and the mixed probablity
#' @examples
#' \dontrun{
#' rv_mix<-function(N,p){
#' u1<-rnorm(N,0,1)
#' u2<-rnorm(N,2,1)
#' j<-runif(N)
#' k<-as.integer(j>p)
#' u<-k*u1+(1-k)u2
#' return (u)
#' }
#' 
#' p<-0.6
#' N<-1000
#' set.seed(2335)
#' rv<-rv_mix(N,p)
#' p_0<-0.5
#' mu1<-0.5
#' mu2<-1.5
#' sigma1<-sigma2<-sd(rv)
#' EM<-EM_mix(mu1,sigma1,mu2,sigma2,p_0,N,rv)
#' sigma1<-EM[2]
#' sigma2<-EM[4]
#' p_hat<-EM[5]
#' Gibbs_mix(mu1,mu2,sigma1,sigma2,p_hat,N,rv)
#' }
#' @export
Gibbs_mix<-function(mu1,mu2,sigma1,sigma2,p,N,rv){
  X<-matrix(0,1000,2)
  z<-rep(0,N)
  for (i in 1:1000) {
    for (j in 1:N){
      f1<-dnorm(rv[j],mu1,sigma1)
      f2<-dnorm(rv[j],mu2,sigma2)
      phi<-p*f2/((1-p)*f1+p*f2)
      
      u<-runif(1)
      if(u<phi){
        z[j]<-1
      }else{
        z[j]<-0
      }
    }
    mu1<-sum((1-z)*rv)/sum(1-z)
    mu2<-sum(z*rv)/sum(z)
    X[i,1]<-rnorm(1,mu1,sigma1)
    X[i,2]<-rnorm(1,mu2,sigma2)
  }
  return(X)
}
rv_mix<-function(N,p){
  u1<-rnorm(N,0,1)
  u2<-rnorm(N,2,1)
  j<-runif(N)
  k<-as.integer(j>p)
  u<-k*u1+(1-k)*u2
  return (u)
}