#useful functions
######################
base.draw <-function(k,J,m0,s0,ms,vs,v,t){
  a <- rnorm(J,mean=0, sd=sqrt(t))
  b <- matrix(NA, ncol = J, nrow=k)
  for(i in 1:J) b[,i] <- rnorm(k,mean=0,sd=sqrt(v))
  
  mu <- rnorm(k, mean = m0, sd=sqrt(s0))
  l.sig <- rnorm(k, mean = ms, sd=sqrt(vs))
  sig <- exp(l.sig)
  
  return(list(sig=sig, a=a,b=b,mu=mu))
}
######################
#liklihood for class label
py <- function(x,y,param){
  a <- param$a[y]
  b <- param$b[,y]
  num <- exp(a+sum(x*b))
  den <- sum(sapply(1:J, function(i) exp(param$a[i]+sum(x*param$b[,i]))))
  #print(num)
  #print(den)
  res <- num/den
  return(res)
}
#likelihood of covariates
px <- function(x,param){
  mu <- param$mu
  sig <- param$sig
  
  px <- dmvnorm(x,mean=mu,sigma=diag(sig))
  return(px)
}
#likelihood 
f<-function(x,y,param){
  a <- param$a[y]
  b <- param$b[,y]
  
  mu <- param$mu
  sig <- param$sig
  
  px <- dmvnorm(x,mean=mu,sigma=diag(sig))
  
  num <- exp(a+sum(x*b))
  den <- sum(sapply(1:J, function(i) exp(param$a[i]+sum(x*param$b[,i]))))
  py <- num/den
  
  return(px*py)
} 


#classification
#classify a single datapoint based on covariates and results from mcmc
classify <- function(testx, res,comp){
  guess <- rep(NA,3)
  for(j in 1:J){
    #cluster covariate
    lR <- length(res)
    cls <- rep(NA, lR)
    for(s in 1:lR){
      u.clus <- sort(unique(comp[[s]]))
      al <- which.max(sapply(u.clus, function(u) px(x = testx, res[[s]][[u]] )))
      cls[s] <- u.clus[al]
    }
    
    guess[j] <- sum(sapply(1:lR, function(s) py(x = testx, y = j,param = res[[s]][[cls[s]]])))/lR
  }
  which.max(guess)
}

