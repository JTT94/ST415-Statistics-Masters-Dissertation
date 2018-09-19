


N.update <- function(m0,s0,s,data){
  #num of datapoints, dimension
  N<-length(data)
  #mean of data
  x.bar <- mean(data)
  p.m <- (s0^2*x.bar + s^2*m0/N)/(s^2/N+s0^2)
  p.s <- 1/(1/s0^2 +N/s^2)
  m <- rnorm(1,p.m,sqrt(p.s))
  return(m)
} 

base.draw<-function(m0,s0){
  m <- rnorm(1,m0,s0)  
  return(m)
}
#####################

#Initiate
#known variance
#prior param
m0<-0; s0<-1; alpha<-1
s<-0.1



#####################
##methods
#conjugate normal update for mean  
#initiate clusters and component params
clus <-  sample(1:9,size=9, replace=T)
nsim <- 20000; B<-100

param <-  rep(NA, N) #unique params 
n <- N.update(m0,s0,s,data[i])
for(u in 1:9){ 
  param[u]<-n
}


ks <- c()
t1 <- c()
###################### MCMC Iterations
q0 <- function(x,s0,s,m0) dnorm(x,mean=m0,sd=sqrt(s0^2+s^2))



#initial calculations
#marginal Integral G0 * F(xi|param) vector
q <- sapply(data, function(x) q0(x,s0,s,m0))
#post.m <- t(apply(y,1, function(x) m.update(m0,x)))
####################

##Gibbs iterations
draw <- numeric( N)

for( t in 1:nsim){
  
  print(t)
  for(i in 1:N){
    
    
    draw[1] <- N.update(m0 ,s0 ,s ,data[i])
    draw[-1] <- param[i]
    
    lik <- sapply(param[-i], function(x) dnorm(data[i], mean = x, sd = s) )
    prob <- c(alpha*q[i],lik)
    prob <- prob/sum(prob)
    #update parameter for datapoint i
    param[i] <- sample(draw,size=1,prob=prob)
  } 

  #record
  if(t>B){
    t1<- c(t1, param[1])
    ks <- c(ks, length(unique(param)))
  }
}

#results
auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(t1)
auto.c(ks)
plot(t1,type="l")
