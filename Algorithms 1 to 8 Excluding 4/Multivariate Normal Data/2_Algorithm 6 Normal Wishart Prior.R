#Algo 6


#load packages
#---------------------------------------------------
install.packages("mvtnorm")
library(mvtnorm)
#---------------------------------------------------

#create data
#---------------------------------------------------
#number of data points
N<-3*10

#data
tm1 <- rmvnorm(1,mean = c(10,-10))
tm2<- rmvnorm(1,mean = c(0,0))
tm3 <- rmvnorm(1,mean = c(10,10))

tm <- rbind(tm1,tm2,tm3)
sig <- matrix(c(1,0.5,0.5,1),2,2)
y <- matrix(NA, N,3) #data
for(i in 1:3){
  y[(1+(i-1)*N/3):(N/3*i),1:2] <- rmvnorm(N/3,mean = tm[i,],sigma =sig)
  y[(1+(i-1)*N/3):(N/3*i),3]<-i
}
data<-y[,1:2]
#---------------------------------------------------

#NW updates
#---------------------------------------------------
NW.update <- function(m0,b0,a0,B0,data){
  
  if(is.vector(data)){
    data <- t(as.matrix(data))
  }else{
    data<-as.matrix(data)
  }
  
  #num of datapoints, dimension
  N<-dim(data)[1]
  d<-dim(data)[2]
  
  #mean of data
  x.bar <- apply(data,2,sum)/N
  
  #sig
  xmm<-matrix(0,d,d)
  for(i in 1:N){
    x<-as.matrix(data[i,])
    xmm<-xmm+(x-x.bar)%*%t(x-x.bar)
  }
  sig.bar<- (1/N)*xmm
  
  #degrees of freedom
  a.N<-a0+N/2
  #prior precision of mean relative to param
  b.N <- b0+N
  #mean
  m.N <-(b0*m0+N*x.bar)/b.N
  #scale matrix  
  B.N <- B0+(N/2)*(sig.bar+(b0/b.N)*(x.bar-m0)%*%t(x.bar-m0))
  
  #Posterior draws
  w<-rWishart(1,a.N,B.N)
  w<-matrix(w,d,d)
  m<-rmvnorm(1,mean=m.N,sigma=b.N*solve(w))
  
  return(list(mu=m,sig=solve(w)))
}  
#---------------------------------------------------

# Prior draw
#---------------------------------------------------
base.draw<-function(m0,b0,a0,B0){
  w<-rWishart(1,a0,B0)
  w<-matrix(w,d,d)
  m<-rmvnorm(1,mean=m0,sigma=b0*solve(w))
  return(list(mu=m,sig=solve(w)))
}

#---------------------------------------------------

#initiate params
#---------------------------------------------------
#prior param
m0<-c(0,0)
b0<-1
a0<-2
B0<-matrix(c(10,5,5,10),2,2)
d<-length(data[1,])
N<-length(data[,1])
alpha<-1
R <- 4

#unique params
param <-  list() 
for(s in 1:N){ 
  param[[s]]<-base.draw(m0 ,b0 ,a0  ,B0  )
}
#---------------------------------------------------

#mcmc
#---------------------------------------------------
#Num of iterations
warm.up <-100
iter<-10000+warm.up

ks<- c()
t1 <- c()

#---------------------------------------------------

#MCMC Iterations
#---------------------------------------------------
for(t in 1:iter){
  print(t)
  for(i in 1:N){
    for(r in 1:R){
      #proposal
      x <- (1:(N+1))[-i]
      prob <- c(rep(1,(N-1)),alpha)
      prob <- prob/sum(prob)
      
      prop <- sample(x, size=1, prob=prob)
      if(prop ==N+1){
        prop.par <- base.draw(m0, b0,a0, B0)
      }else{
        prop.par <- param[[prop]]
      }
      
      #MH step
      lik.new <- dmvnorm(data[i,],mean=prop.par$mu,sigma=prop.par$sig)
      lik.old <- dmvnorm(data[i,],mean=param[[i]]$mu,sigma=param[[i]]$sig)
      
      if(lik.old==0){
        acc.rat <- 1
      }else{
        acc.rate <- lik.new/lik.old
      }
      
      if(runif(1) < acc.rate){
        param[[i]] <- prop.par
      }   
    }
  }
  if(t>warm.up){
    ks <- c(ks, length(unique(clus)))
    t1<- rbind(t1, param[[clus[1]]]$mu)
  }
  
}
#---------------------------------------------------

#Analysis
#---------------------------------------------------
auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(ks)
auto.c(t1)
#---------------------------------------------------
