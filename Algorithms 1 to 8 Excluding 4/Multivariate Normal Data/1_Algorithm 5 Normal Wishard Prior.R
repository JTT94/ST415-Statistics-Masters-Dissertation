#algo 5 self contained

#create data

#number of data points
#---------------------------------------------------
N<-3*10
#---------------------------------------------------

#data
#---------------------------------------------------
tm1 <- rmvnorm(1,mean = c(10,-10))
tm2 <- rmvnorm(1,mean = c(0,0))
tm3 <- rmvnorm(1,mean = c(10,10))

tm <- rbind(tm1,tm2,tm3)

sig <- matrix(c(1,0.5,0.5,1),2,2)

y <- matrix(NA, N,3)
for(i in 1:3){
  y[(1+(i-1)*N/3):(N/3*i),1:2] <- rmvnorm(N/3,mean = tm[i,],sigma =sig)
  y[(1+(i-1)*N/3):(N/3*i),3]<-i
}
data<-y[,1:2]
#---------------------------------------------------

#NW updates, sample from posterior
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
  
  #sample from posterior
  w<-rWishart(1,a.N,B.N)
  w<-matrix(w,d,d)
  m<-rmvnorm(1,mean=m.N,sigma=b.N*solve(w))
  
  return(list(mu=m,sig=solve(w)))
  
}  
#---------------------------------------------------

#Prior draw
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
B0<-matrix(c(1,0.5,0.5,1),2,2)
d<-length(data[1,])
N<-length(data[,1])
alpha<-1
R <- 4

#initiate clusters
clus <- c(rep(1,10),rep(2,10),rep(3,10))
u.clus <- sort(unique(clus)) #unique clusters
param <-  list() #unique params 

for(s in u.clus){ 
  param[[s]]<-NW.update(m0,b0,a0,B0,data[which(clus==s),])
}

#---------------------------------------------------

#mcmc
#---------------------------------------------------

#Num of iterations
warm.up <- 100
iter<-20000+warm.up

ks<- c()
t1 <- c()

#number of data points in each cluster
num<-list()
for(i in u.clus){
  num[[i]]<-length(which(clus==i))
  
}

#MCMC Iterations

for(t in 1:iter){
  print(t)
  #update allocations
  for(i in 1:N){
    for( r in 1:R){
      u.clus <- sort(unique(clus))
      #number of data points in each cluster
      num <- sapply(u.clus, function(x) length(which(clus[-i]== x)))
      
      #conditional prior draw
      prob <- c(num, alpha)
      prob <- prob/sum(prob)
      n.clus <- (max(u.clus)+1)
      
      prop <- sample(c(u.clus,n.clus),size=1, prob=prob)
      
      if (prop==n.clus){
        param[[prop]] <- base.draw(m0,b0,a0,B0)
      }
      
      #MH step
      lik.old <- dmvnorm(data[i,],mean=param[[clus[i]]]$mu,sigma=param[[clus[i]]]$sig)
      lik.new <- dmvnorm(data[i,],mean=param[[prop]]$mu,sigma=param[[prop]]$sig)
      
      if(lik.old==0){
        acc.rat <- 1
      }else{
        acc.rate <- lik.new/lik.old
      }
      
      if(runif(1) < acc.rate) clus[i] <- prop
    }
  }
  
  #print("Allocations Done")
  
  #update parameters
  u.clus <- sort(unique(clus))
  for(s in u.clus){ 
    param[[s]]<-NW.update(m0,b0,a0,B0,data[which(clus==s),])
  }
  
  u.clus <- sort(unique(clus))
  t.param <- list()
  for(i in 1:length(u.clus)){
    clus[which(clus==u.clus[i])] <- i
    t.param[[i]] <- param[[u.clus[i]]] 
  }
  param <- t.param

 #Record results
  if(t>warm.up){
    ks <- c(ks, length(u.clus))
    t1<- rbind(t1, param[[clus[1]]]$mu)
  }
}
#---------------------------------------------------

#Analysis
#---------------------------------------------------
#Look at results
auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(ks)
auto.c(t1)
#---------------------------------------------------
