#algo 8 self contained
library(mvtnorm)
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
xs <- matrix(NA, N,3) #data
for(i in 1:3){
  xs[(1+(i-1)*N/3):(N/3*i),1:2] <- rmvnorm(N/3,mean = tm[i,],sigma =sig)
  xs[(1+(i-1)*N/3):(N/3*i),3]<-i
}
data <- xs[,1:2]
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
  
  #Posterior draw
  w<-rWishart(1,a.N,B.N)
  w<-matrix(w,d,d)
  m<-rmvnorm(1,mean=m.N,sigma=b.N*solve(w))
  
  return(list(mu=m,sig=solve(w)))
  
}  
#---------------------------------------------------

#Base draw
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

#initiate clusters
clus <- c(rep(1,10),rep(2,10),rep(3,10))
u.clus <- sort(unique(clus)) #unique clusters
param <-  list() #unique params 
ks<- c()
t1 <- c()

for(s in u.clus){ 
  param[[s]]<-NW.update(m0,b0,a0,B0,data[which(clus==s),])
}
#---------------------------------------------------

#algo 8
#---------------------------------------------------
#initiate
warm_up <- 100
nsim <- 10000+warm_up
m<- 3
res <- list() #store params at each iteration
comp<-list() #store clus allocation at each iteration


#mcmc
#---------------------------------------------------
for( t in 1:nsim){
  print(t)
  for(i in 1:N){
    #compute likelihoods
    lik <- c()
    
    #existing clusters
    u.clus <- sort(unique(clus))
    for(u in u.clus){
      num<-length(which(clus[-i]==u))
      tmp <- num*dmvnorm(data[i,],mean=param[[u]]$mu,sigma=param[[u]]$sig)
      lik <- c(lik, tmp )
      
    }
    
    #draw new parameters for m
    m.phi<-list()
    for(r in 1:m){
      m.phi[[r]]<- base.draw(m0,b0,a0,B0)
      tmp <-alpha*dmvnorm(data[i,],mean=m.phi[[r]]$mu,sigma=m.phi[[r]]$sig)/m
      lik <- c(lik,   tmp) 
    }
    
    #normalise
    p <-lik/sum(lik)
    
    #sample component allocation
    ind <- c(u.clus,(1:m))
    draw <- sample(1:length(ind),size=1,replace = F,prob =p)
    
    # If new component, add new component index and associated parameter
    if(draw > length(u.clus)){
      clus[i]<- max(u.clus)+1
      param[[clus[i]]] <- m.phi[[ind[draw]]]
      
    }else{
      # if existing cluster, set index and associated parameter to be the same
      clus[i]<- ind[draw]
    }
    
  }
  
  #update parameters
  #---------------------------------------------------
  #existing clusters
  #---------------------------------------------------
  u.clus <- sort(unique(clus))
  for(u in u.clus){ 
    param[[u]]<-NW.update(m0,b0,a0,B0,data[which(clus==u),])
  }
  #---------------------------------------------------
  
  #record
  #---------------------------------------------------
  if(t>warm_up){
    ks <- c(ks, length(unique(clus)))
    t1<- rbind(t1, param[[clus[1]]]$mu)
  }
}

#Analysis
#---------------------------------------------------
auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(ks)
auto.c(t1)
