#install.packages("coda")
#library(coda)
#Neal Data
data <- c(-1.48,-1.40,-1.16,-1.08,-1.02,0.14,0.51,0.53,0.78)

N<-9;d<-1



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
clus <- sample(1:9,size=9, replace=T)
nsim <- 20000; B<-100
u.clus <- sort(unique(clus)) #unique clusters
param <-  list() #unique params 
for(cs in u.clus){ 
  param[[cs]]<-N.update(m0,s0,s,data[which(clus==cs)])
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

#mcmc
for( t in 1:nsim){
  print(t)
  for(i in 1:N){

    #draw new cluster allocation
    #likelihood  component of probability of data point i being in each of the unique clusters
    lik <- sapply(u.clus, function(x) dnorm(data[i], mean = param[[x]], sd = s) )
    
    #number of data points in each cluster
    num <- sapply(u.clus, function(x) length(which(clus[-i]== x)))
                  #as.numeric (table(clus[-1]) ) will do this job too
  
    #prob of new cluster, up to constant of proportioinality
    p.new <- alpha*q[i]
    
    #prob of existing cluster, this is a vector, up to multiplicative constant
    p.prev <- num*lik
  
    #draw new cluster
    prob <- c(p.prev,p.new)/ sum(c(p.new,p.prev))
    n.clus <- (max(u.clus)+1)
    clus[i] <- sample(c(u.clus,n.clus), size = 1,prob = prob)
    u.clus <- unique(clus)
    if(clus[i] == n.clus){
      param[[clus[i]]] <- N.update(m0,s0,s,data[which(clus==clus[i])])
    }
  }
  
  #draw parameters for clusters
  #update parameters
  u.clus <- sort(unique(clus))
  for(u in u.clus){ 
    param[[u]]<-N.update(m0,s0,s,data[which(clus==u)])
  }
  
  #renumber clusters
  
  t.param <- list()
  for(i in 1:length(u.clus)){
    clus[which(clus==u.clus[i])] <- i
    t.param[[i]] <- param[[u.clus[i]]] 
  }
  param <- t.param
  
  
  #record
  if(t>B){
    t1<- c(t1, param[[clus[1]]])
    ks <- c(ks, length(u.clus))
  }
  u.clus <- sort(unique(clus))
}

auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(ks)
auto.c(t1)

##################

