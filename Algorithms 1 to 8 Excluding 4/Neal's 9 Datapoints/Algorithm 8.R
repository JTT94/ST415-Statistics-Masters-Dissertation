#Neal data
#install.packages("coda")
#library(coda)

#Neal Data
data <- c(-1.48,-1.40,-1.16,-1.08,-1.02,0.14,0.51,0.53,0.78)

N<-9;d<-1

#Initiate
#known variance
#prior param
m0<-0; s0<-1; alpha<-1
s<-0.1
#normal updates, known variance
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

#initiate clusters and component params
clus <- rep(1,9)
u.clus <- sort(unique(clus)) #unique clusters
param <-  list() #unique params 

for(cs in u.clus){ 
  param[[cs]]<-N.update(m0,s0,s,data[which(clus==cs)])
}

#algorithm 8 
#algo 8
######################
#initiate
nsim <- 20000; m<- 1; B <- 100
####################
# #mcmc
res <- matrix(NA, ncol=9,nrow=(nsim-B))
ks <- c()
t1 <- c()
#mcmc
res1<- c()
res2 <-c()

  ks <- c()
  t1 <- c()
for( t in 1:nsim){
  print(t)
  for(i in 1:N){
    #compute likelihoods
    lik <- c()
    
    #existing clusters
    u.clus <- sort(unique(clus))
    
    for(u in u.clus){
      num<-length(which(clus[-i]==u))
      tmp <- num*dnorm(data[i],mean=param[[u]],sd=(s))
      lik <- c(lik, tmp )
    }
    
    #draw new parameters for m
    m.phi<-list()
    for(r in 1:m){
      m.phi[[r]]<- base.draw(m0,s0)
      tmp <-alpha*dnorm(data[i],mean=m.phi[[r]],sd=(s))/m  
      lik <- c(lik,   tmp) 
    }
    
    #normalise
    p <-lik/sum(lik)
    
    #sample component allocation
    ind <- c(u.clus,(1:m))
    draw <- sample(1:length(ind),size=1,replace = F,prob =p)
    
    if(draw > length(u.clus)){
      clus[i]<- max(u.clus)+1
      param[[clus[i]]] <- m.phi[[ind[draw]]]
    }else{
      clus[i]<- ind[draw]
    }
    
  }
  
  #update parameters
  #existing clusters
  u.clus <- sort(unique(clus))
  for(u in u.clus){ 
    param[[u]]<-N.update(m0,s0,s,data[which(clus==u)])
  }

  u.clus <- sort(unique(clus))
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
    
    res[t-B,] <- sapply(1:N, function(i) param[[clus[i]]])
  }
  
}

##########################
#empirical mean of results




auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))

plot(data)

fit8<-sapply(1:9, function(x) mean(res[,x]))
points(fit, ,pch="x", col="red")


