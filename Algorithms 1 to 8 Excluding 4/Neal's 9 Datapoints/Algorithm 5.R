
#Neal Data
data <- c(-1.48,-1.40,-1.16,-1.08,-1.02,0.14,0.51,0.53,0.78)

N<-9;d<-1

count1<-0
tot1<-0


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
#initiate clusters and component params
clus <- sample(1:9, size=9, replace=T)
u.clus <- sort(unique(clus)) #unique clusters
param <-  list() #unique params 

for(cs in u.clus){ 
  param[[cs]]<-N.update(m0,s0,s,data[which(clus==cs)])
}
count1/tot1
ks <- c()
t1 <- c()
###################### MCMC Iterations
iter <- 20000; R<-4;B <- 100
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
      new.p <- 0
      if (prop==n.clus){
        param[[prop]] <- base.draw(m0,s0)
        tot1<-tot1+1
        new.p <-1
      }
      
      #MH step
      lik.old <- dnorm(data[i],mean=param[[clus[i]]],sd=s)
      lik.new <- dnorm(data[i],mean=param[[prop]],sd=s)
      
      if(lik.old==0){
        acc.rat <- 1
      }else{
        acc.rate <- lik.new/lik.old
      }
      
      if(runif(1) < acc.rate){
        clus[i] <- prop
        count1<- count1+1*new.p
      } 
    }
  }
  
  #print("Allocations Done")
  #update parameters
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
  }
  
  
}

count1/tot1
auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(ks)
auto.c(t1)
mean(t1)
