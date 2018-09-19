

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
clus <- sample(1:9,size=9, replace=T)
u.clus <- sort(unique(clus)) #unique clusters
param <-  list() #unique params 

for(cs in u.clus){ 
  param[[cs]]<-N.update(m0,s0,s,data[which(clus==cs)])
}


######################
#initiate
nsim <- 20000; m<- 2; B <- 100
ks <- c()
t1 <- c()
for(u in u.clus){
  num[[u]]<-length(which(clus==u))
  
}

count1 <- 0
tot1 <- 0
count2 <- 0
tot2 <- 0

###################### MCMC Iterations

for(t in 1:nsim){
  print(t)
  #proposal 1
  for(i in 1:N){
    if(num[[clus[i]]]>1){
      phi<- base.draw(m0,s0)
      
      x1<-dnorm(data[i],mean=phi,sd=s)
      y1<-dnorm(data[i],mean=param[[clus[i]]],sd=s)
      
      if(y1==0){
        acc.rat <- 1
      }else{
        acc.rat <- (alpha/(N-1))*x1/y1
      }
      tot1 <- tot1+1
      if(runif(1)<acc.rat){
        count1 <- count1+1
        clus[i]<- max(clus)+1
        param[[clus[i]]]<-phi
      }
    }else{
      tot2 <- tot2+1
      allocation<-sample(sort(u.clus),size=1,replace=F,prob=table(clus)/(N-1))
     
      phi<- param[[allocation]]
      
      x1<-dnorm(data[i],mean=phi,sd=s)
      y1<-dnorm(data[i],mean=param[[clus[i]]],sd=s)
      
      
      if(y1==0){
        acc.rat <- 1
      }else{
        acc.rat <- ((N-1)/alpha)*x1/y1
      }
      
      
      if(runif(1)<acc.rat){
        count2 <- count2+1
        clus[i]<- allocation
      } 
    }
    u.clus<- sort(unique(clus))
    #number of data points in each cluster
    for(u in u.clus){
      num[[u]]<-length(which(clus==u))
      
    }
  }
  
  #proposal 2
  for(i in 1:N){
    
    if(num[[clus[i]]]>1){
      lik<- lapply(u.clus, function(x) dnorm(data[i],mean=param[[x]],sd=s))
      lik<-unlist(lik)
      p<- table(clus)*lik*alpha/(N-1)
      p<-p/sum(p)

      if(length(u.clus)==1){clus[i]  <- u.clus
                            }else{
                              clus[i]<-sample(sort(u.clus),size=1,replace=F,prob=p)
                            }

    }
    u.clus<- sort(unique(clus))
    #number of data points in each cluster
    for(u in u.clus){
      num[[u]]<-length(which(clus==u))
      
    }
  }
  
  
  
  #update parameters
  for(u in u.clus){ 
    param[[u]]<-N.update(m0 ,s0,s,data[which(clus==u)])
  }
  
  #record
  if(t>B){
    t1<- c(t1, param[[clus[1]]])
    ks <- c(ks, length(u.clus))
  }
  
  u.clus <- sort(unique(clus))
  t.param <- list()
  for(i in 1:length(u.clus)){
    clus[which(clus==u.clus[i])] <- i
    t.param[[i]] <- param[[u.clus[i]]] 
  }
  param <- t.param
  u.clus <- sort(unique(clus))
}

auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(t1)
auto.c(ks)
length(t1)
