#algo 7 

#initiate params
#---------------------------------------------------
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

#mcmc
#---------------------------------------------------

#Num of iterations
warm.up <-100
iter<-10000+warm.up

res <- list() #store params at each iteration
comp<-list() #store clus allocation at each iteration

#number of data points in each cluster
num<-list()
for(i in u.clus){
  num[[i]]<-length(which(clus==i))
  
}

#MCMC Iterations
#---------------------------------------------------

for(t in 1:iter){
  print(t)
  #proposal 1
  #---------------------------------------------------
  for(i in 1:N){
    
    #If non singleton
    if(num[[clus[i]]]>1){
      #Proposal draw from base dist
      phi<- base.draw(m0,b0,a0,B0)
      x1<-dmvnorm(data[i,],mean=phi$mu,sigma=phi$sig)
      y1<-dmvnorm(data[i,],mean=param[[clus[i]]]$mu,sigma=param[[clus[i]]]$sig)
      
      #Accept / Reject
      if(y1==0){
        acc.rat <- 1
      }else{
        acc.rat <- (alpha/(N-1))*x1/y1
      }
      
      if(runif(1)<acc.rat){
        clus[i]<- max(clus)+1
        param[[clus[i]]]<-phi
      }
      
    }else{
      #If singleton, draw from existing clusters
      allocation<-sample(sort(u.clus),size=1,replace=F,prob=table(clus)/(N-1))
      phi<- param[[allocation]]
      
      x1<-dmvnorm(data[i,],mean=phi$mu,sigma=phi$sig)
      y1<-dmvnorm(data[i,],mean=param[[clus[i]]]$mu,sigma=param[[clus[i]]]$sig)
      
      #Accept / Reject
      if(y1==0){
        acc.rat <- 1
      }else{
        acc.rat <- ((N-1)/alpha)*x1/y1
      }
      
      if(runif(1)<acc.rat){
        clus[i]<- allocation
      } 
    }
    u.clus<- sort(unique(clus))
    #number of data points in each cluster
    for(s in u.clus){
      num[[s]]<-length(which(clus==s))
      
    }
  }
  #---------------------------------------------------
  
  #proposal 2
  #---------------------------------------------------
  for(i in 1:N){
    
    # move non-singletons about
    if(num[[clus[i]]]>1){
      lik<- lapply(u.clus, function(x) dmvnorm(data[i,],mean=param[[x]]$mu,sigma=param[[x]]$sig))
      lik<-unlist(lik)
      p<- table(clus)*lik*alpha/(N-1)
      p<-p/sum(p)
      clus[i]<-sample(sort(u.clus),size=1,replace=F,prob=p)
    }
    
    u.clus<- sort(unique(clus))
    #number of data points in each cluster
    for(s in u.clus){
      num[[s]]<-length(which(clus==s))
    }
    
  }
  
  #---------------------------------------------------
  
  
  #update parameters
  #---------------------------------------------------
  for(s in u.clus){ 
    param[[s]]<-NW.update(m0,b0,a0,B0,data[which(clus==s),])
  }
  #---------------------------------------------------
  
  # plot groupings
  #---------------------------------------------------
  if(t%%100==0){
    plot(data)
    gr <- tapply(1:N,clus,function(x) x, simplify=F)
    for(j in 1:dim(gr)){
      po<- data[unlist(gr[j]),]
      if(is.vector(po)) po <- t(matrix(po))
      points(po,col=j)
    } 
  }
  #---------------------------------------------------
  
  #Record results
  #---------------------------------------------------
  if(t>warm.up){
    ks <- c(ks, length(unique(clus)))
    t1<- rbind(t1, param[[clus[1]]]$mu)
  }
  #---------------------------------------------------
  
}

#Analysis
#---------------------------------------------------
auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(ks)
auto.c(t1)
#---------------------------------------------------
