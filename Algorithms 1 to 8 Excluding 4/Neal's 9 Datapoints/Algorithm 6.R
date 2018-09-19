#Algo 6
data <- c(-1.48,-1.40,-1.16,-1.08,-1.02,0.14,0.51,0.53,0.78)
N<-9;d<-1
count1 <- 0
tot1 <- 0
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

res <- matrix(NA, ncol=9,nrow=(nsim-B))
param <-  rep(NA, N) #unique params 
for(u in 1:9){ 
  param[u]<-N.update(m0,s0,s,data[i])
}

ks <- c()
t1 <- c()
R<-4
###################### MCMC Iterations
for(t in 1:nsim){
  print(t)
  for(i in 1:N){
    for(r in 1:R){
      #proposal
      x <- (1:(N+1))[-i]
      prob <- c(rep(1,(N-1)),alpha)
      prob <- prob/sum(prob)
      
      prop <- sample(x, size=1, prob=prob)
      
      if(prop ==N+1){
        new.p <- 1
        tot1 <- tot1+1
        prop.par <- base.draw(m0, s0 )
      }else{
        prop.par <- param[prop]
        new.p <- 0
      }
      
      #MH step
      lik.new <- dnorm(data[i],mean=prop.par,sd=s)
      lik.old <- dnorm(data[i],mean=param[i],sd=s)
      
      if(lik.old==0){
        acc.rat <- 1
      }else{
        acc.rate <- lik.new/lik.old
      }
      
      if(runif(1) < acc.rate){
        param[i] <- prop.par
        count1 <- count1+1*new.p
      }

    }
  }
    
  #record
  if(t>B){
    ks <- c(ks, length(unique(param)))
    t1 <- c(t1,param[1])
    for(i in 1:9){
      res[t-B,i]<- param[i]
    }
  } 
}

auto.c <- function(x) length(x)/(effectiveSize(mcmc(x)))
auto.c(ks)

for(i in 1:9) print(auto.c(res[,i]))


plot(data)
fit6<-sapply(1:9, function(x) mean(res[,x]))
points(fit6, ,pch="x", col="red")
legend("topleft", legend=c("Neal's 9 Data-points","Mean of Recorded  Parameters")
       ,pch=c("o","x"), col=1:2)
count1/tot1
