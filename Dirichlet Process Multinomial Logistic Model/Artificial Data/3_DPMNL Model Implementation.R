#load data
#load cl which is correct labels
#---------------------------------------------------

#split into test and training data
#---------------------------------------------------
n <- dim(data)[1]
r.ind <- sample(1:n,size=n, replace=F)
#split 1:n into k equal groups
k.cuts <- cut(r.ind, breaks = 3, labels=F)
test.ind <- which(k.cuts==1)
train.ind <- which(k.cuts !=1)
#covariates xs
xs <- data[train.ind, ]
#correct labels
cl <- (cl[train.ind])

#---------------------------------------------------

#initialise
#---------------------------------------------------
alpha<-0.5

#number of datapoints
N <- dim(xs)[1]

#dimension of covariates
NX <-k <- dim(xs)[2]
m0 <- 0
s0 <- 1
ms <- 0
vs <- 2

#dimension of number of classes
J <- length(unique(cl))
v <- 1
t <- 1

#set initial clusters as class labels ones
clus<-cl
u.clus <- sort(unique(clus)) #unique clusters
#algo 8
m<- 3
R <-5
params <-  list() #unique params 
for(s in u.clus){ 
  params[[s]] <- base.draw(k = NX,J = 3,m0 = m0,s0 = s0,ms = ms,vs = vs,v = v,t = t)
}
#---------------------------------------------------

#algo 8
#initiate
nsim <- 10000
res <- list() #store params at each iteration
comp<-list() #store clus allocation at each iteration
count <- 0

#mcmc
#---------------------------------------------------

for( iter in 200:nsim){
  print(iter)
  
  #component allocations
  #---------------------------------------------------
  for(i in 1:N){
    #print(i)
    #compute likelihoods
    lik <- c()
    
    #existing clusters
    u.clus <- sort(unique(clus))
    for(u in u.clus){
      num <-length(which(clus[-i]==u))
      tmp <- num*f(x = xs[i,],y = cl[i],param = params[[u]])
      lik <- c(lik, tmp )  
    }
    
    #draw new parameters for m
    m.phi<-list()
    for(r in 1:m){
      m.phi[[r]] <- base.draw(k = NX,J = J,m0 = m0,s0 = s0,ms = ms,vs = vs,v = v,t = t)
      tmp <-alpha*f(x = xs[i,],y = cl[i],param = m.phi[[r]])/m
      lik <- c(lik,   tmp) 
    }
    
    #normalise
    p <-lik/sum(lik)
    
    #sample component allocation
    ind <- c(u.clus,(1:m))
    draw <- sample(1:length(ind),size=1,replace = F,prob =p)

    if(draw > length(u.clus)){
      clus[i]<- max(u.clus)+1
      params[[clus[i]]] <- m.phi[[ind[draw]]]
    }else{
      clus[i]<- ind[draw]
    }
    
  }
  #---------------------------------------------------
  
  #update parameters, MH
  #---------------------------------------------------
  for(q in 1:R){
    #existing clusters
    u.clus <- sort(unique(clus))
    for(u in u.clus){ 
      prop <- base.draw(k = NX,J = J,m0 = m0,s0 = s0,ms = ms,vs = vs,v = v,t = t)
      old.p <- params[[u]]
      
      #mh
      acc.rate1 <- 1
      acc.rate2 <- 1
      for(j in which(clus==u)){
        #update mu and sig
        num1 <- prod(sapply(1:NX, function(i) dnorm(xs[j,i], mean= prop$mu[i], sd=sqrt(prop$sig[i]))/prop$sig[i]))
        den1 <-  prod(sapply(1:NX, function(i) dnorm(xs[j,i], mean= old.p$mu[i], sd=sqrt(old.p$sig[i]))/old.p$sig[i]))
        temp1 <- num1/den1
        acc.rate1 <- acc.rate1 *temp1
        
        #update a nd b
        num2 <- py(x = xs[j,],y = cl[j],param = prop)
        den2 <- py(x = xs[j,],y = cl[j],param = old.p)
        temp2 <- num2/den2
        acc.rate2 <- acc.rate2 *temp2
      }
      
      if(runif(1)<acc.rate1){
        params[[u]]$mu <- prop$mu
        params[[u]]$sig <- prop$sig
        #print("MH1")
        #count <- count +1
      }
      if(runif(1)<acc.rate2){
        params[[u]]$a <- prop$a
        params[[u]]$b <- prop$b 
        #print("MH2")
      }
      
    }
  }
  #---------------------------------------------------
  
  #record
  #---------------------------------------------------
  if(iter>10){
    res[[iter-10]]<- params
    comp[[iter-10]]<- clus
  }
}
#---------------------------------------------------

#Analysis
#---------------------------------------------------

#test data
test.xs <- data[test.ind, ]
#actual class labels of test data
test.cl <- as.numeric(cl.all[test.ind])
#class results
test.r <- sapply(1:length(test.ind), function(t) classify(test.xs[t,], res,comp))
#accuracy
mean(test.r==test.cl)

#2 dimensional data charts
par(mfrow=c(1,2))
plot(xs, xlab="Dimension 1", ylab="Dimension 2")
for(i in 1:J) points (xs[which(cl==i),], col=i)
plot(test.xs,xlab="Dimension 1", ylab="Dimension 2")
dim(test.xs)
for(i in 1:J) points(test.xs[which(test.cl==i),], col=i)
for(i in 1:J) points(test.xs[which(test.r==test.cl & test.cl==i),],col=i, pch="x")
for(i in 1:J) points(test.xs[which(pred==test.cl & test.cl==i),],col=i, pch="x")
#---------------------------------------------------


