#load wine data
#dataset attributes
data <- w.data

#split data into test and train
#random permutation of indicies
n <- dim(data)[1]
r.ind <- sample(1:n,size=n, replace=F)
#split 1:n into k equal groups
k.cuts <- cut(r.ind, breaks = 3, labels=F)
#use position of each element of each group in random permutation as new indicies
k.ind <- lapply(1:k, function(x) which(k.cuts==x))
test.ind <- which(k.cuts!=1)
train.ind <- which(k.cuts ==1)

#covariates
xs <- data[train.ind,]
#classes
cl <- w.class[train.ind]


#initialise
alpha<-1
NX <-k <- dim(xs)[2]
m0 <- 0
s0 <- 1
ms <- 0
vs <- 1
J <- length(unique(cl))
m<- 3
N <- dim(xs)[1]
R <- 5
count <- 0
v <- 1
t <- 1
clus <- cl
u.clus <- sort(unique(clus)) #unique clusters


params <-  list() #unique params 
for(s in u.clus){ 
  params[[s]] <- base.draw(k = NX,J = 3,m0 = m0,s0 = s0,ms = ms,vs = vs,v = v,t = t)
}


#algo 8
######################
#initiate
nsim <- 10000
res <- list() #store params at each iteration
comp<-list() #store clus allocation at each iteration
####################

# #mcmc
#mcmc
for( iter in 1:nsim){
  print(t)
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
  
  #update parameters
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
        num1 <- prod(sapply(1:NX, function(x) dnorm(xs[j,x], mean= prop$mu[x], sd=sqrt(prop$sig[x]))/prop$sig[x]))
        den1 <-  prod(sapply(1:NX, function(x) dnorm(xs[j,x], mean= old.p$mu[x], sd=sqrt(old.p$sig[x]))/old.p$sig[x]))
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
  
  
  #record
  if(t>1){
    res[[t-1]]<- params
    comp[[t-1]]<- clus
  }
}

#analyse results
test.ind <- which(k.cuts!=1)#actual labels
test.act <- (w.class[test.ind])
test.xs <- w.data[test.ind,]#test covariates
#model class predictions
test.r <- sapply(1:dim(test.xs)[1], function(t) classify(test.xs[t,], res,comp))
#accuracy
mean(test.act==test.r)

##Comparison to other predictors
#naive bayes
m.nb <- naiveBayes(as.factor(w.class[train.ind])~ .
                     , data=w.data[train.ind,]
)


#random forest
m.rf <- randomForest(as.factor(w.class[train.ind])~ .
                     , data=w.data[train.ind,]
)


#svm
m.svml <- svm(w.class[train.ind]~ .
             , data=data[train.ind,], kernel = "linear",
             type ="C", scale =T, gamma =1, cost=1
)
#svm
m.svm <- svm(w.class[train.ind]~ .
             , data=data[train.ind,], kernel = "radial",
             type ="C", scale =T, gamma =1, cost=1
)
pred <- predict(m.svml, newdata=data[test.ind,])
mean(pred==test.act)

pred <- predict(m.svm, newdata=data[test.ind,])
mean(pred==test.act)

pred <- predict(m.nb, newdata=w.data[test.ind,])
mean(pred==test.act)

pred <- predict(m.rf, newdata=w.data[test.ind,])
mean(pred==test.act)

mean(test.r==test.act)

#repeat the following for each "pred" for each model
mean(sapply((1:3) , function(x) precision(pred = pred,act = test.act,cl = x)))
mean(sapply((1:3) , function(x) recall(pred = pred,act = test.act,cl = x)))
