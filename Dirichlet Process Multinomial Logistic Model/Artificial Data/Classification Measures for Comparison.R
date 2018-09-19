#classification measures
#multi class classification measures

#true positives
tp <- function(pred,act, cl){
  ind <- which(act==cl & pred==cl)
  length(ind)
} 

#true negatives
tn <- function(pred, act, cl){
  ind <- which(act!=cl & pred != cl)
  length(ind)
}

#false negatives
fn <- function(pred, act, cl){
  ind <- which(act==cl & pred != cl)
  length(ind)
}

#false positives
fp <- function(pred,act,cl){
  ind <- which(pred==cl & act != cl)
  length(ind)
}

#recall
recall<- function(pred,act,cl){
  num <- tp(pred,act,cl)
  den <- (tp(pred,act,cl)+fn(pred,act,cl))
  if(den == 0) den <-1
  return(num/den)
}

#precision
precision <- function(pred,act,cl){
  num <- tp(pred,act,cl)
  den <- (tp(pred,act,cl)+fp(pred,act,cl))
  if(den ==0) den <- 1
  return(num/den)
}

#micro precision
micro.prec <- function(pred, act,cls){
  cls <- sort(unique(cls)) 
  num <- sum(sapply(cls, function(x) tp(pred,act,x)))
  den <- sum(sapply(cls, function(x) fn(pred,act,x)+tp(pred,act,x)))
  if(den==0) den <- 1
  #print(den)
  return(num/den)
}

#micro recall
micro.rec <- function(pred, act,cls){
  cls <- sort(unique(cls))
  
  num <- sum(sapply(cls, function(x) tp(pred,act,x)))
  
  den <- sum(sapply(cls, function(x) fp(pred,act,x)+tp(pred,act,x)))
  #print(den)
  if(den==0) den <- 1
  return(num/den)
}

#specificity
spec <- function(pred, act,cl){
  num <- tn(pred,act,cl)
  den <- tn(pred,act,cl) + fp(pred,act,cl)
  
  return(num/den)
}