#useful functions
#---------------------------------------------------

#likeli
#number of words w across all documents in component k
r.w <-function(w,k) sum(sapply(1:J, function(tj) length(which(tk[[tj]][tab[[tj]][which(data[[tj]]==w)]]==k))))

#likelihood for a single data-pint
f <- function(k,i,j,alpha,V){
  w <- data[[j]][i]
  n.w <- r.w(w,k)
  num <- alpha/V+n.w
  den <- alpha+n.w
  return(num/den)
}

#likelihood a table
f.t <- function(k,t,j, alpha,V){
  
  w.t <- data[[j]][(which(tab[[j]]==t))]
  u.w <- sort(unique(w.t))
  ns <- sapply(1:length(u.w), function(x) length(which(w.t==u.w[x])))
  n <- length(w.t)
  rs <- sapply(1:length(u.w), function(x) r.w(u.w[x],k))
  R <- sum(rs)
  
  if(tk[[j]][t]==k){
    num <- (sapply(1:length(u.w), function(x) gamma(rs[x]+alpha/V)))
    den <-  (sapply(1:length(u.w), function(x) gamma(rs[x]+alpha/V-ns[x])))
    return(prod(num/den)/prod(alpha+R-n-1+1:n))
  }else{
    num <- (sapply(1:length(u.w), function(x) gamma(rs[x]+alpha/V+ns[x])))
    den <-  (sapply(1:length(u.w), function(x) gamma(rs[x]+alpha/V)))
    return(prod(num/den)/prod(alpha+R-1+1:n))
  }
}

#likelihood for new table
f.t.new <- function(t,j, alpha,V){
  w.t <- data[[j]][(which(tab[[j]]==t))]
  u.w <- sort(unique(w.t))
  ns <- sapply(1:length(u.w), function(x) length(which(w.t==u.w[x])))
  n <- length(w.t)
  
  p21 <- (sapply(1:length(u.w), function(x) gamma(ns[x]+alpha/V)))
  p22 <- (sapply(1:length(u.w), function(x) gamma(alpha/V)))
  
  p1 <- gamma(alpha)/gamma(alpha+n)
  p2 <-prod(p21/p22)
  
  return(p1*p2)
}
#---------------------------------------------------

# Initialise, notation heavy
#---------------------------------------------------
# list of vectors, one for each document
# the length of tab[[j]] is the number of words in document j, 
# the ith entry of tab[[j]], is the table allocation of the ith word in doc j
tab <- list()
# list of vectors one for each document, length is the number of tables in each corresponding document
# the ith table in document j has component allocation tk[[j]][i]
tk <- list()

#load twitter data

#initiate tables and components
for(j in 1:J){
  tab[[j]] <- sample(1:5,size=length(data[[j]]),replace=T)
  tabs <- sort(unique(tab[[j]])) 
  tk[[j]] <- tabs
  for(u in 1:max(tabs)){
    tk[[j]][u]<- j#sample(1:10,size=1)
  }
} 

#initiate hyperparameters and number of simulations
alpha <- 0.5 # concentraion paramater for Dp G_j
nsim <- 100
g<-0.5 #concentraion param for G_0

#lists to store results
res.tk <- list()
res.tab <- list()
#---------------------------------------------------

#run MCMC simulations
for(iter in 1:nsim){
  print(iter)
  #for document J
  #---------------------------------------------------
  for(j in 1:J){
   # print("Document")
    #print(j)
    
    #for word i in document j
    #---------------------------------------------------
    for(i in 1:length(data[[j]])){
      #component allocations
      ks <- unlist(sapply(1:J,function(j) sapply(tab[[j]], function(x) tk[[j]][x])))
      u.ks <- sort(unique(ks))
      ms <- sapply(1:length(u.ks), function(x) length(which(ks==u.ks[x])))
      
      #draw new table allocation
      u.t <- sort(unique(tab[[j]]))
      num <- sapply(1:length(u.t), function(x) length(which(tab[[j]][-i]==u.t[x])))
      lik <- sapply(1:length(u.t), function(x) f(k = tk[[j]][x],i =i,j = j,alpha = alpha,V = V ))
      
      choices <- 1:length(u.t)
      choices <- c(choices,length(u.t)+1)
      
      lik.k <- sapply(1:length(u.ks),function(x) f(k = x,i = i,j = j,alpha = alpha,V = V) )
      p.new <- (g/V + sum(ms*lik.k))/(sum(ms)+g)
      prob <- c(num*lik,  alpha*p.new)
      prob <- prob/sum(prob)
      
      draw <- sample(choices, prob=prob, size=1)
      if(draw > length(u.t)){
        tab[[j]][i]<- max(u.t)+1
        #new component allocation
        prob <- c((ms*lik.k),g/V)
        prob <- prob/sum(prob)  
        choices <- c(u.ks, max(u.ks)+1)
        draw <- sample(choices, size=1, prob =prob)
        tk[[j]] <- c(tk[[j]], draw)
      }else{
        tab[[j]][i]<- u.t[draw]
      }
    }
    #---------------------------------------------------
    
   
    #component allocations for each table
    #---------------------------------------------------
    #count component allocations
    ks <- unlist(sapply(1:J,function(j) sapply(tab[[j]], function(x) tk[[j]][x])))
    u.ks <- sort(unique(ks))
    ms <- sapply(1:length(u.ks), function(x) length(which(ks==u.ks[x])))
    u.t <- sort(unique(tab[[j]]))
    for(t in 1:length(u.t)){
      t.k <- tk[[j]][t]
      n.t <- length(which(tab[[j]]==t))
      #component allocations
      ks <- unlist(sapply(1:J,function(x) sapply(tab[[x]], function(x2) tk[[x]][x2])))
      u.ks <- sort(unique(ks))
      ms <- sapply(1:length(u.ks), function(x) length(which(ks==u.ks[x])))
      ms[which(u.ks==t.k)] <- ms[which(u.ks==t.k)]-n.t
      
      lik.k <- sapply(1:length(u.ks), function(x) f.t(k = u.ks[x],t = u.t[t],j = j,alpha = alpha, V = V))
      p.new <- f.t.new(t = u.t[t],j = j,alpha = alpha,V = V)*g
      prob <- c(ms*lik.k, p.new)
      prob <- prob/sum(prob)
      choices <- c(u.ks, max(u.ks)+1)
      tk[[j]][u.t[t]] <- sample(choices, size=1, prob=prob)
    }
    #---------------------------------------------------
    
  }
  
  #store results
  res.tk[[iter]]<- tk
  res.tab[[iter]] <- tab
  #print table of components
  allk <-c()
  for(i in 1:J) allk<-c(allk,(res.tk[[iter]][[i]][unique(res.tab[[iter]][[i]])]))
  print(table(allk))
}
#---------------------------------------------------

#Analysis
#---------------------------------------------------

#post mcmc functions
#words belonging to component k in iteration it
wk <- function(k,it){
  res <- c()
  for(j in 1:J){
    tmp <- rownames(tms)[data[[j]][which(res.tk[[it]][[j]][res.tab[[it]][[j]]]==k)]]
    if(length(tmp)>0)
      res <- c(res,tmp)
    
  } 
  return(res)
} 
#number of documents with words in component k
doc.count <- function(k,it){
  res <- 0
  for(j in 1:J){
    tmp <- rownames(tms)[data[[j]][which(res.tk[[it]][[j]][res.tab[[it]][[j]]]==k)]]
    if(length(tmp)>0)
      res <- res+1
    
  } 
  return(res)
} 
#number of words associated to component k 
word.count <- function(k,it){
  res <- 0
  for(j in 1:J){
    tmp <- rownames(tms)[data[[j]][which(res.tk[[it]][[j]][res.tab[[it]][[j]]]==k)]]
    if(length(tmp)>0)
      res <- res+length(tmp)
    
  } 
  return(res)
} 

#table of for counts of words in each component underiteration it
allk <-c()
for(i in 1:J) allk<-c(allk,(res.tk[[it]][[i]][unique(res.tab[[it]][[i]])]))
(table(allk))


#post mcmc analysis
for(i in comps) print(doc.count(i,it))
for(i in comps) print(word.count(i,it))
par(mfrow=c(1,1))
for(i in comps) wordcloud(wk(i,it),random.order=F,min.freq=0)
sum(sapply(comps, function(x) doc.count(x,it)))
#---------------------------------------------------