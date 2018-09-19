

#number of data pointa
N<-3*100
######################
#data
tm1 <-  c(1,-1)
tm2<- c(0,0)
tm3 <-  c(2,2)

tm <- rbind(tm1,tm2,tm3)
sig <- matrix(c(1,0.5,0.5,1),2,2)
data <- matrix(NA, N,2) #data
#cl <-numeric(N)
for(i in 1:3){
  data[(1+(i-1)*N/3):(N/3*i),1:2] <- rmvnorm(N/3,mean = tm[i,],sigma =sig)
  #cl[(1+(i-1)*N/3):(N/3*i),3]<-i
}
#true classes
cl.all<-c(rep(1,100),rep(2,100),rep(3,100))
#set clusters to be true classes
cl <-cl.all
#set covariates as xs
xs<- data
#plot test data
par(mfrow=c(1,2))
plot(xs, xlab="Dimension 1", ylab="Dimension 2")
for(i in 1:3) points(xs[which(cl==i),],col=i)


