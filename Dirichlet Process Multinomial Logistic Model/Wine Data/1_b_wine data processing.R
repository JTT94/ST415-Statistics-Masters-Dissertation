#load wine data- pre processing
wine <- read.csv("C:/Users/James/Dropbox/Dissertation/R-code/DPMNL/Wine/wine.data", header=FALSE)
w.class <- wine[,1]
w.data <- wine[,-1]

#normalise data
for(att in 1:dim(w.data)[2]){
  w.data[,att] <- (w.data[,att]-mean(w.data[,att]))/sqrt(var(w.data[,att]))
}
w.data<- w.data[,c(1,10,7,5)]

#summary of wine data
names(w.data) <- c("Alcohol", "Colour Intensity", "Flavanoids", "Magnesium")
plot(w.data)
par(mfrow=c(2,3))
plot(w.data[,c(1,2)], pch="")
for(i in 1: length(sort(unique(w.class)))){
  points(w.data[which(w.class==i),c(1,2)], col=i)
}
plot(w.data[,c(2,3)], pch="")
for(i in 1: length(sort(unique(w.class)))){
  points(w.data[which(w.class==i),c(2,3)], col=i)
}
plot(w.data[,c(3,4)], pch="")
for(i in 1: length(sort(unique(w.class)))){
  points(w.data[which(w.class==i),c(3,4)], col=i)
}
plot(w.data[,c(1,3)], pch="")
for(i in 1: length(sort(unique(w.class)))){
  points(w.data[which(w.class==i),c(1,3)], col=i)
}
plot(w.data[,c(1,4)], pch="")
for(i in 1: length(sort(unique(w.class)))){
  points(w.data[which(w.class==i),c(1,4)], col=i)
}
plot(w.data[,c(2,4)], pch="")
for(i in 1: length(sort(unique(w.class)))){
  points(w.data[which(w.class==i),c(2,4)], col=i)
}
