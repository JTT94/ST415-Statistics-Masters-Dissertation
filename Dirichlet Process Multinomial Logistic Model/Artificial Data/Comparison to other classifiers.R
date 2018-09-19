#naive bayes
#---------------------------------------------------
m.nb <- naiveBayes(class~p1+p2
                   , data=data.frame(p1=data[train.ind,1], p2=data[train.ind,2], class =factor(cl.all[train.ind])),
                   laplace =1
)
#---------------------------------------------------

#random forest
#---------------------------------------------------
m.rf <- randomForest(class~p1+p2
                     , data=data.frame(p1=data[train.ind,1], p2=data[train.ind,2], 
                                       class =factor(cl.all[train.ind]))
)
#---------------------------------------------------

#svm, RBF kernel
#---------------------------------------------------
m.svm <- svm(class~p1+p2
             , data=data.frame(p1=data[train.ind,1], p2=data[train.ind,2], class =factor(cl.all[train.ind])),
             type= "C"
)
#---------------------------------------------------

#svm, linear kernel
#---------------------------------------------------
m.svml <- svm(class~p1+p2
             , data=data.frame(p1=data[train.ind,1], p2=data[train.ind,2], class =factor(cl.all[train.ind])),kernel="linear"
             ,type= "C"
)
#---------------------------------------------------

#Some predction, accurary
#---------------------------------------------------
pred <- predict(m.nb, 
                newdata=data.frame(p1=data[test.ind,1], p2=data[test.ind,2]),
                type="class"
                )
mean(pred==cl.all[test.ind])

pred <- predict(m.svm,
                newdata=data.frame(p1=data[test.ind,1], p2=data[test.ind,2])
)
pred <- predict(m.svml,
                newdata=data.frame(p1=data[test.ind,1], p2=data[test.ind,2])
)
mean(pred==cl.all[test.ind])

pred <- predict(m.rf,
                newdata=data.frame(p1=data[test.ind,1], p2=data[test.ind,2]), type="class"
)
mean(pred==cl.all[test.ind])




