#Authentication procedure has been omitted due to containing confidential secure keys and passwords

#retrieve and clean twitter data
#---------------------------------------------------
tweets <- homeTimeline( n=2000)

#RETRIEVE TWEETS
tweets<-homeTimeline( n=2000)

#CONVERT TO DF
tweets.df <- twListToDF(tweets)

#build corpus
myCorpus <- Corpus(VectorSource(tweets.df$text))

#convert to lower case
myCorpus <- tm_map(myCorpus, content_transformer(tolower))

# remove URLs
removeURL <- function(x) gsub("http[^[:space:]]*", "", x)

# tm v0.6
myCorpus <- tm_map(myCorpus, content_transformer(removeURL))

# remove anything other than English letters or space
removeNumPunct <- function(x) gsub("[^[:alpha:][:space:]]*", "", x)
myCorpus <- tm_map(myCorpus, content_transformer(removeNumPunct))

#stop words
myStopwords <- c(stopwords('english'),"rt", "indy","guardiansport","fastft")

#remove stopwords from corpus
myCorpus <- tm_map(myCorpus,  removeWords, myStopwords)

#remove extra space
myCorpus <- tm_map(myCorpus, stripWhitespace)

#copy corpus for dictionary
myCorpusCopy <- myCorpus
#---------------------------------------------------

# Term doc matrix, initial inspection, limit to interesting data
#---------------------------------------------------
tdm <- TermDocumentMatrix(myCorpus,control=list(wordLengths=c(1,Inf)))
findFreqTerms(tdm, lowfreq = 30)

tms<-as.matrix(tdm)

#index of words of interest
pop.terms <-  which(rownames(tms) %in% 
                      c("trump","corbyn","china", "sex"
                        ,"europe","migrant","putin","cameron"
                        ,"sanders","eu"
                      )
)
#---------------------------------------------------

#Output data matrix of tokenised words by tweet doc. containing only useful words
#---------------------------------------------------

#which documents contain these words
twts <- c() #index of tweets containing pop.terms
for(w in pop.terms) twts <- c(twts, which(tms[w,]>0))
twts <- unique(twts)
#matrix restricted to tweets containing popular terms
tms <- tms[,twts]
J <- dim(tms)[2]

#data list. A vector for each tweet containing words used
data <- list()
for(j in 1:J){
  ws <- which(tms[,j]>0)
  d<- c()
  for(w in ws){
    d <- c(d, rep(w,tms[w,j]))
  }
  data[[j]]<-d
}

#number of unique words, vocab
V <- length(which(rowSums(tms)>0))
#---------------------------------------------------

