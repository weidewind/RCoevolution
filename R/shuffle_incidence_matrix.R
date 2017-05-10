source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("BiRewire"))

data <- read.table("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_incidence_matrix", sep="", colClasses='character')
dat <-sapply(data, function(e) {
  splitter <- strsplit(e,"")
})
dat <- sapply(dat, function (e) {
  as.numeric(unlist(e))
})
dat <- t(dat)



for (i in 1:10000){
  sink(file=paste("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/MockShuffledMatrices/", i, sep=""))
  m2<-birewire.rewire.bipartite(dat,verbose=FALSE)
  subs_on_node <-apply(m2, 1, function(e){
    which(e>0)
  })
  sapply(subs_on_node, function(e){
    cat (paste(e, collapse=","), "\n") 
  })
  sink()
}