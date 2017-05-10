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
                              #print (paste(e, sep=","),  quote = FALSE)
                              #write (paste(e, sep=","), "")
                               cat (paste(e, collapse=","), "\n") 
                            })
    sink()
}


## shuffle without constraints
for (i in 1:100){
  sink(file=paste("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/CompletelyShuffledMatrices/", i, sep=""))
  m2<-dat[sample.int(nrow(dat)),]
  m2<-m2[,sample.int(ncol(m2))]
  subs_on_node <-apply(m2, 1, function(e){
    which(e>0)
  })
  sapply(subs_on_node, function(e){
    #print (paste(e, sep=","),  quote = FALSE)
    #write (paste(e, sep=","), "")
    cat (paste(e, collapse=","), "\n") 
  })
  sink()
}

