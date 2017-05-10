library(KernSmooth)
library(ggplot2)
dist_file <- "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h3_nsyn_nobins"
ints_file <- "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/n2_ints_nobins"

file <- read.csv("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/n2_shuffler_10000_1_trans", header=F)

dist_vectors <- read.csv (dist_file,  header = F, blank.lines.skip=T, stringsAsFactor=F)
ints_vectors <- read.csv (ints_file,  header = T, blank.lines.skip=T, stringsAsFactor=F)

cleaned_dist<-dist_vectors[dist_vectors$V1>0,]
cleaned_ints<-ints_vectors[ints_vectors$X95_percentile < 19499 & ints_vectors$mean != 0,]
cleaned_ints$distance <- as.integer(cleaned_ints$distance)
cleaned_dist$V1 <- as.integer(cleaned_dist$V1)


diff<-cleaned_dist$V3-cleaned_dist$V2
fit.loess <-loess( diff~cleaned_dist$V1, span = 0.5, degree=1)
lupper <- predict (fit.loess, data.frame(x=c(1:450))) + predict (fit.loess, data.frame(x=c(1:450)),se=TRUE)$se.fit*qnorm(1-0.05/2)
llower <- predict (fit.loess, data.frame(x=c(1:450))) - predict (fit.loess, data.frame(x=c(1:450)),se=TRUE)$se.fit*qnorm(1-0.05/2)


#Kernel smooth all bootstrap @observations@, and olny then compute confidence intervals for the smoothed values.

#kernSmooth smoothing
#lp<-locpoly(c(2:401), file[[2]], bandwidth=25, range.x=c(2,401))
#lps<-lapply(file[1:500], function(elm) {lp<-locpoly(c(2:401), elm, bandwidth=25, range.x=c(2,401)) 
#                                        lp$y})
lp<-c(1:400)
lp$x<-c(1:400)
lps<-lapply(file[1:500], function(elm) {predict(loess( elm~c(2:460), span = 0.75, degree=1,family="gaussian"), 
                                                data.frame(x=c(2:460)))})

#this is a vector representing all smoothed points at distance 1 (2.14 in this case): lps[[i]][1]# (i>1) 
intervals<-c(1:400) #was 1:400
quantiles <- sapply(intervals, function(i){ 
  dist<-sapply(lps, function(elm) {elm[i]})
  quantile(dist, probs = c(0.05, 0.95))
})


plotdata <- data.frame(distance=c(rep(c(1:450),3), rep(lp$x[1:400],2)),
                       y=c(predict (fit.loess, data.frame(x=c(1:450))), llower, lupper, quantiles[1,], quantiles[2,]),
                       group=c(rep("loess", 450),rep("lower", 450),rep("upper", 450),rep("bootstrap_lower", 400),rep("bootstrap_upper", 400)))
ggplot(data=plotdata, aes(x=as.integer(distance), y=y, group=group, color=group))+geom_line()
