library(KernSmooth)
dist_file <- "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_nsyn_nobins"
ints_file <- "C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_ints_nobins"
dist_vectors <- read.csv (dist_file,  header = F, blank.lines.skip=T, stringsAsFactor=F)
ints_vectors <- read.csv (ints_file,  header = T, blank.lines.skip=T, stringsAsFactor=F)
#cleaned_dist<-dist_vectors[dist_vectors$V3>0,]
cleaned_dist<-dist_vectors[dist_vectors$V1>0,]
cleaned_ints<-ints_vectors[ints_vectors$X95_percentile < 19499 & ints_vectors$mean != 0,]
#h<-dpill(cleaned_dist$V1, cleaned_dist$V3-cleaned_dist$V2) # is a bandwidth
lp<-locpoly(cleaned_dist$V1, cleaned_dist$V3-cleaned_dist$V2, bandwidth=10, range.x=c(1,430))

cleaned_ints$distance <- as.integer(cleaned_ints$distance)
cleaned_dist$V1 <- as.integer(cleaned_dist$V1)

data <-merge(cleaned_ints, cleaned_dist, by.x="distance", by.y="V1", all.y=T)
data<-data.frame(data, data$V3-data$V2)

neu<-data.frame(distance=c(rep(data$distance, 3), lp$x),
                y=c(data$data, data$X05_percentile, data$X95_percentile, lp$y),
                group=c(rep("V", nrow(data)),rep("L", nrow(data)), rep("U", nrow(data)), rep("S", length(lp[[1]]))))
#L - lower, U- upper, V- observed value, S - smoothed
ggplot(data=neu, aes(x=as.integer(distance), y=y, group=group, color=group))+geom_line()

#lp<-locpoly(cleaned_ints$distance, cleaned_ints$X95, bandwidth=25, range.x=c(1,430))
#smoothed upper limit of the confidence interval

upper <- predict(fit.loess, data.frame(x=new.x))+
  predict(fit.loess, data.frame(x=new.x), se=TRUE)$se.fit*qnorm(1-.05/2)

#with smoothed bootstrap CI
lp<-locpoly(cleaned_dist$V1, cleaned_dist$V3-cleaned_dist$V2, bandwidth=10, range.x=c(1,430))
lpl<-locpoly(cleaned_ints$distance, cleaned_ints$X05, bandwidth=10, range.x=c(1,430))
lpu<-locpoly(cleaned_ints$distance, cleaned_ints$X95, bandwidth=10, range.x=c(1,430))
neuer<-data.frame(distance=c(rep(data$distance, 3), lp$x, lpu$x, lpl$x),
                  y=c(data$data, data$X05_percentile, data$X95_percentile, lp$y, lpu$y, lpl$y),
                  group=c(rep("V", nrow(data)),rep("L", nrow(data)), rep("U", nrow(data)), rep("S", length(lp[[1]])), rep("SU", length(lpu[[1]])), rep("SL", length(lpl[[1]]))))
ggplot(data=neuer, aes(x=as.integer(distance), y=y, group=group, color=group))+geom_line()


# loess and its confidence intervals!!
diff<-cleaned_dist$V3-cleaned_dist$V2
fit.loess <-loess( diff~cleaned_dist$V1, span = 0.75, degree=1,family="gaussian")
lupper <- predict (fit.loess, data.frame(x=c(1:450))) + predict (fit.loess, data.frame(x=c(1:450)),se=TRUE)$se.fit*qnorm(1-0.05/2)
llower <- predict (fit.loess, data.frame(x=c(1:450))) - predict (fit.loess, data.frame(x=c(1:450)),se=TRUE)$se.fit*qnorm(1-0.05/2)
neuer<-data.frame(distance=c(c(1:450), c(1:450), c(1:450)),
y=c(predict (fit.loess, data.frame(x=c(1:450))), llower, lupper),
group=c(rep("loess", 450),rep("lower", 450),rep("upper", 450)))
ggplot(data=neuer, aes(x=as.integer(distance), y=y, group=group, color=group))+geom_line()