con<-file('C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_shuffler_10000_1')
open(con)
t<-read.table(con, skip=1, nrow=1)
v<-vapply(c(1:450), function(elm) {read.table(con, skip=10000, nrow=1)[[1]]}, numeric(1))
intlp<-locpoly(c(0:449), v, bandwidth=25, range.x=c(1,449))

data<-data.frame("distance"=intlp$x)

for (i in 2:100){
seek(con, where=0)
t<-read.table(con, skip=i, nrow=1)
v<-vapply(c(1:450), function(elm) {read.table(con, skip=10000, nrow=1)[[1]]}, numeric(1))
intlp<-locpoly(c(0:449), v, bandwidth=25, range.x=c(1,449))
data<-data.frame(data, intlp$y)

}
close(con)


#Kernel smooth all bootstrap @observations@, and olny then compute confidence intervals for the smoothed values.
file <- read.csv("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_shuffler_10000_1_trans")
lp<-locpoly(c(2:401), file[[2]], bandwidth=25, range.x=c(2,401))
lps<-lapply(file[1:500], function(elm) {lp<-locpoly(c(2:401), elm, bandwidth=25, range.x=c(2,401)) 
                            lp$y})
#lps[[i]][1]# (i>1) is a vector representing all smoothed points at distance 1 (2.14 in this case)
intervals<-c(1:400) #was 1:400
quantiles <- sapply(intervals, function(i){ 
dist<-sapply(lps, function(elm) {elm[i]})
quantile(dist, probs = c(0.05, 0.95))
})
plot(lp$x[1:400], quantiles[1,]) #was [1:400]

plotdata <- data.frame(distance=c(c(1:450), rep(lp$x[1:400],2)),
                       y=c(predict (fit.loess, data.frame(x=c(1:450))), quantiles[1,], quantiles[2,]),
                       group=c(rep("loess", 450),rep("bootstrap_lower", 400),rep("bootstrap_upper", 400)))
ggplot(data=plotdata, aes(x=as.integer(distance), y=y, group=group, color=group))+geom_line()


## Kernel smooth all bootstrap @observations@, and olny then compute confidence intervals for the smoothed values.
## First two observations (radius = 0, 1) are trimmed.
file <- read.csv("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_shuffler_10000_1_trans")
lp<-locpoly(c(2:400), file[[2]][3:401], bandwidth=25, range.x=c(2,400)) # get new set of distances (points after smoothing)
# 2- first distance with some observations. 3 - number of the corresponding line in the file (skip 1 and 2)
lps<-lapply(file[1:500], function(elm) {lp<-locpoly(c(2:400), elm, bandwidth=25, range.x=c(2,400)) # take only 500 observations from 10000
                                        lp$y})
#lps[[i]][1]# (i>1) is a vector representing all smoothed points at distance 1 (2.14 in this case)
intervals<-c(1:(length(lp$x)+1))
quantiles <- sapply(intervals, function(i){ 
  dist<-sapply(lps, function(elm) {elm[i]})
  quantile(dist, probs = c(0.05, 0.95))
})
plot(lp$x, quantiles[1,])