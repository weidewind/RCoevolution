library(ggplot2)
draw_graph <- function(prot){

    realdata <- read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/",prot,"_r_150_halfdepth.csv"), collapse=""))
    data <- read.csv(paste(c("C:/Users/weidewind/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/",prot,"_noneed_simulation_150.csv"), collapse=""))
#_noneed_simulation_150.csv
    means = lapply(data[grep('obs',names(data))], function(e) {mean(e,na.rm = TRUE)})
    means <- means[grep('X$',names(means), invert=TRUE)]
  # means <- means[grep('32',names(means), invert=TRUE)]
    means <- unlist(means)
    names(means)<-NULL

    quantiles <- lapply(data, function(e) {
                             quantile(e, 
                                      probs = c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95,0.99),
                                      na.rm = TRUE)
    })
    obs_quantiles <- quantiles[grep('obs',names(quantiles))]
 #   obs_quantiles <- obs_quantiles[grep('32',names(obs_quantiles), invert=TRUE)]
    obs_quantiles <- obs_quantiles[grep('X$',names(obs_quantiles), invert=TRUE)]
#plot(seq(10,310, by=10),lapply(obs_quantiles, function(e) {e["5%"]}), type="l")

    sim_5 <- lapply(obs_quantiles, function(e) {e["5%"]})
    sim_95 <- lapply(obs_quantiles, function(e) {e["95%"]})
    y=c(c(realdata$obs,rep(0,32-length(realdata$obs))), c(realdata$exp,rep(0,32-length(realdata$exp))), c(sim_5,rep(0,32-length(sim_5))),  c(means,rep(0,32-length(means))),  c(sim_95,rep(0,32-length(sim_95))))
    y <- unlist(y)
    names(y)<-NULL

    plotdata <- data.frame(distance=rep(seq(10,320, by=10),5),
                          y=y,
                          group=c(rep("real_obs", 32),rep("real_exp", 32),rep("5%", 32),rep("mean", 32),rep("95%", 32)))
    ggplot(data=plotdata, aes(x=as.integer(distance), y=y, group=group, color=group)) + geom_line(aes(size = group))  + scale_size_manual(values = c( 0.1,0.1,2,0.1, 2))  + geom_point() + ggtitle(paste(c(prot, "_old_simulation_150"),collapse=""))

}
