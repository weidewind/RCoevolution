library("fdrtool")
#h1pvals <- read.csv("../input/h1pval.csv")
prot <- "n2"
h1pvals <- read.table(paste("../input/", prot, "_sites", sep = ""), sep="\t")
#h1pvals <-as.numeric(h1pvals)


fdrs <- lapply(levels(factor(h1pvals$V1)), function(restriction){
  dfr <- h1pvals[h1pvals$V1 == restriction,]
  fdrV5 <- fdrtool(dfr$V5,statistic="pvalue")
  fdrV6 <- fdrtool(dfr$V6,statistic="pvalue")
  fdrV7 <- fdrtool(dfr$V7,statistic="pvalue")
  fdrV8 <- fdrtool(dfr$V8,statistic="pvalue")
  newdfr <- data.frame(dfr$V1,dfr$V2, dfr$V3, dfr$V4, dfr$V5, dfr$V6, dfr$V7, dfr$V8, fdrV5$qval, fdrV6$qval, fdrV7$qval,fdrV8$qval)
  names(newdfr) <- c("depth_restriction",  "ancestor_site_node",	"number_of_mutations_in_the_subtree",	"maxdepth",	"epistasis_(median)_pvalue",  "epistasis_(mean)_pvalue",  "env_(median)_pvalue", "env_(mean)_pvalue", "epimed_FDR","epimean_FDR", "envmed_FDR","envmean_FDR")
  write.table(newdfr, paste( prot, "fdrs", restriction, sep="_"), quote=FALSE, row.names = F)
})


#cat(fdrs$qval, sep="\n")
