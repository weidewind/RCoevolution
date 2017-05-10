
cs <- read.csv("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/environment_or_epistasis/sanityCheck/h1_all_fdr")
r100 <- cs[cs$restr == 100, ]
r100median <- r100[r100$stat == "median_stat", ]
r100mean <- r100[r100$stat == "mean_stat", ]
hist(r100median$epi)
hist(r100median$env)
hist(r100mean$epi)
hist(r100mean$env)



