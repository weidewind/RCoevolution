graphs <- function (stats, simstats){
  
  lsim <- hist(simstats$maxdepth)
  lfake <- hist(stats$maxlength)
  #plot( lsim, col=rgb(0,0,1,1/4), xlim=c(0,350), freq = FALSE, ylim = c(0, 0.01), main = "maxdepths" )  # first histogram
  #plot( lfake, col=rgb(1,0,0,1/4), xlim=c(0,350), freq = FALSE,  ylim = c(0, 0.01), add=T)  # second
  old.par <- par(mfrow=c(2, 2))
  hist(stats[stats$mutations <=3, "pvalue_epistasis.median."], main = "median, epi, <=3")
  hist(stats[stats$mutations <=3, "pvalue_environment.median."], main = "median, env, <=3")
  hist(stats[stats$mutations >3, "pvalue_epistasis.median."], main = "median, epi, >3")
  hist(stats[stats$mutations >3, "pvalue_environment.median."], main = "median, env, >3")
  par <- old.par
  old.par <- par(mfrow=c(2, 2))
  hist(stats[stats$mutations <=3, "pvalue_epistasis.mean."], main = "mean, epi, <=3")
  hist(stats[stats$mutations <=3, "pvalue_environment.mean."], main = "mean, env, <=3")
  hist(stats[stats$mutations >3, "pvalue_epistasis.mean."], main = "mean, epi, >3")
  hist(stats[stats$mutations >3, "pvalue_environment.mean."], main = "mean, env, >3")
  par <- old.par
  #  plot(stats$pvalue_epistasis.median., stats$pvalue_epistasis.mean.)
  old.par <- par(mfrow=c(2, 2))
  hist(stats$pvalue_epistasis.median, main = "median, epi")
  hist(stats$pvalue_environment.median, main = "median, env")
  hist(stats$pvalue_epistasis.mean, main = "mean, epi")
  hist(stats$pvalue_environment.mean, main = "mean, env")
  par <- old.par
}

maxlen_graphs<- function (stats, simstats){
  hist(stats[stats$mutations >3 & stats$maxlength >150, "pvalue_environment.median."])
  hist(stats[stats$mutations >3 & stats$maxlength < 90, "pvalue_environment.median."])
  hist(stats[stats$mutations >3 & stats$maxlength >150, "pvalue_environment.mean."])
  hist(stats[stats$mutations >3 & stats$maxlength <90, "pvalue_environment.mean."])
}

print_adjusted <-function(expmetastats, nominals, output){
  epi_med_adjusted <- singles_fdr_plot( type = "median", value = "epistasis", nominal_pvalues = nominals$epi_median, expmetastats = expmetastats, onecol =  TRUE)
  epi_mean_adjusted <- singles_fdr_plot( type = "mean", value = "epistasis", nominal_pvalues = nominals$epi_mean, expmetastats = expmetastats, onecol =  TRUE)
  env_med_adjusted <- singles_fdr_plot( type = "median", value = "environment", nominal_pvalues = nominals$env_median, expmetastats = expmetastats, onecol =  TRUE)
  env_mean_adjusted <- singles_fdr_plot( type = "mean", value = "environment", nominal_pvalues = nominals$env_mean, expmetastats = expmetastats, onecol =  TRUE)
  matr <- cbind (epi_med,epi_mean,env_med,env_mean)
  write.csv(matr, output)
}

hists <- function (maxdepth = 50, group = "all",  value = "pvalue", breaks = "Sturges"){
  old.par <- par(mfrow=c(2, 2), mar = c(2, 4, 5, 2) + 0.1)
  hist(fake_fdr[fake_fdr$maxdepth == maxdepth & fake_fdr$group == group & fake_fdr$type == "mean_stat", paste(c("epi_", value), collapse = "")], main ="mean, epi",       breaks = breaks)
  hist(fake_fdr[fake_fdr$maxdepth == maxdepth & fake_fdr$group == group & fake_fdr$type == "mean_stat", paste(c("env_", value), collapse = "")], main ="mean, env",       breaks = breaks)
  hist(fake_fdr[fake_fdr$maxdepth == maxdepth & fake_fdr$group == group & fake_fdr$type == "median_stat", paste(c("epi_", value), collapse = "")], main ="median, epi",   breaks = breaks)
  hist(fake_fdr[fake_fdr$maxdepth == maxdepth & fake_fdr$group == group & fake_fdr$type == "median_stat", paste(c("env_", value), collapse = "")], main ="median, env",   breaks = breaks)
  title(paste (group, value, sep = ", "), outer=TRUE, line = -1)
  par(old.par)
}