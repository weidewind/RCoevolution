install.packages("Hmisc")
library(foreign)
library(ggplot2)
library(MASS)
library(Hmisc)
library(reshape2)

    prot = "h1"
    ## beautifully explained at http://www.ats.ucla.edu/stat/r/dae/ologit.htm

    mytable <- read.csv(paste(c("C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/",prot,"_trend_table.csv"), collapse=""))
    
    ## fit ordered logit model and store results 'model'
    model <- polr(trend ~ surface + koel + epitope + trailing + leading + number, data = mytable, Hess=TRUE)
    
    ## view a summary of the model
    summary(model)
    
    ## Some people are not satisfied without a p-value. 
    ## One way to calculate a p-value in this case is by comparing the t-value against the standard normal distribution, like a z test. 
    ##  Of course this is only true with infinite degrees of freedom, but is reasonably approximated by large samples,
    ## becoming increasingly biased as sample size decreases. 
    ## This approach is used in other software packages such as Stata and is trivial to do.
    ## First we store the coefficient table, then calculate the pvalues and combine back with the table.
    
    ## store table
    (ctable <- coef(summary(model)))
    
    ## calculate and store p values
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    
    ## combined table
    (ctable <- cbind(ctable, "p value" = p))
    
    ## We can also get confidence intervals for the parameter estimates. 
    ## These can be obtained either by profiling the likelihood function 
    ## or by using the standard errors and assuming a normal distribution.
    ## Note that profiled CIs are not symmetric (although they are usually close to symmetric).
    ## If the 95% CI does not cross 0, the parameter estimate is statistically significant.
    
    (ci <- confint(model)) # default method gives profiled CIs
    confint.default(model) # CIs assuming normality



