###################
#Jake Rich
#Writing Project Code
#Spring 2018
###################
###Setting Up R Environment
#Installing necessary packages if not already installed.
list.of.packages <- c("lme4","MCMCglmm","coda","glmm","ggplot2","RCurl")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#Loading in dataset(s)
owl_dat <- read.table("D:/Google Drive/17 - STAT 532 Bayesian Statistics/Owls.txt",header = TRUE)
























































