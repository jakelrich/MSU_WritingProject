###################
#Jake Rich
#Writing Project Code
#Spring 2018
###################
###Setting Up R Environment
#Installing necessary packages if not already installed.
######
#NOTE: You must have JAGS installed before using any of the Bayesian analyses
#in this script.
######
list_of_packages <- c("lme4","MCMCglmm","coda","glmm","ggplot2","RCurl", "simr","devtools","rjags","mcemGLM","MASS","glmmML")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, library, character.only = TRUE)
#Installing GLMMmisc from Github. GLMMmisc is a package of GLMM diagnostic tools
#and other functions that are useful with GLMMs. 
if (!require("GLMMmisc")) devtools::install_github("pcdjohnson/GLMMmisc")
library(GLMMmisc)
#Loading in dataset(s)
owl_dat <-read.table(text=getURL("https://raw.githubusercontent.com/jakelrich/MSU_WritingProject/master/Owls.txt"), header=T)
head(owl_dat)
###Exploratory Data Analysis

###Writing Simulation Function
##Overview:
#glmm.sim - a wrapper function to simulate data from a given glmm model object
#for a specified number of iterations and return a iter x mod matrix of model 
#scores, where iter is the number of iterations and mod is the number of model fitting
#methods used.
##Inputs:
#model - a previously defined glmer object used to simulate data (this will hopefully
#become a user defined function eventually.)
#iter - the number of simulation iterations to run. For each iteration, a new dataset is generated. 
##Outputs:
#scores - a matrix (list?) of model fitting "scores" for each model. 
###########
#Notes: simulate.merMod (lme4) will simulate from a glmer object.
###########
glmm.sim <- function(model, iter){
  
  for(i in 1:length(iter)){
    
  }
  
  return(scores)
}

###########################################################
###########################################################
#Scratch Pad
set.seed(54177)
###Likelihood-based Analysis using glmmPQL
library(MASS)
owl_dat$Calls <- owl_dat$SiblingNegotiation
owl_glmmPQL <- glmmPQL(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime, random = ~ 1|Nest,
                       data = owl_dat, family = poisson)
summary(owl_glmmPQL)
###Likelihood-based Analysis using glmer
library(lme4)
#Laplace Approximation
owl_glmer_LA <- glmer(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime + (1|Nest),
                      data = owl_dat, family = poisson, nAGQ = 1)
#Adaptive Gaussian-Hermite Quadrature
owl_glmer_AGHQ <- glmer(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime + (1|Nest),
                   data = owl_dat, family = poisson, nAGQ = 25)
summary(owl_glmer_LA)
summary(owl_glmer_AGHQ)
###Likelihood-based Analysis using glmmML
owl_glmmML <- glmmML(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime, cluster = Nest,
                     data = owl_dat, family = poisson, n.points = 25)
summary(owl_glmmML)

###Bayesian GLMM Analysis with JAGS
library(rjags)
#Setting up the data for JAGS
X <- model.matrix(owl_glmer)
K <- ncol(X)
Nre <- length(unique(owl_dat$Nest))
win.data1 <- list(Y = owl_dat$Calls, X = X, N = nrow(owl_dat), re = owl_dat$Nest,
                  b0 = rep(0, K), B0 = diag(0.0001, K), a0 = rep(0, Nre), A0 = diag(Nre))
#Writing JAGS Model Code
sink("GLMM.txt")
cat("
model { 
  #1. Diffuse priors for regression parameters.
  beta ~ dmnorm(b0[], B0[,])
  
  #Diffuse priors for randoms effect hive
  a ~ dmnorm(a0, tau.re*A0[,])
  num ~ dnorm(0, 0.0016)
  denom ~ dnorm(0,1)
  sigma.re <- abs(num/denom)
  tau.re <- 1 / (sigma.re*sigma.re)

  #2. Likelihood
  for (i in 1:N){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- eta[i]
    eta[i] <- inprod(beta[], X[i,]) + a[re[i]]
    
    #3. Discrepancy Measures
    YNew[i] ~ dpois(mu[i]) #New Data
    ExpY[i] <- mu[i]
    VarY[i] <- mu[i]
    PRes[i] <- (Y[i] - ExpY[i])/ sqrt(VarY[i])
    PResNew[i] <- (YNew[i] - ExpY[i]) / sqrt(VarY[i])
    D[i] <- pow(PRes[i], 2)
    DNew[i] <- pow(PResNew[i], 2)
  }

  Fit <- sum(D[1:N])
  FitNew <- sum(DNew[1:N])

}", fill = TRUE)
sink()

inits1 <- function() {
  list(beta = rnorm(K, 0, 0.01), a = rnorm(Nre, 0, 2), num = runif(1, 0, 25), denom = runif(1, 0, 1))
}

params1 <- c('beta', 'a', 'sigma.re', 'Fit', 'FitNew')

J0 <- jags.model(file = "GLMM.txt", data = win.data1, inits = inits1,
                 n.chains = 4, n.adapt = 30000)

mcmc.samples <- coda.samples(J0, params1, n.iter = 250000, thin = 25)

bayes_stats <- summary(mcmc.samples)$statistics
#First, using the basic Booth and Hobert dataset
#to fit a glmm with a logistic link, one variance component,
#one fixed effect, and an intercept of 0. The Monte Carlo
#sample size is 100 to save time.
library(glmm)
data(BoothHobert)
set.seed(1234)
mod.mcml1<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
                family.glmm=bernoulli.glmm,m=100,doPQL=TRUE)
nest_names <- as.character(unique(owl_dat$Nest))
var.eq <- c(rep(0,27))
owl_mc <- glmm(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime,
            random = list(~0 + Nest), varcomps.names = nest_names, data = owl_dat,
            family.glmm = poisson.glmm, varcomps.equal = var.eq, m = 100, debug = TRUE)
mod.mcml1$beta
mod.mcml1$nu
summary(mod.mcml1)
coef(mod.mcml1)
#Next, a model setting two variance components equal.
data(Booth2)
set.seed(1234)
mod.mcml3<-glmm(y~0+x1,list(y~0+z1,~0+z2),varcomps.names=c("z"),
                varcomps.equal=c(1,1), data=Booth2,family.glmm=bernoulli.glmm,
                m=100,doPQL=FALSE)
mod.mcml3$beta
mod.mcml3$nu
summary(mod.mcml3)
#Now, a model with crossed random effects. There are two distinct
#variance components. To get more accurate answers for this model,
#use a larger Monte Carlo sample size, such as m=10^4 or 10^5
#and doPQL=TRUE.
set.seed(1234)
data(salamander)
m<-10
sal<-glmm(Mate~0+Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"),
          data=salamander,family.glmm=bernoulli.glmm,m=m,debug=TRUE,doPQL=FALSE)
summary(sal)


library(mcemGLM)
ptm<-proc.time()
mcem <- mcemGLMM(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime, 
         random  = list(~0 + Nest), owl_dat, family = "poisson",
         vcDist = "normal")
proc.time() - ptm
mcemGLMMext(mcem, it = 20)
