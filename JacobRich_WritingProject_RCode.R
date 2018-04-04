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
list_of_packages <- c("lme4","MCMCglmm","coda","glmm","ggplot2","RCurl", "simr","devtools",
                      "rjags","mcemGLM","MASS","glmmML","GLMMmisc","nlme")
new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (!require("GLMMmisc")) devtools::install_github("pcdjohnson/GLMMmisc")
if(length(new.packages)) install.packages(new.packages)
lapply(list_of_packages, library, character.only = TRUE)
#Installing GLMMmisc from Github. GLMMmisc is a package of GLMM diagnostic tools
#and other functions that are useful with GLMMs. 

#Loading in dataset(s) and cleaning
owl_dat <-read.table(text=getURL("https://raw.githubusercontent.com/jakelrich/MSU_WritingProject/master/Owls.txt"), header=T)
owl_dat$Calls <- owl_dat$SiblingNegotiation
#Saving a copy of the dataset for simulation
sim_dat <- owl_dat
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
glmm.sim <- function(model_object, iter, seed){
  ###Initializing storage matrices
  beta_ests <- matrix(0, iter, 12)
  colnames(beta_ests) <- c("PQL_B1", "PQL_B2", "PQL_B3", "LA_B1", "LA_B2", "LA_B3",
                           "AGHQ_B1", "AGHQ_B2", "AGHQ_B3", "JAGS_B1", "JAGS_B2", "JAGS_B3")
  rand_ests <- matrix(0, iter, 108)
  colnames(rand_ests) <- c("PQL_1", "PQL_2", "PQL_3", "PQL_4", "PQL_5", "PQL_6", "PQL_7",
                           "PQL_8", "PQL_9", "PQL_10", "PQL_11", "PQL_12", "PQL_13", "PQL_14",
                           "PQL_15", "PQL_16", "PQL_17", "PQL_18", "PQL_19", "PQL_20", "PQL_21",
                           "PQL_22", "PQL_23", "PQL_24", "PQL_25", "PQL_26", "PQL_27",
                           "LA_1", "LA_2", "LA_3", "LA_4", "LA_5", "LA_6", "LA_7",
                           "LA_8", "LA_9", "LA_10", "LA_11", "LA_12", "LA_13", "LA_14",
                           "LA_15", "LA_16", "LA_17", "LA_18", "LA_19", "LA_20", "LA_21",
                           "LA_22", "LA_23", "LA_24", "LA_25", "LA_26", "LA_27", 
                           "AGHQ_1", "AGHQ_2", "AGHQ_3", "AGHQ_4", "AGHQ_5", "AGHQ_6", "AGHQ_7",
                           "AGHQ_8", "AGHQ_9", "AGHQ_10", "AGHQ_11", "AGHQ_12", "AGHQ_13", "AGHQ_14",
                           "AGHQ_15", "AGHQ_16", "AGHQ_17", "AGHQ_18", "AGHQ_19", "AGHQ_20", "AGHQ_21",
                           "AGHQ_22", "AGHQ_23", "AGHQ_24", "AGHQ_25", "AGHQ_26", "AGHQ_27",
                           "JAGS_1", "JAGS_2", "JAGS_3", "JAGS_4", "JAGS_5", "JAGS_6", "JAGS_7",
                           "JAGS_8", "JAGS_9", "JAGS_10", "JAGS_11", "JAGS_12", "JAGS_13", "JAGS_14",
                           "JAGS_15", "JAGS_16", "JAGS_17", "JAGS_18", "JAGS_19", "JAGS_20", "JAGS_21",
                           "JAGS_22", "JAGS_23", "JAGS_24", "JAGS_25", "JAGS_26", "JAGS_27")
  sigma_re_ests <- matrix(0, iter, 4)
  colnames(sigma_re_ests) <- c("PQL", "LA", "AGHQ", "JAGS")
  
  ###Writing JAGS Model Code
  sink("GLMM.txt")
  cat("
      model { 
      #1. Diffuse priors for regression parameters.
      beta ~ dmnorm(b0[], B0[,])
      
      #Diffuse priors for randoms effect hive
      a ~ dmnorm(a0, tau.re*A0[,])
      num ~ dnorm(0, 0.0016)
      denom ~ dnorm(0,1)
      sigma.re <- abs(num/(denom*denom))
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
    list(beta = rnorm(K, 0, 1), a = rnorm(Nre, 0, 2), num = runif(1, 0, 25), denom = runif(1, 0, 1))
  }
        
  params1 <- c('beta', 'a', 'sigma.re', 'Fit', 'FitNew')
  
  ###Running simulations
  for(i in 1:length(iter)){
    
    ###Setting simulation seed
    set.seed(seed+iter)
    
    ###Generating new responses
    #use simr 
    sim_dat$Calls <- doSim(model_object)
    
    ###Fitting Classical Models
    #Penalized Quasi-likelihood
    owl_glmmPQL <- glmmPQL(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime, random = ~ 1|Nest,
                           data = sim_dat, family = poisson, verbose=FALSE)
    #Laplace Approximation
    owl_glmer_LA <- glmer(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime + (1|Nest),
                          data = sim_dat, family = poisson, nAGQ = 1)
    #Adaptive Gaussian-Hermite Quadrature
    owl_glmer_AGHQ <- glmer(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime + (1|Nest),
                            data = sim_dat, family = poisson, nAGQ = 25)
    #Storing summaries for later use
    PQL <- summary(owl_glmmPQL)
    LA <- summary(owl_glmer_LA)
    AGHQ <- summary(owl_glmer_AGHQ)
    
    ###Fitting Bayesian Model with JAGS
    X <- model.matrix(owl_glmer_AGHQ)
    K <- ncol(X)
    Nre <- length(unique(sim_dat$Nest))
    win.data1 <- list(Y = sim_dat$Calls, X = X, N = nrow(sim_dat), re = sim_dat$Nest,
                      b0 = rep(0, K), B0 = diag(0.0001, K), a0 = rep(0, Nre), A0 = diag(Nre))
    J0 <- jags.model(file = "GLMM.txt", data = win.data1, inits = inits1,
                     n.chains = 4, n.adapt = 125000)
    mcmc.samples <- coda.samples(J0, params1, n.iter = 75000, thin = 5)
    bayes_stats <- summary(mcmc.samples)$statistics
        
    ###Storing beta estimates
    beta_ests[i,1:3] <- PQL$coefficients$fixed
    beta_ests[i,4:6] <- LA$coefficients[,1]
    beta_ests[i,7:9] <- AGHQ$coefficients[,1]
    beta_ests[i,10:12] <- bayes_stats[30:32, 1]
    
    ###Storing random intercept estimates
    rand_ests[i,1:27] <- c(t(PQL$coefficients$random$Nest))
    rand_ests[i,28:54] <- t(coef(owl_glmer_LA)$Nest['(Intercept)'])
    rand_ests[i,55:81] <- t(coef(owl_glmer_AGHQ)$Nest['(Intercept)'])
    rand_ests[i,82:108] <- bayes_stats[3:29, 1]
    
    ###Storing random effects variance estimates
    sigma_re_ests[i,1] <- sqrt(c(getVarCov(owl_glmmPQL)))
    sigma_re_ests[i,2] <- sqrt(c(LA$varcor$Nest))
    sigma_re_ests[i,3] <- sqrt(c(AGHQ$varcor$Nest))
    sigma_re_ests[i,4] <- bayes_stats[33,1]
    
    ###(Hopefully) print the iteration of simulation
    print(iter)
  }
  
  ###Saving as a list for returning    
  sim_results <- list(beta_ests, rand_ests, sigma_re_ests)
  return(sim_results)
}
###Implementing Simulation
param_ests <- glmer(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime + (1|Nest),
                    data = owl_dat, family = poisson, nAGQ = 25)
system.time(test <- glmm.sim(param_ests, 1, 54177))
simulations <- test
###Saving simulations as RData file (just in case).
saveRDS(simulations, file = "D:/Google Drive/18 - Writing Project/sim.RData")





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
c(getVarCov(owl_glmmPQL))
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

PQL <- summary(owl_glmmPQL)
AGHQ <- summary(owl_glmer_AGHQ)
LA <- summary(owl_glmer_LA)

beta_ests <- matrix(0, 10, 12)
colnames(beta_ests) <- c("PQL_B1", "PQL_B2", "PQL_B3", "LA_B1", "LA_B2", "LA_B3",
                 "AGHQ_B1", "AGHQ_B2", "AGHQ_B3", "JAGS_B1", "JAGS_B2", "JAGS_B3")
beta_ests[1,1:3] <- AGHQ$coefficients[,1]
beta_ests[1,4:6] <- LA$coefficients[,1]
beta_ests[1,7:9] <- PQL$coefficients$fixed
beta_ests[1,10:12] <- bayes_stats[30:32, 1]

###Likelihood-based Analysis using glmmML
owl_glmmML <- glmmML(Calls ~ offset(BroodSize) + FoodTreatment + ArrivalTime, cluster = Nest,
                     data = owl_dat, family = poisson, n.points = 25)
summary(owl_glmmML)

###Bayesian GLMM Analysis with JAGS
library(rjags)
#Setting up the data for JAGS
X <- model.matrix(owl_glmer_AGHQ)
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
  sigma.re <- abs(num/(denom*denom))
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
  list(beta = rnorm(K, 0, 1), a = rnorm(Nre, 0, 2), num = runif(1, 0, 25), denom = runif(1, 0, 1))
}

params1 <- c('beta', 'a', 'sigma.re', 'Fit', 'FitNew')

J0 <- jags.model(file = "GLMM.txt", data = win.data1, inits = inits1,
                 n.chains = 4, n.adapt = 125000)

mcmc.samples <- coda.samples(J0, params1, n.iter = 75000, thin = 5)

plot(mcmc.samples)
summary(mcmc.samples)
bayes_stats <- summary(mcmc.samples)$statistics
bayes_stats[3:29]

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
