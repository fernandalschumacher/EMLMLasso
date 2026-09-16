############################################################
# Code for an example of data application 
############################################################
#Packages
if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)
if(!require(mvtnorm)) install.packages("mvtnorm")
library(mvtnorm)
if(!require(splmm)) install.packages("splmm") # for comparison
library(splmm)
if(!require(glmmLasso)) install.packages("glmmLasso")  # for comparison
library(glmmLasso)
if(!require(glmnet)) install.packages("glmnet")
library(glmnet)

#sourcing main file
source("main-EMLMLasso.R")

# Simulating a data set
N      <- 200 # subjects
nj1    <- 5
p      <- 50
pstar  <- 10 # important variables 

set.seed(2022)
X       <- data.matrix(matrix(rnorm(N*nj1*p,0,1), ncol=p))
X          <- data.matrix(X)
beta      <- c( rep(1,pstar), rep(0, (p-pstar)) )

z       <- cbind(1,rep((seq(1:nj1)),N))
time    <- z[,2]
ID      <- as.factor(rep(1:N, each = nj1))

q       <- ncol(z)
nj      <- rep(nj1,N)

D       <- matrix(c(1,.25,.25,1),2,2)
sigmae  <- 1

y          <- matrix(0,N*nj1,1)
for (j in 1:N)
{
  x1       <- matrix(X[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
  z1       <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q)
  y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))] <- x1%*%beta + z1%*%t(rmvnorm(1,rep(0,q),D)) + t(rmvnorm(1,rep(0,nj[j]), sigmae*diag(nj[j])))
}

######################
#EMLMLasso
## using "lambda.1se"
est_EMLM <- EMLasso_cvLMM(y = y,x = X,z = z,nj = nj,folds = 20, lmethod = "lambda.1se",
                          precisao=0.00001, MaxIter=10000)
est_EMLM$beta
est_EMLM$lambda


######################
### REGULARIZATION BY glmmLasso
dat <- data.frame(y,X, ID, time)

## Search best lambda - long processing time
lambda_glmm <- seq(100,0,by=-1)
family <- gaussian(link = "identity")
BIC_glmm <- rep(Inf,length(lambda_glmm))
for(j in 1:length(lambda_glmm)){
  print(paste("Iteration of glmmLasso:", j,sep=""))
  glm1 <- try(glmmLasso(y~X-1, 
                        rnd = list(ID=~ time),  
                        family = family, data = dat, lambda=lambda_glmm[j], 
                        control=list(center=F)), silent=TRUE)
  if(!inherits(glm1, "try-error")){
    BIC_glmm[j]<-glm1$bic
  }
}
# final model
lambda_glmm[which.min(BIC_glmm)]
glm1 <- glmmLasso(y~X-1, rnd = list(ID=~ time),family = family, 
                  data = dat, lambda=lambda_glmm[which.min(BIC_glmm)], 
                      control=list(center=FALSE))
glm1$coefficients

######################
### REGULARIZATION BY splmm
## Search best lambda
lambda_splmm <- seq(0.1,0.5,0.004) #fixed effects
lam2 <- 0.1 #random effects
Lasso_splmm  <- list()
BIC_splmm <- matrix(Inf,length(lambda_splmm),1)
Lasso_splmm <- try(splmmTuning(x=data.frame(1,X), y=y, z=z, grp=ID, lam1.seq=lambda_splmm,
                              lam2.seq=lam2, penalty.b="lasso", 
                              penalty.L="lasso"), silent=TRUE)
if(!inherits(Lasso_splmm, "try-error")) {
  BIC_splmm  <- round(Lasso_splmm$BIC.lam1,3)
}
# final model
Lasso_splmm$best.fit$coefficients

#Lambda Estimates
#EMLMLasso
est_EMLM$lambda
#0.0653

#glmmLasso
lambda_glmm[which.min(BIC_glmm)]
#18

#splmm
lambda_splmm[which.min(BIC_splmm)]
#0.168
