############################################################
# Code for cholesterol data application 
############################################################
#Packages
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(mvtnorm)) install.packages("mvtnorm")
library(mvtnorm)
if(!require(splmm)) install.packages("splmm")
library(splmm)
if(!require(glmmLasso)) install.packages("glmmLasso")
library(glmmLasso)
if(!require(glmnet)) install.packages("glmnet")
library(glmnet)


if(!require(qrLMM)) install.packages("qrLMM")

#sourcing main file
source("main-EMLMLasso.R")

#Data
data(Cholesterol, package = "qrLMM")
Cholesterol <- Cholesterol %>% arrange(newid)
nj <- Cholesterol %>% count(newid) %>% select(n) %>% unlist()

#tratat   <- as.vector(Cholesterol$newid)
yorig    <- as.vector(Cholesterol$cholst/100)
year     <- as.vector((Cholesterol$year-5)/10)
yearf    <- as.factor(Cholesterol$year)
ID       <- as.factor(Cholesterol$newid)
D        <- matrix(c(1,0.5,0.5,1),2,2)
q        <- 2

#V5: simulated covariate of the Bernoulli distribution with p=0.5
set.seed(36)
Cholesterol$V5   <- rep(rbinom(length(nj),1,0.5), nj)

#V7 and V8: simulated covariate of the bivariate normal distribution 
#with mean=c(0,0) and Sigma=matrix(c(1,0.5,0.5,1),2,2)
set.seed(36)
V78           <- rmvnorm(length(nj),rep(0,q),D)
Cholesterol$V7   <- rep(V78[,1],nj)
Cholesterol$V8   <- rep(V78[,2],nj)

V5            <- Cholesterol$V5

V7            <- Cholesterol$V7

V8            <- Cholesterol$V8

sex           <- as.vector(Cholesterol$sex)

age           <- as.vector(Cholesterol$age)

Z             <- cbind(matrix(1,length(yorig),1),year)

#####################################

## For glmmLasso and EMLMLasso - matrix X with randomized covariates
#C1=Sex, C2=Bivariate Normal_2,
#C3=Bivariate Normal_1, C4=Time, C5=Age, C6=Bernoulli(0.5)
#C7=Sex*Age, C8=Sex*Age*Time, C9=Sex*Time, C10=Age*Time 
X               <- cbind(sex, scale(V8), scale(V7), scale(year),
                         scale(age), V5, scale(sex*age), scale(sex*age*year),
                         scale(sex*year), scale(age*year))
colnames(X) <- c("sex", "V8", "V7", "year", "age", "V5", "sex:age",
                 "sex:age:year", "sex:year", "age:year")

p1              <- ncol(X) #no intercept

y               <- scale (yorig)

## For splmm - matrix x with randomized covariates
#C1=1, C2=Sex, C3=Bivariate Normal_2,
#C4=Bivariate Normal_1, C5=Time, C6=Age, C7=Bernoulli(0.5)
#C8=Sex*Age, C9=Sex*Age*Time, C10=Sex*Time, C11=Age*Time 
x               <- cbind(1, sex, scale(V8), scale(V7), scale(year),
                         scale(age), V5, scale(sex*age), scale(sex*age*year),
                         scale(sex*year), scale(age*year))
colnames(x) <- c("(Intercept)", "sex", "V8", "V7", "year", "age", "V5", 
                 "sex:age", "sex:age:year", "sex:year", "age:year")

p               <- ncol(x) #with intercept


######################
### REGULARIZATION BY EMLMLasso
## Search best lambda
lambda_EMLM <- seq(0.001,0.5,length.out=500)
BIC_EMLM <- matrix(0,length(lambda_EMLM),1)
est <- list()
for (k in seq_along(lambda_EMLM))
{
  print(paste("Iteration of EMLMLasso:", k,sep=""))
  est[[k]] <- EMLasso_LMM(y,X,Z,nj,lambda=lambda_EMLM[k],precisao=0.00001,MaxIter=1000)
  BIC_EMLM[k]  <- est[[k]]$BIC
}

######################
### REGULARIZATION BY glmmLasso
## Search best lambda
lambda_glmm <- seq(500,0,by=-1)
family <- gaussian(link = "identity")
BIC_glmm <- rep(Inf,length(lambda_glmm))
for(j in 1:length(lambda_glmm)){
  print(paste("Iteration of glmmLasso:", j,sep=""))
  glm1 <- try(glmmLasso(y~X-1, 
                        rnd = list(ID=~ year),  
                        family = family, data = Cholesterol, lambda=lambda_glmm[j], 
                        control=list(center=FALSE)), silent=TRUE)
  if(!inherits(glm1, "try-error")){
    BIC_glmm[j]<-glm1$bic
  }
}

######################
### REGULARIZATION BY splmm
## Search best lambda
lambda_splmm <- seq(0.1,0.5,0.004) #fixed effects
lam2 <- 0.1 #random effects
Lasso_splmm  <- list()
BIC_splmm <- matrix(Inf,length(lambda_splmm),1)
Lasso_splmm <- try(splmmTuning(x=x, y=y, z=Z, grp=ID, lam1.seq=lambda_splmm,
                              lam2.seq=lam2, penalty.b="lasso", 
                              penalty.L="lasso"), silent=TRUE)
if(!inherits(Lasso_splmm, "try-error")) {
  BIC_splmm  <- round(Lasso_splmm$BIC.lam1,3)
}

#Lambda Estimate
#EMLMLasso
lambda_EMLM[which.min(BIC_EMLM)]
#0.022

#glmmLasso
lambda_glmm[which.min(BIC_glmm)]
#18

#splmm
lambda_splmm[which.min(BIC_splmm)]
#0.156 with y=scale(yorig)

bic_data <- bind_rows(data.frame(lambda = lambda_EMLM, bic = BIC_EMLM, method = "EMLM"),
                      data.frame(lambda = lambda_glmm, bic = BIC_glmm, method = "glmm"),
                      data.frame(lambda = lambda_splmm, bic = BIC_splmm, method = "splmm"))
bic_data %>% group_by(method) %>% filter(bic==min(bic))

bic_data %>% ggplot(aes(x = lambda, y = bic)) + geom_line() + 
  facet_wrap(~method, scales = "free") + theme_bw()

### final models

#EMLMLasso
est_EMLM <- EMLasso_LMM(y,X,Z,nj,
                   lambda=lambda_EMLM[which.min(BIC_EMLM)],
                   precisao=0.00001, MaxIter=1000)
beta_EMLM  <- est_EMLM$beta
round(beta_EMLM,4)

#glmmLasso
est_glmm <- glmmLasso(y~X-1, 
                  rnd = list( ID= ~ year),  
                  data = Cholesterol, lambda=lambda_glmm[which.min(BIC_glmm)], 
                  family = gaussian(link="identity"),
                  final.re=TRUE,
                  control=list(center=FALSE, epsilon=5e-6))
round(summary(est_glmm)$coefficients[,1],4)


#splmm
est_splmm <- splmm(x=x, y=y, z=Z, grp=ID, lam1=lambda_splmm[which.min(BIC_splmm)],
             lam2=0.1, penalty.b="lasso", penalty.L="lasso")
round(est_splmm$fixef,4)
est_splmm %>% summary()
