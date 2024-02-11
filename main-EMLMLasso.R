
logvero <- function(beta,sigmae,D,y,x,z,nj)
{
  
  ################################################################################
  ## beta: regression parameters
  ## sigmae: error variance
  ## D: var-cov matrix of random effects 
  ## y: response vector
  ## x:  fixed effects design matrix
  ## z:  random effects design matrix
  ## nj: vector with number of observations per subject/cluster
  ################################################################################
  
  n     <- length(nj)
  N     <- sum(nj)
  p     <- dim(x)[2]
  q     <- dim(z)[2]
  suma1 <- 0
  
  for (j in 1:n)
  {
    y1    <- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
    x1    <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
    z1    <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q)
    med   <- x1%*%beta
    Psi   <- (z1)%*%(D)%*%t(z1) + sigmae*diag(nj[j])
    suma1 <- suma1 + dmvnorm(y1,med,Psi,log=TRUE)
  }
  
  return(suma1)  ## marginal log-likelihood
}

################################################################################
### IMPLEMENTACAO DO ALGORITMO EM: CASO NORMAL 
################################################################################

EM_LMM <- function(y,x,z,nj,precisao=0.000001,MaxIter=200)
{
  
  ################################################################################
  ## y: response vector
  ## x:  fixed effects design matrix
  ## z:  random effects design matrix
  ## nj: vector with number of observations per subject/cluster
  ################################################################################
  
  n        <- length(nj)
  N        <- sum(nj)
  p        <- dim(x)[2]
  q        <- dim(z)[2]
  
  #### Initial values
  beta     <- solve(t(x)%*%x)%*%t(x)%*%y
  sigmae   <- sum((y - x%*%beta)^2)/(n - p)
  D        <- diag(q)
  loglik   <- logvero(beta,sigmae,D,y,x,z,nj)
  
  ################################################################################
  criterio <- 1
  count    <- 0
  
  while(criterio > precisao)
  {
    count  <- count + 1
    sum1   <- matrix(0,nrow=p,ncol=p)
    sum2   <- matrix(0,nrow=p,ncol=1)
    sum3   <- 0
    sum4   <- matrix(0,nrow=q,ncol=q)
    #suma5  <- matrix(0,nrow=p,ncol=p)
    b      <- matrix(0,nrow=q,ncol=n)
    
    for (j in 1:n )
    {
      y1     <- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      x1     <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1     <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q)
      med    <- x1%*%beta
      Psi    <- (z1)%*%(D)%*%t(z1) + sigmae*diag(nj[j])
      
      mediab <- D%*%t(z1)%*%solve(Psi)%*%(y1 - med)
      Lambda <- solve(solve(D)+t(z1)%*%z1/sigmae)
      
      sum1   <- sum1+t(x1)%*%x1
      sum2   <- sum2+(t(x1)%*%(y1-z1%*%b[,j]))
      sum3   <- sum3+t(y1-med-z1%*%b[,j])%*%(y1-med-z1%*%b[,j])+sum(diag(Lambda%*%t(z1)%*%z1))
      sum4   <- sum4+Lambda+b[,j]%*%t(b[,j])
      #suma5  <- suma5+ (t(x1)%*%solve(Psi)%*%x1)
      b[,j]  <- Lambda%*%t(z1)%*%(y1-med)/sigmae
    }
    
    beta     <- solve(sum1)%*%sum2
    sigmae   <- as.numeric(sum3)/N
    D        <- sum4/n
    loglik1  <- loglik
    loglik   <- logvero(beta,sigmae,D,y,x,z,nj)
    criterio <- abs(1-loglik1/loglik)
    
    if  (count == MaxIter) {criterio=precisao/10}
  }
  
  print(count)
  teta <- c(beta,sigmae,D[upper.tri(D, diag = T)])
  npar <-length(c(teta))
  
  ### selection criteria:
  
  AIC            <- -2*loglik + 2*npar
  BIC            <- -2*loglik + log(N)*npar
  obj.out        <- list(beta = beta, sigmae= sigmae, D = D, b=b,
                         AIC=AIC,BIC=BIC, iter = count)
   class(obj.out) <- "EM_LMM"
  return(obj.out)
}

################################################################################
### IMPLEMENTACAO DO ALGORITMO EM  - LASSO
################################################################################

EMLasso_LMM <- function(y,x,z,nj,lambda, precisao=0.000001,MaxIter=200)
{
  
  ################################################################################
  ## y: response vector
  ## x:  fixed effects design matrix
  ## z:  random effects design matrix
  ## nj: vector with number of observations per subject/cluster
  ################################################################################
  
  n      <- length(nj)
  N      <- sum(nj)
  p      <- dim(x)[2]
  q      <- dim(z)[2]
  
  #### Initial values
  betaI  <- glmnet(x, y, lambda=lambda, family="gaussian", intercept = F, alpha=1)
  beta   <- betaI$beta[,1] 
  sigmae <- sum((y-x%*%beta)^2)/n
  D      <- diag(q)
  
  ################################################################################
  
  criterio <- 1
  count    <- 0
  
  loglik   <- logvero(beta,sigmae,D,y,x,z,nj)
  
  while(criterio > precisao)
  {
    count <- count + 1
    sum3  <- 0
    sum4  <- matrix(0,nrow=q,ncol=q)
    #suma5 <- matrix(0,nrow=p,ncol=p)
    b     <- matrix(0,nrow=q,ncol=n)
    zm    <- NULL
    yhat  <- matrix(0,nrow=N)
    
    for (j in 1:n)
    {
      y1      <- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]
      x1      <- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1      <- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q)
      med     <- x1%*%beta
      Psi     <- (z1)%*%(D)%*%t(z1) + sigmae*diag(nj[j])
      IPsi    <- solve(Psi,tol=1e-30)
      mediab  <- D%*%t(z1)%*%IPsi%*%(y1-med)
      Lambda  <- solve(solve(D)+t(z1)%*%z1/sigmae,tol=1e-30)
      b[,j]   <- Lambda%*%t(z1)%*%(y1-med)/sigmae
      y1m     <- y1 - z1%*%b[,j]
      yhat[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))]     <- med + z1%*%b[,j]
      #sum2   <- sum2 + (t(x1)%*%(y1-z1%*%b[,j]))
      sum3    <- sum3 + t(y1m-med)%*%(y1m-med) + sum(diag(Lambda%*%t(z1)%*%z1))
      sum4    <- sum4 + Lambda + b[,j]%*%t(b[,j])
      #suma5   <- suma5+(t(x1)%*%IPsi%*%x1)
      zm      <- rbind(zm,y1m)
    }
    
    #beta     <- solve(sum1)%*%sum2
    lambda1   <- lambda*sigmae
    la.eq     <- glmnet(x, zm, lambda=lambda1, family="gaussian", intercept = F, alpha=1)
    beta      <- la.eq$beta[,1]  
    sigmae    <- as.numeric(sum3)/N
    D         <- sum4/n
    
    loglik1   <- loglik
    loglik    <- logvero(beta,sigmae,D,y,x,z,nj)
    criterio  <- abs(1-loglik1/loglik)
    
    if  (count == MaxIter) {criterio=precisao/10}
  }
  
  print(count)
  teta           <- c(beta,sigmae,D[upper.tri(D, diag = T)])
  npar           <- length(c(teta))-length(which(beta == 0))
  
  AIC            <- -2*loglik +2*npar
  BIC            <- -2*loglik +log(N)*npar
  obj.out        <- list(beta = beta, sigmae= sigmae, D = D, b=b,
                         AIC=AIC, BIC=BIC, yhat=yhat, iter = count)
  class(obj.out) <- "EMLasso_LMM"
  return(obj.out)
}
