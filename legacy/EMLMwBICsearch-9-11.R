########################################################################################
# revisions from March 2026 motivated by the lack of monotonicity in cv function
########################################################################################
library(glmnet)
library(mvtnorm)
library(inflection)
library(tidyverse)


logveropen <- function(beta, sigmae, D, lambda, y, x, z, nj)
{
  n     <- length(nj)
  N     <- sum(nj)
  p     <- dim(x)[2]
  q     <- dim(z)[2]
  suma1 <- 0
  
  for (j in 1:n)
  {
    y1    <- y[(sum(nj[1:j-1])+1):(sum(nj[1:j]))]
    x1    <- matrix(x[(sum(nj[1:j-1])+1):(sum(nj[1:j])),], ncol=p)
    z1    <- matrix(z[(sum(nj[1:j-1])+1):(sum(nj[1:j])),], ncol=q)
    med   <- x1 %*% beta
    Psi   <- z1 %*% D %*% t(z1) + sigmae * diag(nj[j])
    suma1 <- suma1 + dmvnorm(y1, med, Psi, log=TRUE)
  }
  
  # return(suma1 - lambda*sum(abs(beta)))
  return(suma1 - lambda*sum(abs(beta)))
}

##
EMLasso_LMM_testv5 <- function(y, x, z, nj, lambda, precisao=1e-6, MaxIter=200, 
                               penalty_factor=NULL, verbose=interactive())
{
  n <- length(nj); N <- sum(nj); p <- dim(x)[2]; q <- dim(z)[2]
  
  if(is.null(penalty_factor)) {
    penalty_factor <- rep(1, p)
  }
  
  # --- Initialization ---
  betaI <- glmnet(x, y, lambda=lambda, family="gaussian", intercept=FALSE, alpha=1, penalty.factor=penalty_factor)
  beta <- betaI$beta[,1]
  
  sigmae <- sum((y - x%*%beta)^2)/N
  sigmae <- max(sigmae, 1e-4) # Variance floor protection
  
  D <- diag(q)
  loglik <- logveropen(beta, sigmae, D, lambda, y, x, z, nj)
  loglikvec <- numeric(length = MaxIter)
  loglikvec[1] <- loglik
  teta <- c(beta,sigmae,D[upper.tri(D, diag=TRUE)])
  
  criterio <- 1; count <- 0
  
  end_indices <- cumsum(nj)
  start_indices <- c(1, end_indices[-n] + 1)
  
  while (criterio > precisao && count < MaxIter)
  {
    count <- count + 1
    sum3 <- 0; sum4 <- matrix(0,q,q); suma5 <- matrix(0,p,p)
    b <- matrix(0,q,n); zm <- matrix(0, nrow = N, ncol = 1)
    
    # --- E-Step ---
    for (j in 1:n)
    {
      idx <- start_indices[j]:end_indices[j]
      y1 <- y[idx]
      x1 <- x[idx, , drop = FALSE] #matrix(x[idx,], ncol=p)
      z1 <- z[idx, , drop = FALSE] #matrix(z[idx,], ncol=q)
      
      med <- x1%*%beta
      Psi <- z1%*%D%*%t(z1) + sigmae*diag(nj[j])
      IPsi <- solve(Psi, tol=1e-30)
      Lambda <- solve(solve(D)+t(z1)%*%z1/sigmae, tol=1e-30)
      
      b[,j] <- Lambda %*% t(z1) %*% (y1 - med) / sigmae
      y1m <- y1 - z1%*%b[,j]
      
      sum3 <- sum3 + t(y1m-med)%*%(y1m-med) + sum(diag(Lambda%*%t(z1)%*%z1))
      sum4 <- sum4 + Lambda + b[,j]%*%t(b[,j])
      suma5 <- suma5 + t(x1)%*%IPsi%*%x1
      zm[idx, 1] <- y1m
    }
    
    # --- M-Step ---
    la.eq <- glmnet(x, zm, lambda=lambda, family="gaussian", intercept=FALSE, alpha=1, thresh=1e-10, penalty.factor=penalty_factor)
    beta <- la.eq$beta[,1]
    
    calculated_sigmae <- as.numeric(sum3)/N
    sigmae <- max(calculated_sigmae, 1e-4) # Variance floor protection
    D <- sum4/n
    
    # Optional diagnostic printing
    if (verbose && count %% 10 == 0) {
      cat("  Iter:", count, "| Sigma_e2:", round(sigmae, 6), "| Lambda:", round(lambda, 6), "\n")
    }
    
    # --- Convergence ---
    loglik1 <- loglik
    loglik <- logveropen(beta,sigmae,D,lambda,y,x,z,nj)
    
    if (is.infinite(loglik) || is.nan(loglik)) break
    
    loglikvec[count] <- loglik
    teta1 <- teta
    teta <- c(beta,sigmae,D[upper.tri(D, diag=TRUE)])
    criterio <- sqrt(sum((teta1 - teta)^2))
    
    if (count == MaxIter) criterio <- precisao/10
  }
  
  teta <- c(beta,sigmae,D[upper.tri(D, diag=TRUE)])
  npar <- length(teta)-length(which(beta==0))
  AIC <- -2*loglik + 2*npar
  BIC <- -2*loglik + log(N)*npar
  
  obj <- list(beta=beta, sigmae=sigmae, D=D, b=b, suma5=suma5,
              lambda=lambda, AIC=AIC, BIC=BIC, iter=count, 
              loglik = loglikvec[1:count])
  class(obj) <- "EMLasso_LMM"
  return(obj)
}

EMLasso_BIC_Search <- function(y, x, z, nj, n_lambdas = 30, precisao = 0.0001, MaxIter = 200, 
                               penalty_factor = NULL, verbose = interactive()) {
  
  p <- dim(x)[2]; N<-nrow(y)
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, p)
  }
  
  # 1. Generate the Data-Driven Grid
  init_mod <- glmnet(x, y, family="gaussian", intercept=FALSE, alpha=1, nlambda=n_lambdas, penalty.factor=penalty_factor)
  lambda_grid <- init_mod$lambda#*(N/var(y))
  
  if (verbose) {
    cat("Generated", length(lambda_grid), "data-driven lambdas ranging from", 
        round(max(lambda_grid), 4), "to", round(min(lambda_grid), 4), "\n\n")
  }
  
  results <- list()
  bics <- numeric(length(lambda_grid))
  active_vars <- numeric(length(lambda_grid))
  
  # 2. Run your stable EMLasso over the grid
  for (i in seq_along(lambda_grid)) {
    current_lambda <- lambda_grid[i]
    
    if (verbose) cat("Fitting Model", i, "| Lambda:", round(current_lambda, 5), "\n")
    
    # Pass verbose argument down to the testv5 function
    fit <- EMLasso_LMM_testv5(y = y, x = x, z = z, nj = nj, 
                              lambda = current_lambda, 
                              precisao = precisao, 
                              MaxIter = MaxIter,
                              penalty_factor = penalty_factor,
                              verbose = verbose) 
    
    results[[i]] <- fit
    bics[i] <- fit$BIC
    active_vars[i] <- sum(fit$beta != 0)
    
    if (verbose) cat("--> Active betas:", active_vars[i], "| BIC:", round(bics[i], 2), "\n\n")
  }
  
  # 3. Select the best model
  best_idx <- which.min(bics)
  best_model <- results[[best_idx]]
  
  if (verbose) {
    cat("========================================\n")
    cat("Optimal Model Found at Index", best_idx, "\n")
    cat("Best Lambda:", lambda_grid[best_idx], "\n")
    cat("Active Variables:", active_vars[best_idx], "\n")
    cat("Lowest BIC:", bics[best_idx], "\n")
  }
  
  return(list(
    best_model = best_model,
    best_lambda = lambda_grid[best_idx],
    lambda_grid = lambda_grid,
    bics = bics,
    active_vars = active_vars,
    all_models = results
  ))
}

