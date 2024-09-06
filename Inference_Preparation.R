
# List of required libraries
libraries <- c("copula", "Matrix", "stats", "ggplot2", "dplyr", "gridExtra", "cubature", "numDeriv", 
               "knitr", "kableExtra")

# Loop to check and install missing libraries
for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib, dependencies = TRUE)
  }
  library(lib, character.only = TRUE)
}

# Flatten the parameters of the multivariate Hawkes processes with exp decay
flatten_theta <- function(theta) {
  return(c(theta[[1]], as.vector(theta[[2]]), theta[[3]]))
}

# Unflatten the vector of the multivariate Hawkes processes with exp decay
unflatten_theta <- function(theta_flat, m) {
  lambda <- theta_flat[1:m]
  alpha <- matrix(theta_flat[(m + 1):(m + m^2)], nrow = m, ncol = m)
  beta <- theta_flat[(m + m^2 + 1):(m + m^2 + m)]
  return(list(lambda, alpha, beta))
}

# Maximum Likelihood Estimation for the multivariate Hawkes processes with exp decay
mutual_exp_mle <- function(times, ids, T, theta_start) {
  m <- length(theta_start[[1]])
  theta_start_flat <- flatten_theta(theta_start)
  
  # Define the loss function (negative log-likelihood)
  # Define the loss function (negative log-likelihood)
  loss <- function(theta_flat) {
    theta <- unflatten_theta(theta_flat, m)
    ll <- -mutual_exp_log_likelihood(times, ids, T, theta)
    
    if (!is.finite(ll)) {
      return(1e10)  # Return a large value if the log-likelihood is not finite
    }
    return(ll)
  }
  
  # Set lower bounds to ensure positivity
  lower_bounds <- rep(0, length(theta_start_flat))
  upper_bounds <- c(10, 10, 10, 1e-8, 1e-8, 10, 10, 10, 10)
  
  # Minimize the negative log-likelihood using the L-BFGS-B method with bounds
  res <- optim(par = theta_start_flat, fn = loss, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)
  
  theta_mle <- unflatten_theta(res$par, m)
  logLik <- -res$value
  
  return(list(theta_mle = theta_mle, logLik = logLik))
}

# Flatten the parameters of the frugal Hawkes process into a single vector
flatten_theta_copula <- function(theta) {
  return(c(theta[[1]], as.vector(theta[[2]]), theta[[3]], theta[[4]]))
}

# Unflatten the vector of the frugal Hawkes process back into the original parameters
unflatten_theta_copula <- function(theta_flat, m) {
  lambda <- theta_flat[1:m]
  alpha <- matrix(theta_flat[(m + 1):(m + m^2)], nrow = m, ncol = m)
  beta <- theta_flat[(m + m^2 + 1):(m + m^2 + m)]
  copula_parameter <- theta_flat[(m + m^2 + m + 1)]
  return(list(lambda, alpha, beta, copula_parameter))
}

# Maximum Likelihood Estimation for the frugal Hawkes process with exponential decay
frugal_mutual_exp_mle <- function(times, ids, T, theta_start, copula) {
  m <- length(theta_start[[1]])
  theta_start_flat <- flatten_theta_copula(theta_start)
  
  # Define the loss function (negative log-likelihood)
  loss <- function(theta_flat) {
    theta <- unflatten_theta_copula(theta_flat, m)
    theta[[4]] <- max(theta[[4]], 1e-2)
    copula <- claytonCopula(param = theta[[4]], dim = 2)
    ll <- -frugal_mutual_log_likelihood(times, ids, T, theta, copula, theta[[4]])
    
    if (!is.finite(ll)) {
      return(1e10)  # Return a large value if the log-likelihood is not finite
    }
    return(ll)
  }
  
  # Set lower bounds to ensure positivity
  lower_bounds <- rep(0, length(theta_start_flat))
  upper_bounds <- c(10, 10, 10, 1e-8, 1e-8, 10, 10, 10, 10, 10)
  
  # Minimize the negative log-likelihood using the L-BFGS-B method with bounds
  res <- optim(par = theta_start_flat, fn = loss, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds)
  
  theta_mle <- unflatten_theta_copula(res$par, m)
  logLik <- -res$value
  
  return(list(theta_mle = theta_mle, logLik = logLik))
}

# 2-stage algorithm for MLE for the frugal Hawkes process with exponential decay
frugal_mutual_exp_mle_2stage <- function(times, ids, T, theta_start, copula_param_start) {
  m <- length(theta_start[[1]])
  
  # Flatten the initial values for lambda, alpha, beta for optimization
  theta_start_flat <- flatten_theta(theta_start)
  
  # Step 1: Optimize lambda, alpha, beta with the partial sum over log mutual_exp_hawkes_intensity
  loss_stage1 <- function(theta_flat) {
    theta <- unflatten_theta(theta_flat, m)
    ll <- -mutual_exp_log_likelihood(times, ids, T, theta)
    
    if (!is.finite(ll)) {
      return(1e10)  # Return a large value if the log-likelihood is not finite
    }
    return(ll)
  }

  # Optimize lambda, alpha, beta
  res_stage1 <- optim(par = theta_start_flat, fn = loss_stage1, method = "L-BFGS-B", lower = rep(0, length(theta_start_flat)), upper = c(10, 10, 10, 1e-8, 1e-8, 10, 10, 10) )
  optimized_theta_stage1 <- unflatten_theta(res_stage1$par, m)

  # Step 2: Optimize copula parameter using the optimized alpha, beta, gamma
  loss_stage2 <- function(copula_param) {
    theta_new <- optimized_theta_stage1
    copula_new <- claytonCopula(param = copula_param, dim = 2)
    ll <- -frugal_mutual_log_likelihood(times, ids, T, theta_new, copula_new, copula_param) - log(copula_param)*R/10
    
    if (!is.finite(ll)) {
      return(1e10)  # Return a large value if the log-likelihood is not finite
    }
    return(ll)
  }
  
  res_stage2 <- optim(par = copula_param_start, fn = loss_stage2, method = "L-BFGS-B", lower = 0.1, upper = 10)
  optimized_copula_param <- res_stage2$par
  
  # Step 3: Re-optimize lambda, alpha, beta with the optimized copula parameter
  loss_final <- function(theta_flat) {
    theta <- unflatten_theta(theta_flat, m)
    copula_optimized <- claytonCopula(param = optimized_copula_param, dim = 2)
    ll <- -frugal_mutual_log_likelihood(times, ids, T, theta, copula_optimized, optimized_copula_param) 
    
    if (!is.finite(ll)) {
      return(1e10)  # Return a large value if the log-likelihood is not finite
    }
    
    return(ll)
  }
  
  res_final <- optim(par = flatten_theta(optimized_theta_stage1), fn = loss_final, method = "L-BFGS-B", lower = rep(0, length(theta_start_flat)), upper = c(10, 10, 10, 1e-8, 1e-8, 10, 10, 10) )         
  
  final_theta <- unflatten_theta(res_final$par, m)
  final_param <- c(final_theta, optimized_copula_param)
  
  final_logLik <- -res_final$value
  
  return(list(theta_mle = final_param, logLik = final_logLik))
}







