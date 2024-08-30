
source("~/Desktop/Code/Simulation_Preparation.R")
source("~/Desktop/Code/Simulation_Thinning.R")

# Example setup
# Initialize parameters and Clayton copula parameter
lambda <- c(1.2, 1)   # Baseline intensities for both dimensions
alpha <- matrix(c(0.1, 0,   # Excitation effects from dim 1 to 1 and 2
                  0, 0.1),  # Excitation effects from dim 2 to 1 and 2
                nrow = 2, byrow = TRUE)
beta <- c(1.5, 1.5)     # Decay rates for both dimensions
copula_parameter <- 4
copula <- claytonCopula(param = copula_parameter, dim = 2)
theta <- list(lambda, alpha, beta)
theta_vector <- unlist(theta)

# Example parameters for the Hawkes process
lambda_start <- c(1.5, 1.5)   # Baseline intensities for both dimensions
alpha_start <- matrix(c(0.5, 0,   # Excitation effects from dim 1 to 1 and 2
                        0, 0.5),  # Excitation effects from dim 2 to 1 and 2
                      nrow = 2, byrow = TRUE)
beta_start <- c(2, 2)     # Decay rates for both dimensions
copula_parameter_start <- 2
copula <- claytonCopula(param = copula_parameter, dim = 2)

# Initialize vectors to store the L-infinity distances
Distance_copula <- numeric(0)
Distance_marginal <- numeric(0)

for (T in 10:15) {
  # Simulate the process until time T using copulas
  list <- simulate_until_T_copulas(T, theta, copula)
  times <- list[[1]]
  ids <- list[[2]]
  
  # MLE for frugal Hawkes processes with copula
  theta_start_copula <- list(lambda_start, alpha_start, beta_start, copula_parameter_start)
  result_copula <- frugal_mutual_exp_mle(times, ids, T, theta_start_copula, copula)
  
  # Compute the L-infinity distance for MLE with copula
  theta_mle_vector_copula <- unlist(result_copula$theta_mle)
  theta_vector_copula <- c(rep(1, length(theta_vector) - 1), theta_mle_vector_copula[length(theta_mle_vector_copula)])
  L2_distance_copula <- sqrt(sum((theta_mle_vector_copula - theta_vector_copula)^2))
  Distance_copula <- c(Distance_copula, L2_distance_copula)
  
  # MLE for standard Hawkes processes
  theta_start <- list(lambda_start, alpha_start, beta_start)
  result <- mutual_exp_mle(times, ids, T, theta_start)
  
  # Compute the L-infinity distance for MLE without copula
  theta_mle_vector <- unlist(result$theta_mle)
  L2_distance <- sqrt(sum((theta_mle_vector - theta_vector)^2))
  Distance_marginal <- c(Distance_marginal, L2_distance)
}

Distance_copula 

Distance_marginal





