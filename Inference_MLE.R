
source("Simulation_Preparation.R")
source("Simulation_Methods.R")
source("Inference_Preparation.R")

# Example usage
{
  lambda <- c(0.2, 0.23)   # Baseline intensities for both dimensions
  alpha <- matrix(c(0.32, 0,   # Excitation effects from dim 1 to 1 and 2
                    0, 0.54),  # Excitation effects from dim 2 to 1 and 2
                  nrow = 2, byrow = TRUE)
  beta <- c(4.9, 4.8)     # Decay rates for both dimensions
  copula_parameter <- 1.4
  copula <- claytonCopula(param = copula_parameter, dim = 2)
  theta <- list(lambda, alpha, beta)
  
  theta <- list(lambda, alpha, beta)
  theta_vector <- unlist(theta)
}

# Parameters for the starting point of the MLE
{
  lambda_start <- c(1.5, 1.5)
  alpha_start <- matrix(c(0.5, 0, 
                          0, 0.5), nrow = 2, byrow = TRUE)
  beta_start <- c(2, 2)
  copula_parameter_start <- 2
  copula_start <- claytonCopula(param = copula_parameter_start, dim = 2)
  theta_start <- list(lambda_start, alpha_start, beta_start)
}

# Function to run simulations and calculate MLEs with two-stage optimization and original method
run_simulation <- function(N, theta_start, copula_start, copula_parameter_start) {
  # Simulate data
  sim_data <- simulate_until_N_copulas_batch(N, theta, copula)
  times <- sim_data[[1]]
  ids <- sim_data[[2]]
  
  # Estimate parameters for the standard Hawkes process
  result_standard <- mutual_exp_mle(times, ids, max(times), theta_start)
  
  # Estimate parameters for the frugal Hawkes process with original method
  result_frugal_original <- frugal_mutual_exp_mle(times, ids, max(times), c(theta_start, copula_parameter_start), copula_start)
  
  # Estimate parameters for the frugal Hawkes process with two-stage optimization
  result_frugal_2stage <- frugal_mutual_exp_mle_2stage (times, ids, max(times), theta_start, copula_parameter_start)
  
  # Compute differences
  difference_standard <- abs(unlist(result_standard$theta_mle) - theta_vector)
  difference_frugal_original <- abs(unlist(result_frugal_original$theta_mle) - c(theta_vector, copula_parameter))
  difference_frugal_2stage <- abs(unlist(result_frugal_2stage$theta_mle) - c(theta_vector, copula_parameter))
  
  return(list(
    standard = difference_standard, 
    frugal_original = difference_frugal_original, 
    frugal_2stage = difference_frugal_2stage, 
    logLik_standard = result_standard$logLik, 
    logLik_frugal_original = result_frugal_original$logLik, 
    logLik_frugal_2stage = result_frugal_2stage$logLik
  ))
}

# # Running simulations for some N 
# N_values <- c(3)
# results_inference_200 <- lapply(N_values, function(N) run_simulation(N, theta_start, copula_start, copula_parameter_start))
# results_inference_200

# Function to run the simulation and compute SD, SE, and CP
run_simulation_repeat <- function(N_values, theta_start, copula_start, copula_parameter_start, 
                                  true_params_standard, true_params_frugal, num_reps = 200) {
  results_list <- list()
  
  for (N in N_values) {
    # Initialize storage for repetitions
    repetitions_standard <- matrix(0, nrow = num_reps, ncol = length(true_params_standard))
    repetitions_frugal_original <- matrix(0, nrow = num_reps, ncol = length(true_params_frugal))
    repetitions_frugal_2stage <- matrix(0, nrow = num_reps, ncol = length(true_params_frugal))
    
    # Repeat the simulation num_reps times
    for (rep in 1:num_reps) {
      sim_result <- run_simulation(N, theta_start, copula_start, copula_parameter_start)
      repetitions_standard[rep, ] <- sim_result$standard
      repetitions_frugal_original[rep, ] <- sim_result$frugal_original
      repetitions_frugal_2stage[rep, ] <- sim_result$frugal_2stage
      print(sim_result$standard)
      print(sim_result$frugal_original)
      print(sim_result$frugal_2stage)
      print(rep)
    }
    
    # Compute SD
    sd_standard <- apply(repetitions_standard, 2, sd)
    sd_frugal_original <- apply(repetitions_frugal_original, 2, sd)
    sd_frugal_2stage <- apply(repetitions_frugal_2stage, 2, sd)
    
    print( sd_standard )
    
    # Compute SE (Standard Error)
    se_standard <- sd_standard / sqrt(num_reps)
    se_frugal_original <- sd_frugal_original / sqrt(num_reps)
    se_frugal_2stage <- sd_frugal_2stage / sqrt(num_reps)
    
    # Compute CP (Coverage Probability)
    ci_standard <- (repetitions_standard >= (true_params_standard - 1.96 * sd_standard)) & (repetitions_standard <= (true_params_standard + 1.96 * sd_standard))
    ci_frugal_original <- (repetitions_frugal_original >= (true_params_frugal - 1.96 * sd_frugal_original)) & (repetitions_frugal_original <= (true_params_frugal + 1.96 * sd_frugal_original))
    ci_frugal_2stage <- (repetitions_frugal_2stage >= (true_params_frugal - 1.96 * sd_frugal_2stage)) & (repetitions_frugal_2stage <= (true_params_frugal + 1.96 * sd_frugal_2stage))
    
    cp_standard <- colMeans(ci_standard)
    cp_frugal_original <- colMeans(ci_frugal_original)
    cp_frugal_2stage <- colMeans(ci_frugal_2stage)
    
    # Store results
    results_list[[as.character(N)]] <- list(
      sd_standard = sd_standard,
      sd_frugal_original = sd_frugal_original,
      sd_frugal_2stage = sd_frugal_2stage,
      se_standard = se_standard,
      se_frugal_original = se_frugal_original,
      se_frugal_2stage = se_frugal_2stage,
      cp_standard = cp_standard,
      cp_frugal_original = cp_frugal_original,
      cp_frugal_2stage = cp_frugal_2stage
    )
  }
  
  return(results_list)
}

# Example usage
true_params_standard <- theta_vector  # Replace with actual true parameter values for standard method
true_params_frugal <- c(theta_vector, copula_parameter)  # Replace with actual true parameter values for frugal method
results_summary <- run_simulation_repeat(N_values = c(300), theta_start = theta_start, 
                                         copula_start = copula_start, copula_parameter_start = copula_parameter_start, 
                                         true_params_standard = true_params_standard, true_params_frugal = true_params_frugal)

results_summary





