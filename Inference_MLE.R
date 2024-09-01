
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

# Running simulations for differnt N 
N_values <- c(100,200)

results_inference_100_200 <- lapply(N_values, function(N) run_simulation(N, theta_start, copula_start, copula_parameter_start))

N_values <- c(300,400)

results_inference_300_400 <- lapply(N_values, function(N) run_simulation(N, theta_start, copula_start, copula_parameter_start))

N_values <- c(500,600)

results_inference_500_600 <- lapply(N_values, function(N) run_simulation(N, theta_start, copula_start, copula_parameter_start))

# construct a describing table
{
  # Define a function to compute Mean, SD, and SE for a list of numeric vectors
  compute_statistics <- function(values) {
    mean_val <- mean(values)
    sd_val <- sd(values)
    se_val <- sd_val / sqrt(length(values))
    return(c(Mean = mean_val, SD = sd_val, SE = se_val))
  }
  
  # Initialize empty lists to store statistics
  standard_stats <- list()
  frugal_original_stats <- list()
  frugal_2stage_stats <- list()
  
  # Calculate statistics for each N
  for (i in seq_along(N_values)) {
    # Extract results for the current N
    current_results <- results_inference[[i]]
    
    # Calculate statistics for standard, frugal original, and frugal 2-stage
    standard_stats[[i]] <- compute_statistics(current_results$standard)
    frugal_original_stats[[i]] <- compute_statistics(current_results$frugal_original)
    frugal_2stage_stats[[i]] <- compute_statistics(current_results$frugal_2stage)
  }
  
  # Combine the results into data frames
  standard_df <- do.call(rbind, standard_stats)
  frugal_original_df <- do.call(rbind, frugal_original_stats)
  frugal_2stage_df <- do.call(rbind, frugal_2stage_stats)
  
  # Add N values and labels for identification
  standard_df$N <- N_values
  frugal_original_df$N <- N_values
  frugal_2stage_df$N <- N_values
  
  # Combine all into a single data frame for easy comparison
  results_table <- rbind(
    data.frame(Method = "Standard", standard_df),
    data.frame(Method = "Frugal Original", frugal_original_df),
    data.frame(Method = "Frugal 2-Stage", frugal_2stage_df)
  )
  
  # View the results table
  print(results_table)
}

