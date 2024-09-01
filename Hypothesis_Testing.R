
# Load necessary simulation and inference preparation scripts
source("Simulation_Preparation.R")
source("Simulation_Methods.R")
source("Inference_Preparation.R")


# Example usage
lambda <- c(0.2, 0.23)   # Baseline intensities for both dimensions
alpha <- matrix(c(0.32, 0,   # Excitation effects from dim 1 to 1 and 2
                  0, 0.54),  # Excitation effects from dim 2 to 1 and 2
                nrow = 2, byrow = TRUE)
beta <- c(4.9, 4.8)     # Decay rates for both dimensions
copula_parameter <- 1.4
copula <- claytonCopula(param = copula_parameter, dim = 2)
theta <- list(lambda, alpha, beta)
timeline <- 20:150
num <- 207:250

for (N in num) {
  list_copula_batch <- simulate_until_N_copulas_batch(N, theta, copula)
  times <- list_copula_batch[[1]]
  ids <- list_copula_batch[[2]]
  # print(ids)
  rescaled_times <- residual_analysis(times, ids, theta, copula)
  # Goodness-of-fit tests
  ks_result1 <- ks_test_rescaled_times(rescaled_times[[1]])
  ks_result2 <- ks_test_rescaled_times(rescaled_times[[2]])
  
  cat("Copula Batch Resampling Method: N =", N, "For process 1 KS Stat =", ks_result1$statistic, "p-value =", ks_result1$p.value, "\n")
  cat("Copula Batch Resampling Method: N =", N, "For process 2 KS Stat =", ks_result2$statistic, "p-value =", ks_result2$p.value, "\n")
}

















# Define parameters
n_simulations <- 100
N <- 100  # Number of events
alpha <- 0.05  # Significance level

# Define copulas
independent_copula <- indepCopula(dim = 2)
dependent_copula <- claytonCopula(param = 2, dim = 2)

# Store results
wrong_tests_indep <- 0
wrong_tests_dep <- 0

# Function to perform hypothesis test for independence
test_independence <- function(times, ids) {
  estimated_copula <- fitCopula(independent_copula, cbind(times, ids), method = "mpl")$copula
  gof_result <- gofCopula(estimated_copula, pobs(cbind(times, ids)), test = "CvM")
  return(gof_result$p.value)
}

# Simulate and test for independent copula
for (i in 1:n_simulations) {
  sim_data <- simulate_until_N_copulas(N, theta, independent_copula)
  p_value <- test_independence(sim_data[[1]], sim_data[[2]])
  if (p_value < alpha) wrong_tests_indep <- wrong_tests_indep + 1
}

# Simulate and test for dependent copula
for (i in 1:n_simulations) {
  sim_data <- simulate_until_N_copulas(N, theta, dependent_copula)
  p_value <- test_independence(sim_data[[1]], sim_data[[2]])
  if (p_value >= alpha) wrong_tests_dep <- wrong_tests_dep + 1
}

# Output the number of incorrect tests
cat("Number of incorrect tests for independent copula:", wrong_tests_indep, "\n")
cat("Number of incorrect tests for dependent copula:", wrong_tests_dep, "\n")