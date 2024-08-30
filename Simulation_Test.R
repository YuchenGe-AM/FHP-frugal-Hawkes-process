
source("~/Desktop/Code/Simulation_Preparation.R")
source("~/Desktop/Code/Simulation_Methods.R")
source("~/Desktop/Code/Inference_Preparation.R")

# Example usage
lambda <- c(0.2, 0.23)   # Baseline intensities for both dimensions
alpha <- matrix(c(0.32, 0,   # Excitation effects from dim 1 to 1 and 2
                  0, 0.04),  # Excitation effects from dim 2 to 1 and 2
                nrow = 2, byrow = TRUE)
beta <- c(3.9, 4.8)     # Decay rates for both dimensions
copula_parameter <- 2
copula <- claytonCopula(param = copula_parameter, dim = 2)
theta <- list(lambda, alpha, beta)

# Goodness-of-fit test
# Initialize results for KS, p_value methods, and other tests
results <- list()

methods <- c("copula", "copula_batch_resampling", "thinning", "marginal", "adaptive_waiting_time")
for (method in methods) {
  results[[method]] <- list(
    ks_stat1 = numeric(0), ks_stat2 = numeric(0), 
    ks_pval1 = numeric(0), ks_pval2 = numeric(0),
    cvm_stat1 = numeric(0), cvm_stat2 = numeric(0),
    cvm_pval1 = numeric(0), cvm_pval2 = numeric(0),
    ad_stat1 = numeric(0), ad_stat2 = numeric(0),
    ad_pval1 = numeric(0), ad_pval2 = numeric(0)
  )
}

# Function to update results
update_results <- function(method, ks_result1, ks_result2, cvm_result1, cvm_result2, 
                           ad_result1, ad_result2) {
  results[[method]]$ks_stat1 <- c(results[[method]]$ks_stat1, ks_result1$statistic)
  results[[method]]$ks_stat2 <- c(results[[method]]$ks_stat2, ks_result2$statistic)
  results[[method]]$ks_pval1 <- c(results[[method]]$ks_pval1, ks_result1$p.value)
  results[[method]]$ks_pval2 <- c(results[[method]]$ks_pval2, ks_result2$p.value)
  results[[method]]$cvm_stat1 <- c(results[[method]]$cvm_stat1, cvm_result1$statistic)
  results[[method]]$cvm_stat2 <- c(results[[method]]$cvm_stat2, cvm_result2$statistic)
  results[[method]]$cvm_pval1 <- c(results[[method]]$cvm_pval1, cvm_result1$p.value)
  results[[method]]$cvm_pval2 <- c(results[[method]]$cvm_pval2, cvm_result2$p.value)
  results[[method]]$ad_stat1 <- c(results[[method]]$ad_stat1, ad_result1$statistic)
  results[[method]]$ad_stat2 <- c(results[[method]]$ad_stat2, ad_result2$statistic)
  results[[method]]$ad_pval1 <- c(results[[method]]$ad_pval1, ad_result1$p.value)
  results[[method]]$ad_pval2 <- c(results[[method]]$ad_pval2, ad_result2$p.value)
  return(results)
}

# some T for the Thinning method
for (T in 60:100) {
  list_thinning <- simulate_until_T_thinning(T, theta, copula, copula_parameter)
  times <- list_thinning[[1]]
  ids <- list_thinning[[2]]
  rescaled_times <- residual_analysis(times, ids, theta, copula)
  
  # Goodness-of-fit tests
  ks_result1 <- ks_test_rescaled_times(rescaled_times[[1]])
  ks_result2 <- ks_test_rescaled_times(rescaled_times[[2]])
  cvm_result1 <- custom_test_cvm(rescaled_times[[1]])
  cvm_result2 <- custom_test_cvm(rescaled_times[[2]])
  ad_result1 <- custom_test_ad(rescaled_times[[1]])
  ad_result2 <- custom_test_ad(rescaled_times[[2]])
  
  cat("Thinning Method: T =", T, "For process 1 KS Stat =", ks_result1$statistic, "p-value =", ks_result1$p.value, "\n")
  cat("Thinning Method: T =", T, "For process 2 KS Stat =", ks_result2$statistic, "p-value =", ks_result2$p.value, "\n")
  results <- update_results("thinning", ks_result1, ks_result2, cvm_result1, cvm_result2, ad_result1, ad_result2)
}

# some T for Copula method
for (T in 60:100) {
  list_copula <- simulate_until_T_copulas(T, theta, copula)
  times <- list_copula[[1]]
  ids <- list_copula[[2]]
  rescaled_times <- residual_analysis(times, ids, theta, copula)
  
  # Goodness-of-fit tests
  ks_result1 <- ks_test_rescaled_times(rescaled_times[[1]])
  ks_result2 <- ks_test_rescaled_times(rescaled_times[[2]])
  cvm_result1 <- custom_test_cvm(rescaled_times[[1]])
  cvm_result2 <- custom_test_cvm(rescaled_times[[2]])
  ad_result1 <- custom_test_ad(rescaled_times[[1]])
  ad_result2 <- custom_test_ad(rescaled_times[[2]])
  
  cat("Copula Method: T =", T, "For process 1 KS Stat =", ks_result1$statistic, "p-value =", ks_result1$p.value, "\n")
  cat("Copula Method: T =", T, "For process 2 KS Stat =", ks_result2$statistic, "p-value =", ks_result2$p.value, "\n")
  results <- update_results("copula", ks_result1, ks_result2, cvm_result1, cvm_result2, ad_result1, ad_result2)
}

# some T for Copula Batch Resampling method
for (T in 60:100) {
  list_copula_batch <- simulate_until_T_copulas_batch(T, theta, copula)
  times <- list_copula_batch[[1]]
  ids <- list_copula_batch[[2]]
  rescaled_times <- residual_analysis(times, ids, theta, copula)
  
  # Goodness-of-fit tests
  ks_result1 <- ks_test_rescaled_times(rescaled_times[[1]])
  ks_result2 <- ks_test_rescaled_times(rescaled_times[[2]])
  cvm_result1 <- custom_test_cvm(rescaled_times[[1]])
  cvm_result2 <- custom_test_cvm(rescaled_times[[2]])
  ad_result1 <- custom_test_ad(rescaled_times[[1]])
  ad_result2 <- custom_test_ad(rescaled_times[[2]])
  
  cat("Copula Batch Resampling Method: T =", T, "For process 1 KS Stat =", ks_result1$statistic, "p-value =", ks_result1$p.value, "\n")
  cat("Copula Batch Resampling Method: T =", T, "For process 2 KS Stat =", ks_result2$statistic, "p-value =", ks_result2$p.value, "\n")
  results <- update_results("copula_batch_resampling", ks_result1, ks_result2, cvm_result1, cvm_result2, ad_result1, ad_result2)
}

# some T for Adaptive Waiting Time method
for (T in 60:100) {
  list_adaptivewaittime <- simulate_until_T_copulas_adaptivewait(T, theta, copula)
  times <- list_adaptivewaittime[[1]]
  ids <- list_adaptivewaittime[[2]]
  rescaled_times <- residual_analysis(times, ids, theta, copula)
  
  # Goodness-of-fit tests
  ks_result1 <- ks_test_rescaled_times(rescaled_times[[1]])
  ks_result2 <- ks_test_rescaled_times(rescaled_times[[2]])
  cvm_result1 <- custom_test_cvm(rescaled_times[[1]])
  cvm_result2 <- custom_test_cvm(rescaled_times[[2]])
  ad_result1 <- custom_test_ad(rescaled_times[[1]])
  ad_result2 <- custom_test_ad(rescaled_times[[2]])
  
  cat("Adaptive Waiting Time Method: T =", T, "For process 1 KS Stat =", ks_result1$statistic, "p-value =", ks_result1$p.value, "\n")
  cat("Adaptive Waiting Time Method: T =", T, "For process 2 KS Stat =", ks_result2$statistic, "p-value =", ks_result2$p.value, "\n")
  results <- update_results("adaptive_waiting_time", ks_result1, ks_result2, cvm_result1, cvm_result2, ad_result1, ad_result2)
}

# some T for Marginal method
for (T in 60:100) {
  list_marginal <- simulate_until_T_marginal(T, theta)
  times <- list_marginal[[1]]
  ids <- list_marginal[[2]]
  rescaled_times <- residual_analysis(times, ids, theta, copula)
  
  # Goodness-of-fit tests
  ks_result1 <- ks_test_rescaled_times(rescaled_times[[1]])
  ks_result2 <- ks_test_rescaled_times(rescaled_times[[2]])
  cvm_result1 <- custom_test_cvm(rescaled_times[[1]])
  cvm_result2 <- custom_test_cvm(rescaled_times[[2]])
  ad_result1 <- custom_test_ad(rescaled_times[[1]])
  ad_result2 <- custom_test_ad(rescaled_times[[2]])
  
  cat("Marginal Method: T =", T, "For process 1 KS Stat =", ks_result1$statistic, "p-value =", ks_result1$p.value, "\n")
  cat("Marginal Method: T =", T, "For process 2 KS Stat =", ks_result2$statistic, "p-value =", ks_result2$p.value, "\n")
  results <- update_results("marginal", ks_result1, ks_result2, cvm_result1, cvm_result2, ad_result1, ad_result2)
}

# Plot results for each test
# plot_results_dynamic(results, "ks_stat", "KS Statistic")
plot_results_dynamic(results, "ks_pval", "KS p-value")
# plot_results_dynamic(results, "cvm_stat", "CVM Statistic")
# plot_results_dynamic(results, "cvm_pval", "CVM p-value")
# plot_results_dynamic(results, "ad_stat", "AD Statistic")
# plot_results_dynamic(results, "ad_pval", "AD p-value")