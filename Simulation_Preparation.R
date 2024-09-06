
# List of required libraries
libraries <- c("copula", "Matrix", "stats", "ggplot2", "dplyr", "gridExtra", "cubature", "numDeriv", "hawkes",
               "goftest", "cramer", "gridExtra", "mvtnorm", "boot", "readxl", "dplyr", "readxl")

# Loop to check and install missing libraries
for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib, dependencies = TRUE)
  }
  library(lib, character.only = TRUE)
}

# compute intensity for multivariate Hawkes processes with exp decay
mutual_exp_hawkes_intensity <- function(t, times, ids, theta) {
  lambda <- theta[[1]]
  alpha <- theta[[2]]
  beta <- theta[[3]]
  
  # restrict to time until t
  times <- times[times <= t]
  ids <- ids[times <= t]

  # Iterate over event times and corresponding dimensions
  for (i in seq_along(times)) {
    t_i <- times[i]
    d_i <- ids[i]
    lambda <- lambda + alpha[d_i,] * exp(-beta * (t - t_i)) # Update intensity
  }
  
  return(lambda)
}

# compute compensator function for multivariate Hawkes processes with exp decay
mutual_exp_hawkes_compensator <- function(t, times, ids, theta) {
  lambda <- theta[[1]]
  alpha <- theta[[2]]
  beta <- theta[[3]]
  Lambda <- lambda * t # Initialize compensator
  
  # restrict to time until t
  times <- times[times <= t]
  ids <- ids[times <= t]
  
  # Iterate over event times and corresponding dimensions
  for (i in seq_along(times)) {
    t_i <- times[i]
    d_i <- ids[i]
    Lambda <- Lambda + (alpha[d_i,] / beta) * (1 - exp(-beta * (t - t_i))) # Update compensator
  }
  # print(Lambda)
  return(Lambda)
}

# compute log-likelihood function for multivariate Hawkes processes with exp decay
mutual_exp_log_likelihood <- function(times, ids, T, theta) {
  lambda <- theta[[1]]
  alpha <- theta[[2]]
  beta <- theta[[3]]

  ll <- 0 # initialization of loglikelihood
  lambda_x <- lambda
  
  t_prev <- 0
  
  # Calculate the log-likelihood
  for (i in seq_along(times)) {
    t_i <- times[i]
    d_i <- ids[i]
    
    # Update intensity
    lambda_x <- lambda + (lambda_x - lambda) * exp(-beta * (t_i - t_prev))
    
    # Update log-likelihood
    ll <- ll + log(lambda_x[d_i])
    
    # Update lambda_x for the next iteration
    lambda_x <- lambda_x + alpha[d_i,]
    t_prev <- t_i
  }
  
  # Compensator term
  compensator <- mutual_exp_hawkes_compensator(T, times, ids, theta)
  ll <- ll - sum(compensator)
  
  return(ll)
}

# Transformation T in (4.3)
transformation <- function(t, times, ids, theta, copula, copula_parameter) {
  # restrict to time until t
  times <- times[times <= t]
  ids <- ids[times <= t]
  
  if (all(times > t)) {
    t_i <- 0 
  } else {
    t_i = max( times[times <= t] )
  }
  
  u <- mutual_exp_hawkes_compensator(t, times, ids, theta) - mutual_exp_hawkes_compensator(t_i, times, ids, theta)
  Lambda_diff = exp(-u)
  copula_value <- pCopula(Lambda_diff, copula) # Compute D(.)
  
  # Compute D^(i)(.) for Clayton copula
  v <- Lambda_diff
  partial_derivative_1 <- v[1]^{-1-copula_parameter} * ( v[1]^{-copula_parameter} + v[2]^{-copula_parameter} - 1 )^{-1-1/copula_parameter}
  partial_derivative_2 <- v[2]^{-1-copula_parameter} * ( v[1]^{-copula_parameter} + v[2]^{-copula_parameter} - 1 )^{-1-1/copula_parameter}
  partial_derivative <- c(partial_derivative_1, partial_derivative_2)
  
  T_t <- Lambda_diff * partial_derivative /copula_value
  # print(Lambda_diff)
  # print(copula_value)
  # print(partial_derivative)
  return(T_t)
}

# transformation(t=200, times, ids, theta, copula, copula_parameter)
  
# Compute modified intensity for frugal Hawkes processes with univaraite marginals
frugal_mutual_exp_hawkes_intensity <- function(t, times, ids, theta, copula, copula_parameter) {
  intensity <- mutual_exp_hawkes_intensity(t, times, ids, theta)
  T_t <- transformation(t, times, ids, theta, copula, copula_parameter)
  
  intensity_frugal <- intensity * T_t # take component-wise product
  return(intensity_frugal)
}

# Compute modified compensator for frugal Hawkes processes with univaraite marginals
frugal_mutual_exp_hawkes_compensator <- function(t, times, ids, theta, copula, copula_parameter) {
  # Define the function to integrate
  integrand1 <- function(s) {
    relevant_times <- times[times <= s]
    relevant_ids <- ids[times <= s]
    intensity <- frugal_mutual_exp_hawkes_intensity(s, relevant_times, relevant_ids, theta, copula, copula_parameter)[1]
    return(intensity)
  }
  integrand2 <- function(s) {
    relevant_times <- times[times <= s]
    relevant_ids <- ids[times <= s]
    intensity <- frugal_mutual_exp_hawkes_intensity(s, relevant_times, relevant_ids, theta, copula, copula_parameter)[2]
    return(intensity)
  }
  
  # Perform the integration from 0 to t
  result1 <- hcubature(integrand1, lowerLimit = 0, upperLimit = t, tol = 1e-2)$integral
  result2 <- hcubature(integrand2, lowerLimit = 0, upperLimit = t, tol = 1e-2)$integral
  
  # Return the value of the integral
  return(c(result1,result2))
}

# Define the log-likelihood function for frual Hawkes processes with univaraite marginals
frugal_mutual_log_likelihood <- function(times, ids, T, theta, copula, copula_parameter) {
  ll <- 0
  for (i in seq_along(times)) {
    t_i <- times[i]
    d_i <- ids[i]
    
    # Get the history of arrivals before time t_i
    H_i_times <- times[times < t_i]
    H_i_ids <- ids[times < t_i]
    
    # Compute lambda^x for time t_i with history H_i in (3.7)
    if (length(H_i_times) > 0) {
      t_i_1 = max( H_i_times )
    } else { 
      t_i_1 <- 0 
    }
    # print(mutual_exp_hawkes_compensator(t_i, H_i_times, H_i_ids, theta))
    u <- mutual_exp_hawkes_compensator(t_i, H_i_times, H_i_ids, theta) - mutual_exp_hawkes_compensator(t_i_1, H_i_times, H_i_ids, theta)
    v <- exp(-u)
    
    # Return a small positive value if either v[1] or v[2] is missing or too small
    if (is.na(v[1]) || v[1] <= 1e-08 ) {
      v[1] = 1e-8  
    }
    if (is.na(v[2]) || v[2] <= 1e-08) {
      v[2] = 1e-8  
    }
    
    partial_derivative_1 <- v[1]^{-1-copula_parameter} * ( v[1]^{-copula_parameter} + v[2]^{-copula_parameter} - 1 )^{-1-1/copula_parameter}
    partial_derivative_2 <- v[2]^{-1-copula_parameter} * ( v[1]^{-copula_parameter} + v[2]^{-copula_parameter} - 1 )^{-1-1/copula_parameter}
    partial_derivative <- c(partial_derivative_1, partial_derivative_2)
    # print(partial_derivative)
    lambda_x_i <- mutual_exp_hawkes_intensity(t_i, H_i_times, H_i_ids, theta)[d_i] * partial_derivative[d_i] * v[d_i]
    
    # Update the log-likelihood
    ll <- ll + log(lambda_x_i)
  }

  # add the D() term in (3.7)
  u <- mutual_exp_hawkes_compensator(T, times, ids, theta) - mutual_exp_hawkes_compensator(max(times), times, ids, theta)
  Lambda_diff = exp(-u)
  copula_value <- pCopula(Lambda_diff, copula)
  ll <- ll + copula_value
  
  return(ll)
}

# Plot marginal processes based on times, ids, and mutual_exp_hawkes_intensity
Plot_marginal <- function(times, ids) {
  # Define a time vector for plotting (from the first to the last event)
  t_vec <- seq(min(times), max(times) + 0.5, by = 0.0005)
  
  # Calculate the intensity over time for each dimension
  intensity_matrix <- t(sapply(t_vec, function(t) {
    mutual_exp_hawkes_intensity(t, times[times <= t], ids[times <= t], theta)
  }))
  
  # Create a data frame for plotting
  df_intensity <- data.frame(
    Time = rep(t_vec, times = length(lambda)),
    Intensity = as.vector(intensity_matrix),
    Dimension = factor(rep(1:length(lambda), each = length(t_vec)))
  )
  
  # Plot Intensity for Each Subprocess (Dimension) with interconnected points
  ggplot(df_intensity, aes(x = Time, y = Intensity, color = Dimension)) +
    geom_point(size = 0.1) +
    labs(title = "Intensity for Each Marginal Process",
         x = "Time",
         y = "Intensity") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("blue", "red"))
}

# Residual Analysis
# Function to perform residual analysis
residual_analysis <- function(times, ids, theta, copula) {
  residuals1 <- numeric(0)
  residuals2 <- numeric(0)
  for (i in seq_along(times)) {
    t_i <- times[i]
    if (ids[i] == 1) {
      residuals1 <- c(residuals1, frugal_mutual_exp_hawkes_compensator(t_i, times[1:i], ids[1:i], theta, copula, copula_parameter)[1])
    } else if (ids[i] == 2) {
      residuals2 <- c(residuals2, frugal_mutual_exp_hawkes_compensator(t_i, times[1:i], ids[1:i], theta, copula, copula_parameter)[2])
    }
  }
  
  # Apply the time-rescaling theorem
  rescaled_times1 <- diff(c(0, residuals1))
  rescaled_times2 <- diff(c(0, residuals2))
  
  # If the model is correct, these rescaled times should be exponential(1)
  return(list(rescaled_times1,rescaled_times2)) # rescaled_times
}

# Function to perform KS, BMA test on rescaled times
ks_test_rescaled_times <- function(rescaled_times) {
  # Perform the KS test against the exponential distribution with constant rate 
  ks_test <- ks.test(rescaled_times, "pexp", rate = 1)
  return(ks_test)
}

# Use tests from goftest package: cramer.test and ad.test
custom_test_cvm <- function(rescaled_times) {
  goftest::cvm.test(rescaled_times, "pexp")
}

custom_test_ad <- function(rescaled_times) {
  goftest::ad.test(rescaled_times, "pexp")
}

# Function to perform marginal residual analysis
residual_analysis_marginal <- function(times, ids, theta, copula) {
  residuals1 <- numeric(0)
  residuals2 <- numeric(0)
  
  for (i in seq_along(times)) {
    t_i <- times[i]
    if (ids[i] == 1) {
      residuals1 <- c(residuals1, mutual_exp_hawkes_compensator(t_i, times[ids == 1], ids[ids == 1], theta)[1])
    } else if (ids[i] == 2) {
      residuals2 <- c(residuals2, mutual_exp_hawkes_compensator(t_i, times[ids == 2], ids[ids == 2], theta)[2])
    }
  }
  
  # Apply the time-rescaling theorem
  rescaled_times1 <- diff(c(0, residuals1))
  rescaled_times2 <- diff(c(0, residuals2))
  
  # If the model is correct, these rescaled times should be exponential(1)
  return(list(rescaled_times1, rescaled_times2))
}

# Function to create plots dynamically with the correct N values and adaptive handling for many N values
# Function to create plots dynamically with the correct N values and adaptive handling for many N values
plot_results_dynamic_N <- function(results_for_method, metric, ylabel, N_values, title) {
  
  if (paste0(metric, "1") %in% names(results_for_method) && 
      paste0(metric, "2") %in% names(results_for_method)) {
    
    if (length(results_for_method[[paste0(metric, "1")]]) > 0 && 
        length(results_for_method[[paste0(metric, "2")]]) > 0) {
      
      df <- data.frame(
        N = N_values[1:length(results_for_method[[paste0(metric, "1")]])],
        Process1 = -log10(results_for_method[[paste0(metric, "1")]]),
        Process2 = -log10(results_for_method[[paste0(metric, "2")]]),
        ymin_Process1 = -log10(results_for_method[[paste0(metric, "1")]]) - log10(sd(results_for_method[[paste0(metric, "1")]])),
        ymax_Process1 = -log10(results_for_method[[paste0(metric, "1")]]) + log10(sd(results_for_method[[paste0(metric, "1")]])),
        ymin_Process2 = -log10(results_for_method[[paste0(metric, "2")]]) - log10(sd(results_for_method[[paste0(metric, "2")]])),
        ymax_Process2 = -log10(results_for_method[[paste0(metric, "2")]]) + log10(sd(results_for_method[[paste0(metric, "2")]]))
      )
      
      # Define the critical p-value levels
      p_values <- c(0.05, 0.01, 0.001)
      log_p_values <- -log10(p_values)
      
      # Create the first plot for Process 1
      p1 <- ggplot(df, aes(x = N, y = Process1)) +
        geom_line(size = 0.7, color = "black") + 
        geom_ribbon(aes(ymin = ymin_Process1, ymax = ymax_Process1), fill = "grey70", alpha = 0.3) +
        labs(title = paste(title, "Process 1"), y = ylabel, x = "N (number of events)") +
        theme_minimal() +
        scale_y_continuous(sec.axis = sec_axis(~ ., breaks = log_p_values, labels = paste0("p = ", p_values))) +
        geom_hline(yintercept = log_p_values, linetype = "dashed", color = "red") +
        theme(legend.position = "none")
      
      # Create the second plot for Process 2
      p2 <- ggplot(df, aes(x = N, y = Process2)) +
        geom_line(size = 0.7, color = "black") + 
        geom_ribbon(aes(ymin = ymin_Process2, ymax = ymax_Process2), fill = "grey70", alpha = 0.3) +
        labs(title = paste(title, "Process 2"), y = ylabel, x = "N (number of events)") +
        theme_minimal() +
        scale_y_continuous(sec.axis = sec_axis(~ ., breaks = log_p_values, labels = paste0("p = ", p_values))) +
        geom_hline(yintercept = log_p_values, linetype = "dashed", color = "red") +
        theme(legend.position = "none")
      
      # Arrange the combined plots for Process 1 and Process 2
      combined_plot <- grid.arrange(p1) #, p2, ncol = 2)
      
      return(combined_plot)
      
    } else {
      stop("No data available for plotting.")
    }
  } else {
    stop("The specified metric was not found in the results.")
  }
}

# hist(rescaled_times, breaks = 20, probability = TRUE, main = "Rescaled Times Histogram")
# curve(dexp(x, rate = 1), add = TRUE, col = "red")

