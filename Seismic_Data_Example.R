
source("Simulation_Preparation.R")
source("Simulation_Methods.R")
source("Inference_Preparation.R")

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

# Load the data 
{
  # Define file paths for the Excel files
  file1 <- "Data/Ogasawara Transformed.xlsx"
  file2 <- "Data/Noto Peninsula Transformed.xlsx" # Asuncion, Mabula, Tonga, Parque 
  
  data1 <- read_excel(file1)$Value
  data2 <- read_excel(file2)$Value
  
  times <- c(data1, data2)
  ids <- c(rep(1, length(data1)), rep(2, length(data2)))
  
  # Create a data frame combining times and IDs
  combined_data <- data.frame(times = times, ids = ids)
  
  # Sort the data by times
  sorted_data <- combined_data[order(combined_data$times), ]
  
  # Extract the sorted times and ids
  new_times <- sorted_data$times
  new_ids <- sorted_data$ids
}

# Split the real data into subprocess 1 and subprocess 2 based on IDs
real_subprocess1 <- new_times[new_ids == 1]
real_subprocess2 <- new_times[new_ids == 2]

real_subprocess2
# Combine real and simulated data into a long format for plotting
real_data <- data.frame(times = c(real_subprocess1, real_subprocess2),
                        process = rep(c("Subprocess 1", "Subprocess 2"), 
                                      times = c(length(real_subprocess1), length(real_subprocess2))),
                        type = "Real Data")

# simulation of the real data
result_frugal_original <- frugal_mutual_exp_mle(new_times, new_ids, max(new_times), 
                                                c(theta_start, copula_parameter_start), copula_start)

# result_frugal_original
{
  theta <- list(result_frugal_original$theta_mle[[1]], 
                result_frugal_original$theta_mle[[2]], 
                result_frugal_original$theta_mle[[3]])
  copula_parameter <- result_frugal_original$theta_mle[[4]]
  copula <- claytonCopula(param = copula_parameter, dim = 2)
  theta
  copula_parameter
  
  # Simulate data
  simulated_list <- simulate_until_T_copulas_batch(T = 605, theta, copula)
  simulated_times <- simulated_list[[1]]
  simulated_ids <- simulated_list[[2]]
  simulated_times
  
  # Extract simulated times for subprocess 1 and 2
  simulated_subprocess1 <- simulated_times[simulated_ids == 1]
  simulated_subprocess2 <- simulated_times[simulated_ids == 2]
  
  simulated_data <- data.frame(times = c(simulated_subprocess1, simulated_subprocess2),
                               process = rep(c("Subprocess 1", "Subprocess 2"), 
                                             times = c(length(simulated_subprocess1), length(simulated_subprocess2))),
                               type = "Simulated Data")
}

# Plot ECDF and simulation data 

combined_data <- rbind(real_data, simulated_data)

# Plot the ECDF for each subprocess (1 and 2) comparing real vs simulated data
{
  ecdf_plot <- ggplot(combined_data, aes(x = times, color = type)) +
    stat_ecdf(geom = "step", size = 1) +
    facet_wrap(~ process, scales = "free") +
    labs(
      title = "Empirical Distribution Function (ECDF)",
      x = "Transformed Times",
      y = "ECDF",
      color = "Data Type"
    ) +
    theme_minimal()
  
  # Display the plot
  print(ecdf_plot)
}


# Ensure the number of quantiles matches the smaller dataset
min_length1 <- min(length(real_subprocess1), length(simulated_subprocess1))
min_length2 <- min(length(real_subprocess2), length(simulated_subprocess2))

# QQ plot
{
  # QQ plot for subprocess 1
  qqplot_data1 <- data.frame(
    real = quantile(real_subprocess1, probs = seq(0, 1, length.out = min_length1)),
    simulated = quantile(simulated_subprocess1, probs = seq(0, 1, length.out = min_length1))
  )
  
  # QQ plot for subprocess 2
  qqplot_data2 <- data.frame(
    real = quantile(real_subprocess2, probs = seq(0, 1, length.out = min_length2)),
    simulated = quantile(simulated_subprocess2, probs = seq(0, 1, length.out = min_length2))
  )
  
  # Plot QQ plot for subprocess 1
  qqplot1 <- ggplot(qqplot_data1, aes(x = real, y = simulated)) +
    geom_point(color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = "Q-Q Plot for Subprocess 1",
      x = "Quantiles of Real Data",
      y = "Quantiles of Simulated Data"
    ) +
    theme_minimal()
  
  # Plot QQ plot for subprocess 2
  qqplot2 <- ggplot(qqplot_data2, aes(x = real, y = simulated)) +
    geom_point(color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = "Q-Q Plot for Subprocess 2",
      x = "Quantiles of Real Data",
      y = "Quantiles of Simulated Data"
    ) +
    theme_minimal()
}

# Display the QQ plots
print(qqplot1)
print(qqplot2)
copula_parameter


