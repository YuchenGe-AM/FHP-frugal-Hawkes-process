
source("Simulation_Preparation.R")
source("Inference_Preparation.R")

# Simulation of d-variate multivaraite processes assumes: 
# 1. univariate marginal process = Hawkes process with exponential decay
# 2. dependence characterizes by some copula

# Function to generate the next time t_{i+1} using inverse compensator
generate_next_time <- function(U, t_i, times, ids, theta) {
  root_function <- function(t_next, t_i, u, times, ids, theta, dim) {
    compensator_value <- mutual_exp_hawkes_compensator(t_next, times, ids, theta)[dim] - 
      mutual_exp_hawkes_compensator(t_i, times, ids, theta)[dim] - (-log(u))
    return(compensator_value)
  }
  
  find_root <- function(t_i, u, times, ids, theta, dim) {
    uniroot(root_function, c(t_i, t_i + 500), t_i = t_i, u = u, times = times, ids = ids, theta = theta, dim = dim)$root
  }
  
  # Generate the next time for each dimension based on U
  times_new <- sapply(1:length(U), function(i) {
    find_root(t_i, U[i], times, ids, theta, i)
  })
  
  return(times_new)
}

# Simulate events until time T via copulas
simulate_until_T_copulas <- function(T, theta, copula) {
  # Initialize event times and dimensions
  times <- numeric(0)  # Start with no events
  ids <- integer(0)    # Start with no events
  t_i <- 0 # Starting time
  U_all <- rCopula(100000, copula) # Sample many copulas
  # Simulate events until time T 
  while (TRUE) {
    if (t_i > T) {
      break
    }
    # Generate Clayton copula data
    U <- U_all[sample(1:nrow(U_all), 1), ]
    # print(U)
    
    # Compute the next times for each dimension
    next_times <- generate_next_time(U, t_i, times, ids, theta)
    
    # Determine the minimum time and corresponding index
    t_next <- min(next_times)
    min_index <- which.min(next_times)
    
    # Add the new event to times and ids
    times <- c(times, t_next)
    ids <- c(ids, min_index)
    
    # Update t_i for the next iteration
    t_i <- t_next 
  }
  
  return(list(times, ids))
}

# Simulate N events via copulas
simulate_until_N_copulas <- function(N, theta, copula) {
  times <- numeric(0)
  ids <- integer(0)
  t_i <- 0
  
  U_all <- rCopula(100000, copula)  # Sample many copulas
  
  # Simulate until N events
  while (length(times) < N) {
    # Sample a copula
    U <- U_all[sample(1:nrow(U_all), 1), ]
    
    # Generate the next event times for each dimension
    next_times <- generate_next_time(U, t_i, times, ids, theta)
    
    # Determine the minimum time and corresponding index
    t_next <- min(next_times)
    min_index <- which.min(next_times)
    
    # Add the new event to times and ids
    times <- c(times, t_next)
    ids <- c(ids, min_index)
    
    # Update t_i for the next iteration
    t_i <- t_next
  }
  
  return(list(times, ids))
}

# Simulate N events via choose 1st marginal elements of copulas
simulate_until_N_copulas_marginal1 <- function(N, theta, copula) {
  times <- numeric(0)
  ids <- integer(0)
  t_i <- 0
  
  U_all <- rCopula(100000, copula)  # Sample many copulas
  
  # Simulate until N events
  while (length(times) < N) {
    # Sample a copula
    U <- U_all[sample(1:nrow(U_all), 1), ]
    
    # Generate the next event times for each dimension
    next_times <- generate_next_time(U, t_i, times, ids, theta)
    
    # Determine the minimum time and corresponding index
    t_next <- next_times[1] # instead of min(next_times)
    min_index <- 1 # instead of which.min(next_times)
    
    # Add the new event to times and ids
    times <- c(times, t_next)
    ids <- c(ids, min_index)
    
    # Update t_i for the next iteration
    t_i <- t_next
  }
  
  return(list(times, ids))
}

# Simulate events until time T via batch sampling copulas
simulate_until_T_copulas_batch <- function(T, theta, copula, batch_size = 1000) {
  # Initialize event times and dimensions
  times <- numeric(0)  # Start with no events
  ids <- integer(0)    # Start with no events
  t_i <- 0  # Starting time
  
  while (TRUE) {
    if (t_i > T) {
      break
    }
    
    # Sample a batch of copulas as needed
    U_batch <- rCopula(batch_size, copula)
    
    for (j in 1:nrow(U_batch)) {
      U <- U_batch[j, ]
      # Compute the next times for each dimension
      next_times <- generate_next_time(U, t_i, times, ids, theta)
      
      # Determine the minimum time and corresponding index
      t_next <- min(next_times)
      min_index <- which.min(next_times)
      
      # Add the new event to times and ids
      times <- c(times, t_next)
      ids <- c(ids, min_index)
      
      # Update t_i for the next iteration
      t_i <- t_next
      
      # Break if we've exceeded the horizon T
      if (t_i > T) break
    }
  }
  
  return(list(times, ids))
}

# Simulate N events via batch sampling copulas
simulate_until_N_copulas_batch <- function(N, theta, copula) {
  times <- numeric(0)
  ids <- integer(0)
  t_i <- 0
  
  batch_size <- 1000  # decide batch size, eg. 1000, 500
  while (length(times) < N) {
    U_batch <- rCopula(batch_size, copula)
    for (j in 1:nrow(U_batch)) {
      U <- U_batch[j, ]
      
      next_times <- generate_next_time(U, t_i, times, ids, theta)
      t_next <- min(next_times)
      min_index <- which.min(next_times)
      times <- c(times, t_next)
      ids <- c(ids, min_index)
      t_i <- t_next
      
      if (length(times) >= N) break
    }
  }
  return(list(times, ids))
}

# Simulate events until time T via adaptive waiting time
simulate_until_T_copulas_adaptivewait <- function(T, theta, copula) {
  # Initialize event times and dimensions
  times <- numeric(0)  # Start with no events
  ids <- integer(0)    # Start with no events
  t_i <- 0  # Starting time
  
  U_all <- rCopula(1000000, copula)  # Sample many copulas
  
  # Simulate events until time T 
  while (TRUE) {
    if (t_i > T) {
      break
    }
    
    # Shuffle U_all periodically
    if (length(times) %% 100 == 0) {
      U_all <- U_all[sample(1:nrow(U_all)), ]
    }
    
    # Sample a small batch of copulas and use them
    U_batch <- U_all[sample(1:nrow(U_all), 10), ]
    
    for (j in 1:nrow(U_batch)) {
      U <- U_batch[j, ]
      
      # Compute the next times for each dimension
      next_times <- generate_next_time(U, t_i, times, ids, theta)
      
      # Determine the minimum time and corresponding index
      t_next <- min(next_times)
      min_index <- which.min(next_times)
      
      # Add the new event to times and ids
      times <- c(times, t_next)
      ids <- c(ids, min_index)
      
      # Update t_i for the next iteration
      t_i <- t_next + runif(1)/4
      
      # Break the loop if the next event time exceeds T
      if (t_i > T) {
        break
      }
    }
  }
  
  return(list(times, ids))
}

# Simulate N events until via adaptive waiting time
simulate_until_N_copulas_adaptivewait <- function(N, theta, copula) {
  times <- numeric(0)
  ids <- integer(0)
  t_i <- 0
  
  U_all <- rCopula(1000000, copula)  # Sample many copulas
  
  # Simulate until N events
  while (length(times) < N) {
    # Shuffle U_all periodically
    if (length(times) %% 100 == 0) {
      U_all <- U_all[sample(1:nrow(U_all)), ]
    }
    
    # Sample a small batch of copulas and use them
    U_batch <- U_all[sample(1:nrow(U_all), 10), ]
    
    for (j in 1:nrow(U_batch)) {
      U <- U_batch[j, ]
      
      # Compute the next times for each dimension
      next_times <- generate_next_time(U, t_i, times, ids, theta)
      
      # Determine the minimum time and corresponding index
      t_next <- min(next_times)
      min_index <- which.min(next_times)
      
      # Add the new event to times and ids
      times <- c(times, t_next)
      ids <- c(ids, min_index)
      
      # Update t_i for the next iteration
      t_i <- t_next + runif(1)/4
      
      # Break the loop if we have reached N events
      if (length(times) >= N) break
    }
  }
  
  return(list(times, ids))
}

# Simulation of d-variate multivaraite processes assumes: 
# 1. univariate marginal process = Hawkes process with exponential decay
# 2. NO dependence via simulating each marginal process separately

# Simulate events until time T via marginal simulation with hawkes package
simulate_until_T_marginal <- function(T, theta) {
  # Initialize event times and dimensions
  times <- numeric(0)  # Start with no events
  ids <- integer(0)    # Start with no events
  d <- length(theta[[1]])  # Number of dimensions
  
  # Simulate events for each dimension separately
  for (dim in 1:d) {
    # Simulate the entire process up to time T
    hawkes_process <- simulateHawkes(lambda0 = theta[[1]][dim], alpha = theta[[2]][dim, dim], beta = theta[[3]][dim], horizon = T)
    
    # If there are no events, skip this dimension
    if (length(hawkes_process) == 0) {
      next
    }
    
    # Add the new events to times and ids
    hawkes_process_unlist <- unlist(hawkes_process)
    times <- c(times, hawkes_process_unlist)
    ids <- c(ids, rep(dim, length(hawkes_process_unlist)))
  }
  
  # Sort the times and ids by event times
  order_idx <- order(times)
  times <- times[order_idx]
  ids <- ids[order_idx]
  
  return(list(times, ids))
}

# Simulate N events via marginal simulation with hawkes package
simulate_until_N_marginal <- function(N, theta) {
  copula_inside <- normalCopula(param = 0, dim = 2)
  list_inside <- simulate_until_N_copulas(N, theta, copula_inside)
  return(list_inside)
}

# Simulation of d-variate multivaraite processes assumes: 
# 1. we know the conditional intensity (Ogata's thinning method)

# Simulate events until time T via marginal simulation with hawkes package
simulate_until_T_thinning <- function(T, theta, copula, copula_parameter) {
  # Initialize event times and dimensions
  times <- numeric(0)  # Start with no events
  ids <- integer(0)    # Start with no events
  
  t <- 0  # Starting time
  d <- length(theta[[1]])  # Number of dimensions
  
  while (t < T) {
    # Calculate the intensity for the current time t
    lambda_x <- frugal_mutual_exp_hawkes_intensity(t, times, ids, theta, copula, copula_parameter)
    
    # Sum of intensities (upper bound for the thinning process)
    M <- sum(lambda_x)
    
    # Generate waiting time until the next event
    Dt <- rexp(1, rate = M)
    t <- t + Dt
    
    # If t exceeds T, stop the process
    if (t > T) break
    
    # Generate a uniform random variable to decide acceptance
    u <- runif(1)
    
    if (u <= sum(lambda_x)/M) {
      # Determine which dimension the event occurred in
      cumulative_lambda <- cumsum(lambda_x)
      dim <- which(u <= cumulative_lambda/M)[1]
      
      # Add the event to times and ids
      times <- c(times, t)
      ids <- c(ids, dim)
      
      # Update the intensity for the next iteration
      lambda_x[dim] <- lambda_x[dim] + theta[[2]][dim, dim]
    }
  }
  
  return(list(times, ids))
}

# Simulate N events via marginal simulation with hawkes package
simulate_until_N_thinning <- function(N, theta, copula, copula_parameter) {
  times <- numeric(0)
  ids <- integer(0)
  t <- 0  # Starting time
  d <- length(theta[[1]])  # Number of dimensions
  
  while (length(times) < N) {
    lambda_x <- frugal_mutual_exp_hawkes_intensity(t, times, ids, theta, copula, copula_parameter)
    M <- sum(lambda_x)
    
    # Generate waiting time until the next event
    Dt <- rexp(1, rate = M)
    t <- t + Dt
    
    u <- runif(1)
    
    if (u <= sum(lambda_x) / M) {
      cumulative_lambda <- cumsum(lambda_x)
      dim <- which(u <= cumulative_lambda / M)[1]
      
      times <- c(times, t)
      ids <- c(ids, dim)
      
      lambda_x[dim] <- lambda_x[dim] + theta[[2]][dim, dim]
    }
  }
  
  return(list(times, ids))
}

# # Example setup
# {
#   # Initialize parameters and Clayton copula parameter for 3 univariate processes
#   lambda <- c(1.2, 1)   # Baseline intensities 
#   alpha <- matrix(c(1, 0,   # Excitation effects 
#                     0, 1), nrow = 2, byrow = TRUE)
#   beta <- c(1.5, 1.5)     # Decay rates
#   copula_parameter <- 4
#   copula <- claytonCopula(param = copula_parameter, dim = 2)
#   theta <- list(lambda, alpha, beta)
#   
#   # list_copulas <- simulate_until_T_copulas(T=20, theta, copula)
#   # list_marginal <- simulate_until_T_marginal(T=20, theta)
#   # list_thinning <- simulate_until_T_thinning(T=15, theta, copula, copula_parameter)
# }

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
num <- 230:250


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

