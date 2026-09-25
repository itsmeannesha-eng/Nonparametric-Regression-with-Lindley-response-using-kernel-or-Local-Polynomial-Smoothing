# SIMULATION STUDY
# Scenario 1: Constant theta(z) = 2
# Nadaraya-Watson Estimation with Gaussian Kernel
# Lindley Distribution with Covariate-Dependent theta(z)

rm(list = ls())

set.seed(12345)
# 1. SIMULATION SETTINGS

# Number of Monte Carlo replications
B <- 1000

# Sample sizes
sample_sizes <- c(50, 100, 200, 500)

# Required evaluation points
evaluation_grid <- c(0.10, 0.30, 0.50, 0.70, 0.90)

# Candidate bandwidths
bandwidth_grid <- c(
  0.05, 0.075, 0.10, 0.125,
  0.15, 0.175, 0.20, 0.25, 0.30
)

# 2. TRUE THETA FUNCTION

# Scenario 1: Constant theta(z) = 2

theta_fun <- function(z) {
  2
}


# 3. LINDELY RANDOM NUMBER GENERATION

# Lindley density:
#
# f(x; theta) =
# theta^2/(theta+1) * (1+x) * exp(-theta*x)
#
# Mixture representation:
#
# Exponential(theta), probability theta/(theta+1)
# Gamma(2, theta), probability 1/(theta+1)


rlindley <- function(n, theta) {
  
  # If theta is a single value,
  # repeat it n times
  if (length(theta) == 1) {
    theta <- rep(theta, n)
  }
  
  # Check theta length
  if (length(theta) != n) {
    stop("Length of theta must be 1 or equal to n.")
  }
  
  # Check positivity
  if (any(theta <= 0)) {
    stop("Theta must be positive.")
  }
  
  # Uniform random variables
  u <- runif(n)
  
  # Storage for X
  x <- numeric(n)
  
  # Exponential component
 
  ind_exp <- u <= theta / (theta + 1)
  
  if (sum(ind_exp) > 0) {
    
    x[ind_exp] <-
      rexp(
        sum(ind_exp),
        rate = theta[ind_exp]
      )
  }
  
  # Gamma component
  
  ind_gamma <- !ind_exp
  
  if (sum(ind_gamma) > 0) {
    
    x[ind_gamma] <-
      rgamma(
        sum(ind_gamma),
        shape = 2,
        rate = theta[ind_gamma]
      )
  }
  
  return(x)
}

# 4. GAUSSIAN KERNEL

gaussian_kernel <- function(u) {
  
  K <- (1 / sqrt(2 * pi)) *
    exp(-0.5 * u^2)
  
  return(K)
}

# 5. NADARAYA-WATSON ESTIMATOR

nw_estimator <- function(
    x,
    z,
    z0,
    h) {
  
  # Gaussian kernel argument
  u <- (z0 - z) / h
  
  # Gaussian weights
  weights <- gaussian_kernel(u)
  
  # Denominator
  denominator <- sum(weights)
  
  # Numerical safeguard
  if (denominator <= 0) {
    return(mean(x))
  }
  
  # Nadaraya-Watson estimate
  estimate <-
    sum(weights * x) /
    denominator
  
  return(estimate)
}

# 6. ESTIMATE m(z) ON THE REQUIRED GRID
estimate_m_grid <- function(
    x,
    z,
    grid,
    h) {
  
  m_hat <- numeric(length(grid))
  
  for (j in seq_along(grid)) {
    
    m_hat[j] <-
      nw_estimator(
        x = x,
        z = z,
        z0 = grid[j],
        h = h
      )
  }
  
  return(m_hat)
}

# 7. RECOVER THETA(z) FROM m(z)

# Lindley conditional mean:
#
# m(z) =
# [theta(z) + 2] /
# [theta(z) {theta(z) + 1}]
#
# Therefore:
#
# m theta^2 + (m - 1) theta - 2 = 0
#
# Positive root:
#
# theta =
# [1 - m + sqrt((m-1)^2 + 8m)] / (2m)


theta_from_m <- function(m) {
  
  # Avoid zero or negative values
  m <- pmax(m, 1e-10)
  
  theta_hat <-
    (
      1 - m +
        sqrt(
          (m - 1)^2 +
            8 * m
        )
    ) /
    (2 * m)
  
  # Ensure positive theta
  theta_hat <-
    pmax(theta_hat, 1e-10)
  
  return(theta_hat)
}

# 8. ESTIMATE THETA(z)
estimate_theta_grid <- function(
    x,
    z,
    grid,
    h) {
  
  # Estimate conditional mean
  m_hat <-
    estimate_m_grid(
      x = x,
      z = z,
      grid = grid,
      h = h
    )
  
  # Convert estimated mean to theta
  theta_hat <-
    theta_from_m(m_hat)
  
  return(theta_hat)
}

# 9. LEAVE-ONE-OUT CROSS-VALIDATION
cv_nw <- function(
    x,
    z,
    h) {
  
  n <- length(x)
  
  squared_errors <-
    numeric(n)
  
  for (i in seq_len(n)) {
    
    # Training data excluding observation i
    x_train <- x[-i]
    z_train <- z[-i]
    
    # Predict X_i at Z_i
    prediction <-
      nw_estimator(
        x = x_train,
        z = z_train,
        z0 = z[i],
        h = h
      )
    
    # Squared prediction error
    squared_errors[i] <-
      (x[i] - prediction)^2
  }
  
  return(sum(squared_errors))
}

# 10. OPTIMAL BANDWIDTH
select_bandwidth <- function(
    x,
    z,
    bandwidths) {
  
  cv_values <-
    numeric(length(bandwidths))
  
  for (j in seq_along(bandwidths)) {
    
    cv_values[j] <-
      cv_nw(
        x = x,
        z = z,
        h = bandwidths[j]
      )
  }
  
  # Select bandwidth giving minimum CV
  best_index <-
    which.min(cv_values)
  
  best_h <-
    bandwidths[best_index]
  
  return(best_h)
}

# 11. TRAPEZOIDAL INTEGRATION

trapz <- function(x, y) {
  
  if (length(x) != length(y)) {
    stop("x and y must have same length.")
  }
  
  if (length(x) < 2) {
    return(NA)
  }
  
  value <-
    sum(
      diff(x) *
        (
          head(y, -1) +
            tail(y, -1)
        ) / 2
    )
  
  return(value)
}


# 12. POINTWISE PERFORMANCE

calculate_pointwise_performance <- function(
    simulation_data) {
 
  # Mean estimated theta
  
  mean_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = mean
    )
  
  names(mean_estimate)[2] <-
    "mean_theta_hat"
  
  # Variance
  variance_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = var
    )
  
  names(variance_estimate)[2] <-
    "variance"
  
  # True theta
  true_values <-
    unique(
      simulation_data[
        ,
        c("z", "theta_true")
      ]
    )
  
  # Merge results
  result <-
    merge(
      true_values,
      mean_estimate,
      by = "z"
    )
  
  result <-
    merge(
      result,
      variance_estimate,
      by = "z"
    )
  # Bias
  result$bias <-
    result$mean_theta_hat -
    result$theta_true
  # MSE
  result$mse <-
    result$bias^2 +
    result$variance
  
  # RMSE
  result$rmse <-
    sqrt(result$mse)
  
  
  return(result)
}

# 13. INTEGRATED PERFORMANCE

calculate_integrated_performance <- function(
    pointwise_result) {
  
  # Integrated squared error
  ISE <-
    trapz(
      pointwise_result$z,
      pointwise_result$mse
    )
  
  # Integrated MSE
  MISE <- ISE
  
  # Root MISE
  RMISE <- sqrt(MISE)
  
  # Mean Bias over evaluation points
  Mean_Bias <-
    mean(pointwise_result$bias)
  
  result <-
    data.frame(
      ISE = ISE,
      MISE = MISE,
      RMISE = RMISE,
      Mean_Bias = Mean_Bias
    )
  
  return(result)
}

# 14. CREATE RESULT DIRECTORY
output_folder <-
  "Scenario_1_Gaussian_NW_Simulation"

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# 15. STORAGE OBJECTS
pointwise_results <- list()

integrated_results <- list()

bandwidth_results <- list()

counter <- 1

# 16. MAIN SIMULATION LOOP

cat("\n")
cat("\n")
cat("SCENARIO 1: CONSTANT THETA(z) = 2\n")
cat("NADARAYA-WATSON WITH GAUSSIAN KERNEL\n")
cat("\n")


for (n in sample_sizes) {
  
  cat("\n")
  cat("\n")
  cat("Sample size:", n, "\n")
  cat("\n")
  
  
  # Storage for Monte Carlo replications
  simulation_results <-
    vector("list", B)
  
  selected_bandwidths <-
    numeric(B)
  
  # Monte Carlo replications
  for (b in seq_len(B)) {
    
    # Generate covariate Z
    
    Z <-
      runif(
        n,
        min = 0,
        max = 1
      )
    
    # Generate true theta(Z)
    # Scenario 1: theta(z) = 2
    theta_data <-
      theta_fun(Z)
    
    # Generate Lindley response X
    X <-
      rlindley(
        n = n,
        theta = theta_data
      )
    
    # Select optimal bandwidth
    h_opt <-
      select_bandwidth(
        x = X,
        z = Z,
        bandwidths = bandwidth_grid
      )
    
    
    selected_bandwidths[b] <-
      h_opt
    
    # Estimate theta on required evaluation points
    theta_hat <-
      estimate_theta_grid(
        x = X,
        z = Z,
        grid = evaluation_grid,
        h = h_opt
      )
    
    # True theta on evaluation grid
    theta_grid_true <-
      theta_fun(
        evaluation_grid
      )
    
    # Store replication
    
    simulation_results[[b]] <-
      data.frame(
        replication = b,
        z = evaluation_grid,
        theta_true = theta_grid_true,
        theta_hat = theta_hat,
        bandwidth = h_opt
      )
    
    # Progress
    
    if (b %% 100 == 0) {
      
      cat(
        "Replication:",
        b,
        "of",
        B,
        "\n"
      )
    }
  }
  
 
  # 17. COMBINE MONTE CARLO REPLICATIONS

  scenario_data <-
    do.call(
      rbind,
      simulation_results
    )
  
  # 18. POINTWISE PERFORMANCE
  
  pointwise_table_temp <-
    calculate_pointwise_performance(
      scenario_data
    )
  
  
  # Add identification variables
  pointwise_table_temp$scenario <-
    "Scenario_1"
  
  pointwise_table_temp$n <-
    n
  
  pointwise_table_temp$kernel <-
    "Gaussian"
  
  
  # Store pointwise results
  pointwise_results[[counter]] <-
    pointwise_table_temp
  
  # 19. INTEGRATED PERFORMANCE
  
  integrated_table_temp <-
    calculate_integrated_performance(
      pointwise_table_temp
    )
  
  
  integrated_table_temp$scenario <-
    "Scenario_1"
  
  integrated_table_temp$n <-
    n
  
  integrated_table_temp$kernel <-
    "Gaussian"
  
  
  integrated_results[[counter]] <-
    integrated_table_temp
  
  # 20. BANDWIDTH RESULTS
  
  bandwidth_table_temp <-
    data.frame(
      replication = 1:B,
      bandwidth = selected_bandwidths,
      scenario = "Scenario_1",
      n = n,
      kernel = "Gaussian"
    )
  
  
  bandwidth_results[[counter]] <-
    bandwidth_table_temp
  
  
  counter <-
    counter + 1
  
  
  cat(
    "Completed Scenario_1 | n =",
    n,
    "| Gaussian kernel\n"
  )
}

# 21. COMBINE POINTWISE RESULTS

pointwise_table <-
  do.call(
    rbind,
    pointwise_results
  )

rownames(pointwise_table) <-
  NULL

# 22. COMBINE INTEGRATED RESULTS

integrated_table <-
  do.call(
    rbind,
    integrated_results
  )

rownames(integrated_table) <-
  NULL

# 23. COMBINE BANDWIDTH RESULTS

bandwidth_table <-
  do.call(
    rbind,
    bandwidth_results
  )

rownames(bandwidth_table) <-
  NULL

# 24. AVERAGE SELECTED BANDWIDTH

average_bandwidth <-
  aggregate(
    bandwidth ~ scenario + n,
    data = bandwidth_table,
    FUN = mean
  )

names(average_bandwidth)[3] <-
  "average_bandwidth"


# 25. FINAL MISE TABLE

MISE_table <-
  integrated_table[
    ,
    c(
      "scenario",
      "n",
      "kernel",
      "MISE",
      "RMISE",
      "Mean_Bias"
    )
  ]

# 26. DISPLAY POINTWISE RESULTS
cat("\n")
cat("\n")
cat("POINTWISE RESULTS\n")
cat("\n")

print(
  pointwise_table[
    ,
    c(
      "scenario",
      "n",
      "z",
      "theta_true",
      "mean_theta_hat",
      "variance",
      "bias",
      "mse",
      "rmse"
    )
  ]
)

# 27. DISPLAY INTEGRATED RESULTS

cat("\n")
cat("\n")
cat("INTEGRATED RESULTS\n")
cat("\n")

print(
  integrated_table[
    ,
    c(
      "scenario",
      "n",
      "kernel",
      "ISE",
      "MISE",
      "RMISE",
      "Mean_Bias"
    )
  ]
)

# 28. DISPLAY BANDWIDTH RESULTS

cat("\n")
cat("\n")
cat("AVERAGE SELECTED BANDWIDTH\n")
cat("\n")

print(average_bandwidth)

# 29. FINAL MISE TABLE

cat("\n")
cat("\n")
cat("FINAL MISE TABLE\n")
cat("\n")

print(MISE_table)

# 30. SAVE RESULTS
write.csv(
  pointwise_table,
  file.path(
    output_folder,
    "pointwise_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  integrated_table,
  file.path(
    output_folder,
    "integrated_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  bandwidth_table,
  file.path(
    output_folder,
    "bandwidth_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  average_bandwidth,
  file.path(
    output_folder,
    "average_bandwidth.csv"
  ),
  row.names = FALSE
)


write.csv(
  MISE_table,
  file.path(
    output_folder,
    "MISE_table.csv"
  ),
  row.names = FALSE
)

# 31. PLOT: TRUE VS ESTIMATED THETA

# Use n = 50,100,200,500 for illustration

plot_data <-
  subset(
    pointwise_table,
    n == 500
  )


plot(
  plot_data$z,
  plot_data$theta_true,
  type = "l",
  lwd = 2,
  xlab = "z",
  ylab = expression(theta(z)),
  main =
    ""
)

lines(
  plot_data$z,
  plot_data$mean_theta_hat,
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c(
    "True theta(z)",
    "Estimated theta(z)"
  ),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)



# 32. PLOT: POINTWISE RMSE

plot(
  plot_data$z,
  plot_data$rmse,
  type = "b",
  pch = 19,
  xlab = "z",
  ylab = "RMSE",
  main =
    ""
)



# 33. PLOT: MISE VS SAMPLE SIZE

mise_plot_data <-
  subset(
    integrated_table,
    scenario == "Scenario_1"
  )


plot(
  mise_plot_data$n,
  mise_plot_data$MISE,
  type = "b",
  pch = 19,
  xlab = "Sample size",
  ylab = "MISE",
  main =
    ""
)



# 34. DISPLAY OBJECTS


cat("\n")
cat("\n")
cat("SIMULATION COMPLETED SUCCESSFULLY\n")
cat("\n")

cat("\nObjects available:\n")
cat("pointwise_table\n")
cat("integrated_table\n")
cat("bandwidth_table\n")
cat("MISE_table\n")
cat("average_bandwidth\n")

cat("\nResults saved in:\n")
cat(output_folder, "\n")

cat("\nFiles created:\n")
cat("1. pointwise_results.csv\n")
cat("2. integrated_results.csv\n")
cat("3. bandwidth_results.csv\n")
cat("4. average_bandwidth.csv\n")
cat("5. MISE_table.csv\n")


# SIMULATION STUDY
# Scenario 1: Constant theta(z) = 2
# Nadaraya-Watson Estimation with Epanechnikov Kernel
# Lindley Distribution with Covariate-Dependent theta(z)
rm(list = ls())

set.seed(12345)

# 1. SIMULATION SETTINGS
# Number of Monte Carlo replications
B <- 1000

# Sample sizes
sample_sizes <- c(50, 100, 200, 500)

# Required evaluation points
evaluation_grid <- c(0.10, 0.30, 0.50, 0.70, 0.90)

# Candidate bandwidths
bandwidth_grid <- c(
  0.05, 0.075, 0.10, 0.125,
  0.15, 0.175, 0.20, 0.25, 0.30
)
# 2. TRUE THETA FUNCTION
# Scenario 1: Constant theta(z) = 2

theta_fun <- function(z) {
  2
}
# 3. LINDELY RANDOM NUMBER GENERATION

# Lindley density:
#
# f(x; theta) =
# theta^2/(theta+1) * (1+x) * exp(-theta*x)
#
# Mixture representation:
#
# Exponential(theta), probability theta/(theta+1)
# Gamma(2, theta), probability 1/(theta+1)


rlindley <- function(n, theta) {
  
  # If theta is a single value,
  # repeat it n times
  if (length(theta) == 1) {
    theta <- rep(theta, n)
  }
  
  # Check theta length
  if (length(theta) != n) {
    stop("Length of theta must be 1 or equal to n.")
  }
  
  # Check positivity
  if (any(theta <= 0)) {
    stop("Theta must be positive.")
  }
  
  # Uniform random variables
  u <- runif(n)
  
  # Storage for X
  x <- numeric(n)
  
  # Exponential component
  ind_exp <- u <= theta / (theta + 1)
  
  if (sum(ind_exp) > 0) {
    
    x[ind_exp] <-
      rexp(
        sum(ind_exp),
        rate = theta[ind_exp]
      )
  }
  
  # Gamma component
  ind_gamma <- !ind_exp
  
  if (sum(ind_gamma) > 0) {
    
    x[ind_gamma] <-
      rgamma(
        sum(ind_gamma),
        shape = 2,
        rate = theta[ind_gamma]
      )
  }
  
  return(x)
}
# 4. EPANECHNIKOV KERNEL
# Epanechnikov kernel:
#
# K(u) = 3/4 * (1-u^2), |u| <= 1
#        0,              |u| > 1


epanechnikov_kernel <- function(u) {
  
  K <- ifelse(
    abs(u) <= 1,
    (3 / 4) * (1 - u^2),
    0
  )
  
  return(K)
}

# 5. NADARAYA-WATSON ESTIMATOR

nw_estimator <- function(
    x,
    z,
    z0,
    h) {
  
  # Epanechnikov kernel argument
  u <- (z0 - z) / h
  
  # Epanechnikov weights
  weights <- epanechnikov_kernel(u)
  
  # Denominator
  denominator <- sum(weights)
  
  # Numerical safeguard
  if (denominator <= 0) {
    return(mean(x))
  }
  
  # Nadaraya-Watson estimate
  estimate <-
    sum(weights * x) /
    denominator
  
  return(estimate)
}

# 6. ESTIMATE m(z) ON THE REQUIRED GRID

estimate_m_grid <- function(
    x,
    z,
    grid,
    h) {
  
  m_hat <- numeric(length(grid))
  
  for (j in seq_along(grid)) {
    
    m_hat[j] <-
      nw_estimator(
        x = x,
        z = z,
        z0 = grid[j],
        h = h
      )
  }
  
  return(m_hat)
}

# 7. RECOVER THETA(z) FROM m(z)

# Lindley conditional mean:
#
# m(z) =
# [theta(z) + 2] /
# [theta(z) {theta(z) + 1}]
#
# Therefore:
#
# m theta^2 + (m - 1) theta - 2 = 0
#
# Positive root:
#
# theta =
# [1 - m + sqrt((m-1)^2 + 8m)] / (2m)


theta_from_m <- function(m) {
  
  # Avoid zero or negative values
  m <- pmax(m, 1e-10)
  
  theta_hat <-
    (
      1 - m +
        sqrt(
          (m - 1)^2 +
            8 * m
        )
    ) /
    (2 * m)
  
  # Ensure positive theta
  theta_hat <-
    pmax(theta_hat, 1e-10)
  
  return(theta_hat)
}
# 8. ESTIMATE THETA(z)

estimate_theta_grid <- function(
    x,
    z,
    grid,
    h) {
  
  # Estimate conditional mean
  m_hat <-
    estimate_m_grid(
      x = x,
      z = z,
      grid = grid,
      h = h
    )
  
  # Convert estimated mean to theta
  theta_hat <-
    theta_from_m(m_hat)
  
  return(theta_hat)
}
# 9. LEAVE-ONE-OUT CROSS-VALIDATION
cv_nw <- function(
    x,
    z,
    h) {
  
  n <- length(x)
  
  squared_errors <-
    numeric(n)
  
  for (i in seq_len(n)) {
    
    # Training data excluding observation i
    x_train <- x[-i]
    z_train <- z[-i]
    
    # Predict X_i at Z_i
    prediction <-
      nw_estimator(
        x = x_train,
        z = z_train,
        z0 = z[i],
        h = h
      )
    
    # Squared prediction error
    squared_errors[i] <-
      (x[i] - prediction)^2
  }
  
  return(sum(squared_errors))
}

# 10. OPTIMAL BANDWIDTH
select_bandwidth <- function(
    x,
    z,
    bandwidths) {
  
  cv_values <-
    numeric(length(bandwidths))
  
  for (j in seq_along(bandwidths)) {
    
    cv_values[j] <-
      cv_nw(
        x = x,
        z = z,
        h = bandwidths[j]
      )
  }
  
  # Select bandwidth giving minimum CV
  best_index <-
    which.min(cv_values)
  
  best_h <-
    bandwidths[best_index]
  
  return(best_h)
}

# 11. TRAPEZOIDAL INTEGRATION
trapz <- function(x, y) {
  
  if (length(x) != length(y)) {
    stop("x and y must have same length.")
  }
  
  if (length(x) < 2) {
    return(NA)
  }
  
  value <-
    sum(
      diff(x) *
        (
          head(y, -1) +
            tail(y, -1)
        ) / 2
    )
  
  return(value)
}
# 12. POINTWISE PERFORMANCE

calculate_pointwise_performance <- function(
    simulation_data) {
  
  # Mean estimated theta
  mean_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = mean
    )
  
  names(mean_estimate)[2] <-
    "mean_theta_hat"
  
  # Variance
  variance_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = var
    )
  
  names(variance_estimate)[2] <-
    "variance"
  
  # True theta
  true_values <-
    unique(
      simulation_data[
        ,
        c("z", "theta_true")
      ]
    )
  
  # Merge results
  result <-
    merge(
      true_values,
      mean_estimate,
      by = "z"
    )
  
  result <-
    merge(
      result,
      variance_estimate,
      by = "z"
    )
  
  # Bias
  result$bias <-
    result$mean_theta_hat -
    result$theta_true
  
  # MSE
  result$mse <-
    result$bias^2 +
    result$variance
  
  # RMSE
  result$rmse <-
    sqrt(result$mse)
  
  return(result)
}

# 13. INTEGRATED PERFORMANCE

calculate_integrated_performance <- function(
    pointwise_result) {
  
  # Integrated squared error
  ISE <-
    trapz(
      pointwise_result$z,
      pointwise_result$mse
    )
  
  # Integrated MSE
  MISE <- ISE
  
  # Root MISE
  RMISE <- sqrt(MISE)
  
  # Mean Bias over evaluation points
  Mean_Bias <-
    mean(pointwise_result$bias)
  
  result <-
    data.frame(
      ISE = ISE,
      MISE = MISE,
      RMISE = RMISE,
      Mean_Bias = Mean_Bias
    )
  
  return(result)
}

# 14. CREATE RESULT DIRECTORY
output_folder <-
  "Scenario_1_Epanechnikov_NW_Simulation"

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# 15. STORAGE OBJECTS
pointwise_results <- list()

integrated_results <- list()

bandwidth_results <- list()

counter <- 1

# 16. MAIN SIMULATION LOOP
cat("\n")
cat("\n")
cat("SCENARIO 1: CONSTANT THETA(z) = 2\n")
cat("NADARAYA-WATSON WITH EPANECHNIKOV KERNEL\n")
cat("\n")


for (n in sample_sizes) {
  
  cat("\n")
  cat("\n")
  cat("Sample size:", n, "\n")
  cat("\n")
  
  # Storage for Monte Carlo replications
  simulation_results <-
    vector("list", B)
  
  selected_bandwidths <-
    numeric(B)
  
  
  # Monte Carlo replications
  for (b in seq_len(B)) {
    
    # Generate covariate Z
    Z <-
      runif(
        n,
        min = 0,
        max = 1
      )
    
    
    # Generate true theta(Z)
    # Scenario 1: theta(z) = 2
    
    theta_data <-
      theta_fun(Z)
    
    
    # Generate Lindley response X
    
    X <-
      rlindley(
        n = n,
        theta = theta_data
      )
    
    
    # Select optimal bandwidth
    
    h_opt <-
      select_bandwidth(
        x = X,
        z = Z,
        bandwidths = bandwidth_grid
      )
    
    
    selected_bandwidths[b] <-
      h_opt
    
    
    # Estimate theta on required evaluation points
    
    theta_hat <-
      estimate_theta_grid(
        x = X,
        z = Z,
        grid = evaluation_grid,
        h = h_opt
      )
    
    
    # True theta on evaluation grid
    
    theta_grid_true <-
      theta_fun(
        evaluation_grid
      )
    
    
    # Store replication
    
    simulation_results[[b]] <-
      data.frame(
        replication = b,
        z = evaluation_grid,
        theta_true = theta_grid_true,
        theta_hat = theta_hat,
        bandwidth = h_opt
      )
    
    
    # Progress
    
    if (b %% 100 == 0) {
      
      cat(
        "Replication:",
        b,
        "of",
        B,
        "\n"
      )
    }
  }
  
  # 17. COMBINE MONTE CARLO REPLICATIONS
  scenario_data <-
    do.call(
      rbind,
      simulation_results
    )
  
  # 18. POINTWISE PERFORMANCE
  pointwise_table_temp <-
    calculate_pointwise_performance(
      scenario_data
    )
  
  
  # Add identification variables
  
  pointwise_table_temp$scenario <-
    "Scenario_1"
  
  pointwise_table_temp$n <-
    n
  
  pointwise_table_temp$kernel <-
    "Epanechnikov"
  
  
  # Store pointwise results
  
  pointwise_results[[counter]] <-
    pointwise_table_temp
  
  # 19. INTEGRATED PERFORMANCE
  
  integrated_table_temp <-
    calculate_integrated_performance(
      pointwise_table_temp
    )
  
  
  integrated_table_temp$scenario <-
    "Scenario_1"
  
  integrated_table_temp$n <-
    n
  
  integrated_table_temp$kernel <-
    "Epanechnikov"
  
  
  integrated_results[[counter]] <-
    integrated_table_temp
  
  # 20. BANDWIDTH RESULTS
  bandwidth_table_temp <-
    data.frame(
      replication = 1:B,
      bandwidth = selected_bandwidths,
      scenario = "Scenario_1",
      n = n,
      kernel = "Epanechnikov"
    )
  
  
  bandwidth_results[[counter]] <-
    bandwidth_table_temp
  
  
  counter <-
    counter + 1
  
  
  cat(
    "Completed Scenario_1 | n =",
    n,
    "| Epanechnikov kernel\n"
  )
}

# 21. COMBINE POINTWISE RESULTS

pointwise_table <-
  do.call(
    rbind,
    pointwise_results
  )

rownames(pointwise_table) <-
  NULL

# 22. COMBINE INTEGRATED RESULTS

integrated_table <-
  do.call(
    rbind,
    integrated_results
  )

rownames(integrated_table) <-
  NULL


# 23. COMBINE BANDWIDTH RESULTS

bandwidth_table <-
  do.call(
    rbind,
    bandwidth_results
  )

rownames(bandwidth_table) <-
  NULL

# 24. AVERAGE SELECTED BANDWIDTH

average_bandwidth <-
  aggregate(
    bandwidth ~ scenario + n,
    data = bandwidth_table,
    FUN = mean
  )

names(average_bandwidth)[3] <-
  "average_bandwidth"

# 25. FINAL MISE TABLE
MISE_table <-
  integrated_table[
    ,
    c(
      "scenario",
      "n",
      "kernel",
      "MISE",
      "RMISE",
      "Mean_Bias"
    )
  ]

# 26. DISPLAY POINTWISE RESULTS

cat("\n")
cat("\n")
cat("POINTWISE RESULTS\n")
cat("\n")

print(
  pointwise_table[
    ,
    c(
      "scenario",
      "n",
      "z",
      "theta_true",
      "mean_theta_hat",
      "variance",
      "bias",
      "mse",
      "rmse"
    )
  ]
)

# 27. DISPLAY INTEGRATED RESULTS

cat("\n")
cat("\n")
cat("INTEGRATED RESULTS\n")
cat("\n")

print(
  integrated_table[
    ,
    c(
      "scenario",
      "n",
      "kernel",
      "ISE",
      "MISE",
      "RMISE",
      "Mean_Bias"
    )
  ]
)

# 28. DISPLAY BANDWIDTH RESULTS

cat("\n")
cat("\n")
cat("AVERAGE SELECTED BANDWIDTH\n")
cat("\n")

print(average_bandwidth)

# 29. FINAL MISE TABLE

cat("\n")
cat("\n")
cat("FINAL MISE TABLE\n")
cat("\n")

print(MISE_table)


# 30. SAVE RESULTS

write.csv(
  pointwise_table,
  file.path(
    output_folder,
    "pointwise_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  integrated_table,
  file.path(
    output_folder,
    "integrated_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  bandwidth_table,
  file.path(
    output_folder,
    "bandwidth_results.csv"
  ),
  row.names = FALSE
)


write.csv(
  average_bandwidth,
  file.path(
    output_folder,
    "average_bandwidth.csv"
  ),
  row.names = FALSE
)


write.csv(
  MISE_table,
  file.path(
    output_folder,
    "MISE_table.csv"
  ),
  row.names = FALSE
)

# 31. PLOT: TRUE VS ESTIMATED THETA
# Use n = 50,100,200,500 for illustration

plot_data <-
  subset(
    pointwise_table,
    n == 500
  )


plot(
  plot_data$z,
  plot_data$theta_true,
  type = "l",
  lwd = 2,
  xlab = "z",
  ylab = expression(theta(z)),
  main =
    ""
)

lines(
  plot_data$z,
  plot_data$mean_theta_hat,
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c(
    "True theta(z)",
    "Estimated theta(z)"
  ),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

# 32. PLOT: POINTWISE RMSE

plot(
  plot_data$z,
  plot_data$rmse,
  type = "b",
  pch = 19,
  xlab = "z",
  ylab = "RMSE",
  main =
    ""
)

# 33. PLOT: MISE VS SAMPLE SIZE


mise_plot_data <-
  subset(
    integrated_table,
    scenario == "Scenario_1"
  )


plot(
  mise_plot_data$n,
  mise_plot_data$MISE,
  type = "b",
  pch = 19,
  xlab = "Sample size",
  ylab = "MISE",
  main =
    ""
)

# 34. DISPLAY OBJECTS

cat("\n")
cat("\n")
cat("SIMULATION COMPLETED SUCCESSFULLY\n")
cat("\n")

cat("\nObjects available:\n")

cat("pointwise_table\n")
cat("integrated_table\n")
cat("bandwidth_table\n")
cat("MISE_table\n")
cat("average_bandwidth\n")

cat("\nResults saved in:\n")

cat(
  output_folder,
  "\n"
)

cat("\nFiles created:\n")

cat("1. pointwise_results.csv\n")
cat("2. integrated_results.csv\n")
cat("3. bandwidth_results.csv\n")
cat("4. average_bandwidth.csv\n")
cat("5. MISE_table.csv\n")