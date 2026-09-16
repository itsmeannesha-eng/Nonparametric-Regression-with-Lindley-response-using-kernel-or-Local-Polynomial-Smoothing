# SIMULATION STUDY
# Nadaraya-Watson Estimation with Gaussian Kernel
# Lindley Distribution with Covariate-Dependent theta(z)
rm(list = ls())
set.seed(12345)

# 1. SIMULATION SETTINGS
# Number of Monte Carlo replications
#B <- 50
# For final simulation:
 B <- 1000
# Sample sizes
sample_sizes <- c(50, 100, 200, 500)
# Evaluation grid
evaluation_grid <- seq(0.05, 0.95, by = 0.05)
# Candidate bandwidths
bandwidth_grid <- c(
  0.05, 0.075, 0.10, 0.125,
  0.15, 0.175, 0.20, 0.25, 0.30
)
# 2. TRUE THETA FUNCTIONS
theta_functions <- list(
  
  # Scenario 1: Constant
  Scenario_1 = function(z) {
    2
  },
  
  # Scenario 2: Increasing
  Scenario_2 = function(z) {
    1.5 + 2 * z
  },
  
  # Scenario 3: Decreasing
  Scenario_3 = function(z) {
    3.5 - 2 * z
  },
  
  # Scenario 4: Quadratic
  Scenario_4 = function(z) {
    1.5 + 2 * (z - 0.5)^2
  },
  
  # Scenario 5: Sinusoidal
  Scenario_5 = function(z) {
    2 + 0.75 * sin(2 * pi * z)
  },
  
  # Scenario 6: Nonlinear
  Scenario_6 = function(z) {
    1.25 + 1.5 * z +
      0.75 * sin(2 * pi * z)
  }
)
# 3. LINDELY RANDOM NUMBER GENERATION

# Lindley distribution:
#
# f(x;theta) =
# theta^2/(theta+1) * (1+x) * exp(-theta*x)
#
# Mixture representation:
#
# Exponential(theta), probability theta/(theta+1)
# Gamma(2, theta), probability 1/(theta+1)

rlindley <- function(n, theta) {
  
  if (length(theta) == 1) {
    theta <- rep(theta, n)
  }
  
  if (length(theta) != n) {
    stop("Length of theta must be 1 or equal to n.")
  }
  
  if (any(theta <= 0)) {
    stop("Theta must be positive.")
  }
  
  u <- runif(n)
  
  x <- numeric(n)
  
  # Exponential component
  ind_exp <-
    u <= theta / (theta + 1)
  
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
  
  denominator <- sum(weights)
  
  if (denominator <= 0) {
    
    return(mean(x))
  }
  
  estimate <-
    sum(weights * x) /
    denominator
  
  return(estimate)
}
# 6. ESTIMATE m(z) ON A GRID
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
# [theta(z)+2]/
# theta(z)[theta(z)+1]
# Therefore:
# m theta^2 + (m-1)theta - 2 = 0
# Positive root:
# theta =
# [1-m + sqrt((m-1)^2 + 8m)]/
#             2m

theta_from_m <- function(m) {
  
  # Avoid zero values
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
  
  # Convert m_hat to theta_hat
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
    
    # Training data excluding i
    x_train <- x[-i]
    z_train <- z[-i]
    
    # Prediction at Z_i
    prediction <-
      nw_estimator(
        x = x_train,
        z = z_train,
        z0 = z[i],
        h = h
      )
    
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
        (head(y, -1) +
           tail(y, -1)) / 2
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
  
  MISE <- ISE
  
  RMISE <- sqrt(MISE)
  
  result <-
    data.frame(
      ISE = ISE,
      MISE = MISE,
      RMISE = RMISE
    )
  
  return(result)
}
# 14. CREATE RESULT DIRECTORY

if (!dir.exists("NW_Gaussian_Simulation")) {
  
  dir.create("NW_Gaussian_Simulation")
}
# 15. STORAGE OBJECTS

pointwise_results <- list()

integrated_results <- list()

bandwidth_results <- list()

counter <- 1

# 16. MAIN SIMULATION LOOP
for (scenario_name in names(theta_functions)) {
  
  cat("\n")
  cat("\n")
  cat("Scenario:", scenario_name, "\n")
  cat("\n")
  
  
  # Select true theta function
  theta_fun <-
    theta_functions[[scenario_name]]
  
  
  for (n in sample_sizes) {
    
    cat("\n")
    cat("Sample size:", n, "\n")
    
    # Monte Carlo replications
    simulation_results <-
      vector("list", B)
    
    selected_bandwidths <-
      numeric(B)
    
    
    for (b in seq_len(B)) {
      
      # Generate covariate Z
      Z <-
        runif(
          n,
          min = 0,
          max = 1
        )
      
      # Generate true theta(Z)
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
          bandwidths =
            bandwidth_grid
        )
      
      
      selected_bandwidths[b] <-
        h_opt
      
      # Estimate theta on evaluation grid
      
      theta_hat <-
        estimate_theta_grid(
          x = X,
          z = Z,
          grid =
            evaluation_grid,
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
          theta_true =
            theta_grid_true,
          theta_hat =
            theta_hat,
          bandwidth =
            h_opt
        )
      
      
      # Progress
      if (b %% 10 == 0) {
        
        cat(
          "Replication:",
          b,
          "of",
          B,
          "\n"
        )
      }
    }
    
    # COMBINE REPLICATIONS
    
    scenario_data <-
      do.call(
        rbind,
        simulation_results
      )
    
    # POINTWISE PERFORMANCE
    
    pointwise_table <-
      calculate_pointwise_performance(
        scenario_data
      )
    
    
    # Add identification variables
    pointwise_table$scenario <-
      scenario_name
    
    pointwise_table$n <-
      n
    
    pointwise_table$kernel <-
      "Gaussian"
    
    
    # Store
    pointwise_results[[counter]] <-
      pointwise_table
    
    # INTEGRATED PERFORMANCE
    integrated_table_temp <-
      calculate_integrated_performance(
        pointwise_table
      )
    
    
    integrated_table_temp$scenario <-
      scenario_name
    
    integrated_table_temp$n <-
      n
    
    integrated_table_temp$kernel <-
      "Gaussian"
    
    
    integrated_results[[counter]] <-
      integrated_table_temp
    
    # BANDWIDTH RESULTS
    bandwidth_table_temp <-
      data.frame(
        bandwidth =
          selected_bandwidths,
        scenario =
          scenario_name,
        n = n,
        kernel =
          "Gaussian"
      )
    
    
    bandwidth_results[[counter]] <-
      bandwidth_table_temp
    
    
    counter <-
      counter + 1
    
    
    cat(
      "Completed:",
      scenario_name,
      "| n =",
      n,
      "| Gaussian kernel\n"
    )
  }
}

# 17. COMBINE POINTWISE RESULTS
pointwise_table <-
  do.call(
    rbind,
    pointwise_results
  )

rownames(pointwise_table) <-
  NULL


# 18. COMBINE INTEGRATED RESULTS
integrated_table <-
  do.call(
    rbind,
    integrated_results
  )

rownames(integrated_table) <-
  NULL

# 19. COMBINE BANDWIDTH RESULTS
bandwidth_table <-
  do.call(
    rbind,
    bandwidth_results
  )

rownames(bandwidth_table) <-
  NULL

# 20. SAVE RESULTS
write.csv(
  pointwise_table,
  "NW_Gaussian_Simulation/pointwise_results.csv",
  row.names = FALSE
)

write.csv(
  integrated_table,
  "NW_Gaussian_Simulation/integrated_results.csv",
  row.names = FALSE
)

write.csv(
  bandwidth_table,
  "NW_Gaussian_Simulation/bandwidth_results.csv",
  row.names = FALSE
)
# 21. DISPLAY POINTWISE RESULTS
cat("\n")
cat("\n")
cat("POINTWISE RESULTS\n")
cat("\n")

print(
  head(
    pointwise_table,
    20
  )
)
# 22. DISPLAY INTEGRATED RESULTS

cat("\n")
cat("\n")
cat("INTEGRATED RESULTS\n")
cat("=\n")

print(
  integrated_table[
    ,
    c(
      "scenario",
      "n",
      "kernel",
      "ISE",
      "MISE",
      "RMISE"
    )
  ]
)



# 23. DISPLAY BANDWIDTH RESULTS
cat("\n")
cat("\n")
cat("BANDWIDTH RESULTS\n")
cat("\n")

print(
  head(
    bandwidth_table,
    20
  )
)

# 24. AVERAGE SELECTED BANDWIDTH
average_bandwidth <-
  aggregate(
    bandwidth ~ scenario + n,
    data = bandwidth_table,
    FUN = mean
  )

names(average_bandwidth)[3] <-
  "average_bandwidth"

cat("\n")
cat("\n")
cat("AVERAGE SELECTED BANDWIDTH\n")
cat("\n")

print(
  average_bandwidth
)

# 25. FINAL MISE TABLE
MISE_table <-
  integrated_table[
    ,
    c(
      "scenario",
      "n",
      "MISE",
      "RMISE"
    )
  ]

cat("\n")
cat("\n")
cat("FINAL MISE TABLE\n")
cat("\n")

print(MISE_table)


# Save MISE table
write.csv(
  MISE_table,
  "NW_Gaussian_Simulation/MISE_table.csv",
  row.names = FALSE
)

# 26. PLOT:
# TRUE VS ESTIMATED THETA
plot_data <-
  subset(
    pointwise_table,
    scenario == "Scenario_2" &
      n == 200
  )


plot(
  plot_data$z,
  plot_data$theta_true,
  type = "l",
  lwd = 2,
  xlab = "z",
  ylab = expression(theta(z)),
  main =
    "True and Estimated theta(z)"
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

# 27. PLOT:
# POINTWISE RMSE
plot(
  plot_data$z,
  plot_data$rmse,
  type = "l",
  lwd = 2,
  xlab = "z",
  ylab = "RMSE",
  main =
    "Pointwise RMSE"
)
# 28. PLOT:
# MISE VS SAMPLE SIZE
mise_plot_data <-
  subset(
    integrated_table,
    scenario == "Scenario_2"
  )


plot(
  mise_plot_data$n,
  mise_plot_data$MISE,
  type = "b",
  pch = 19,
  xlab = "Sample size",
  ylab = "MISE",
  main =
    "MISE versus Sample Size"
)

# 29. FINAL MESSAGE

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
cat("NW_Gaussian_Simulation/\n")

cat("\nFiles created:\n")
cat("1. pointwise_results.csv\n")
cat("2. integrated_results.csv\n")
cat("3. bandwidth_results.csv\n")
cat("4. MISE_table.csv\n")

