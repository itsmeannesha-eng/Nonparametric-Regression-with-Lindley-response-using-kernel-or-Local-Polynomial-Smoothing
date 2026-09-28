# SIMULATION STUDY
# Scenario 6: Highly Nonlinear Parameter Function
#
# theta_6(z) = 1.25 + 1.5*z + 0.75*sin(2*pi*z)
#
# Nadaraya-Watson Estimation with Gaussian Kernel
# Lindley Distribution with Covariate-Dependent theta(z)


rm(list = ls())

set.seed(12345)


# 1. SIMULATION SETTINGS

B <- 1000

sample_sizes <- c(50, 100, 200, 500)

evaluation_grid <- c(
  0.10, 0.30, 0.50, 0.70, 0.90
)

bandwidth_grid <- c(
  0.05, 0.075, 0.10, 0.125,
  0.15, 0.175, 0.20, 0.25, 0.30
)


# 2. TRUE THETA FUNCTION

# Scenario VI:
#
# theta_6(z) =
# 1.25 + 1.5*z + 0.75*sin(2*pi*z)

theta_fun <- function(z) {
  
  1.25 +
    1.5 * z +
    0.75 * sin(2 * pi * z)
  
}


# 3. LINDELY RANDOM NUMBER GENERATION

rlindley <- function(n, theta) {
  
  if (length(theta) == 1) {
    
    theta <- rep(theta, n)
  }
  
  if (length(theta) != n) {
    
    stop(
      "Length of theta must be 1 or equal to n."
    )
  }
  
  if (any(theta <= 0)) {
    
    stop(
      "Theta must be positive."
    )
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
  
  K <-
    (1 / sqrt(2 * pi)) *
    exp(-0.5 * u^2)
  
  return(K)
}


# 5. NADARAYA-WATSON ESTIMATOR

nw_estimator <- function(
    x,
    z,
    z0,
    h) {
  
  
  u <-
    (z0 - z) / h
  
  
  weights <-
    gaussian_kernel(u)
  
  
  denominator <-
    sum(weights)
  
  
  if (denominator <= 0) {
    
    return(mean(x))
  }
  
  
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
  
  
  m_hat <-
    numeric(
      length(grid)
    )
  
  
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

theta_from_m <- function(m) {
  
  
  m <-
    pmax(
      m,
      1e-10
    )
  
  
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
    pmax(
      theta_hat,
      1e-10
    )
  
  
  return(theta_hat)
}

# 8. ESTIMATE THETA(z)

estimate_theta_grid <- function(
    x,
    z,
    grid,
    h) {
  
  
  m_hat <-
    estimate_m_grid(
      x = x,
      z = z,
      grid = grid,
      h = h
    )
  
  
  theta_hat <-
    theta_from_m(
      m_hat
    )
  
  
  return(theta_hat)
}

# 9. LEAVE-ONE-OUT CROSS-VALIDATION

cv_nw <- function(
    x,
    z,
    h) {
  
  
  n <-
    length(x)
  
  
  squared_errors <-
    numeric(n)
  
  
  for (i in seq_len(n)) {
    
    
    x_train <-
      x[-i]
    
    z_train <-
      z[-i]
    
    
    prediction <-
      nw_estimator(
        x = x_train,
        z = z_train,
        z0 = z[i],
        h = h
      )
    
    
    squared_errors[i] <-
      (
        x[i] -
          prediction
      )^2
  }
  
  
  return(
    sum(squared_errors)
  )
}


# 10. OPTIMAL BANDWIDTH
select_bandwidth <- function(
    x,
    z,
    bandwidths) {
  
  
  cv_values <-
    numeric(
      length(bandwidths)
    )
  
  
  for (j in seq_along(bandwidths)) {
    
    cv_values[j] <-
      cv_nw(
        x = x,
        z = z,
        h = bandwidths[j]
      )
  }
  
  
  best_index <-
    which.min(
      cv_values
    )
  
  
  best_h <-
    bandwidths[
      best_index
    ]
  
  
  return(best_h)
}


# 11. TRAPEZOIDAL INTEGRATION
trapz <- function(
    x,
    y) {
  
  
  if (
    length(x) !=
    length(y)
  ) {
    
    stop(
      "x and y must have same length."
    )
  }
  
  
  if (
    length(x) < 2
  ) {
    
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
  
  
  mean_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = mean
    )
  
  
  names(mean_estimate)[2] <-
    "mean_theta_hat"
  
  
  variance_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = var
    )
  
  
  names(variance_estimate)[2] <-
    "variance"
  
  
  true_values <-
    unique(
      simulation_data[
        ,
        c(
          "z",
          "theta_true"
        )
      ]
    )
  
  
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
  
  
  result$bias <-
    result$mean_theta_hat -
    result$theta_true
  
  
  result$mse <-
    result$bias^2 +
    result$variance
  
  
  result$rmse <-
    sqrt(
      result$mse
    )
  
  
  return(result)
}

# 13. INTEGRATED PERFORMANCE
calculate_integrated_performance <- function(
    pointwise_result) {
  
  
  ISE <-
    trapz(
      pointwise_result$z,
      pointwise_result$mse
    )
  
  
  MISE <-
    ISE
  
  
  RMISE <-
    sqrt(
      MISE
    )
  
  
  Mean_Bias <-
    mean(
      pointwise_result$bias
    )
  
  
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
  "Scenario_6_Highly_Nonlinear_Gaussian_NW_Simulation"


if (
  !dir.exists(
    output_folder
  )
) {
  
  dir.create(
    output_folder
  )
}


# 15. STORAGE OBJECTS

pointwise_results <-
  list()


integrated_results <-
  list()


bandwidth_results <-
  list()


counter <-
  1


# 16. MAIN SIMULATION LOOP

cat("\n")
cat("\n")

cat(
  "SCENARIO 6: HIGHLY NONLINEAR THETA(z)\n"
)

cat(
  "theta(z) = 1.25 + 1.5*z + 0.75*sin(2*pi*z)\n"
)

cat(
  "NADARAYA-WATSON WITH GAUSSIAN KERNEL\n"
)

cat("\n")


for (
  n in sample_sizes
) {
  
  
  cat("\n")
  cat("\n")
  
  cat(
    "Sample size:",
    n,
    "\n"
  )
  
  cat("\n")
  
  
  simulation_results <-
    vector(
      "list",
      B
    )
  
  
  selected_bandwidths <-
    numeric(B)
  
  
  for (
    b in seq_len(B)
  ) {
    
    
    # Generate Z
    
    Z <-
      runif(
        n,
        min = 0,
        max = 1
      )
    
    
    # Generate theta(Z)
    
    theta_data <-
      theta_fun(
        Z
      )
    
    
    # Generate Lindley observations
    
    X <-
      rlindley(
        n = n,
        theta = theta_data
      )
    
    
    # Select bandwidth
    
    h_opt <-
      select_bandwidth(
        x = X,
        z = Z,
        bandwidths = bandwidth_grid
      )
    
    
    selected_bandwidths[b] <-
      h_opt
    
    
    # Estimate theta
    
    theta_hat <-
      estimate_theta_grid(
        x = X,
        z = Z,
        grid = evaluation_grid,
        h = h_opt
      )
    
    
    # True theta
    
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
    
    
    if (
      b %% 100 == 0
    ) {
      
      cat(
        "Replication:",
        b,
        "of",
        B,
        "\n"
      )
    }
  }
  
  
  # Combine replications
  
  scenario_data <-
    do.call(
      rbind,
      simulation_results
    )
  
  
  # Pointwise performance
  
  pointwise_table_temp <-
    calculate_pointwise_performance(
      scenario_data
    )
  
  
  pointwise_table_temp$scenario <-
    "Scenario_6"
  
  
  pointwise_table_temp$n <-
    n
  
  
  pointwise_table_temp$kernel <-
    "Gaussian"
  
  
  pointwise_results[[counter]] <-
    pointwise_table_temp
  
  
  # Integrated performance
  
  integrated_table_temp <-
    calculate_integrated_performance(
      pointwise_table_temp
    )
  
  
  integrated_table_temp$scenario <-
    "Scenario_6"
  
  
  integrated_table_temp$n <-
    n
  
  
  integrated_table_temp$kernel <-
    "Gaussian"
  
  
  integrated_results[[counter]] <-
    integrated_table_temp
  
  
  # Bandwidth results
  
  bandwidth_table_temp <-
    data.frame(
      replication = 1:B,
      bandwidth = selected_bandwidths,
      scenario = "Scenario_6",
      n = n,
      kernel = "Gaussian"
    )
  
  
  bandwidth_results[[counter]] <-
    bandwidth_table_temp
  
  
  counter <-
    counter + 1
  
  
  cat(
    "Completed Scenario_6 | n =",
    n,
    "| Gaussian kernel\n"
  )
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

# 20. AVERAGE SELECTED BANDWIDTH
average_bandwidth <-
  aggregate(
    bandwidth ~ scenario + n,
    data = bandwidth_table,
    FUN = mean
  )

names(
  average_bandwidth
)[3] <-
  "average_bandwidth"

# 21. FINAL MISE TABLE
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

# 22. DISPLAY POINTWISE RESULTS
cat("\n")
cat("\n")

cat(
  "POINTWISE RESULTS\n"
)

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

# 23. DISPLAY INTEGRATED RESULTS
cat("\n")
cat("\n")

cat(
  "INTEGRATED RESULTS\n"
)

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

# 24. DISPLAY AVERAGE BANDWIDTH
cat("\n")
cat("\n")

cat(
  "AVERAGE SELECTED BANDWIDTH\n"
)

cat("\n")

print(
  average_bandwidth
)

# 25. DISPLAY FINAL MISE TABLE
cat("\n")
cat("\n")

cat(
  "FINAL MISE TABLE\n"
)

cat("\n")

print(
  MISE_table
)

# 26. SAVE RESULTS
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


# 27. PLOT: TRUE VS ESTIMATED THETA
#n=50,100,200,500
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
  lty = c(
    1,
    2
  ),
  lwd = 2,
  bty = "n"
)


# 28. PLOT: POINTWISE RMSE


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



# 29. PLOT: MISE VS SAMPLE SIZE

mise_plot_data <-
  subset(
    integrated_table,
    scenario == "Scenario_6"
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


# 30. FINAL MESSAGE

cat("\n")
cat("\n")

cat(
  "SIMULATION COMPLETED SUCCESSFULLY\n"
)

cat("\n")

cat(
  "Objects available:\n"
)

cat(
  "pointwise_table\n"
)

cat(
  "integrated_table\n"
)

cat(
  "bandwidth_table\n"
)

cat(
  "MISE_table\n"
)

cat(
  "average_bandwidth\n"
)

cat("\n")

cat(
  "Results saved in:\n"
)

cat(
  output_folder,
  "\n"
)





# SIMULATION STUDY
# Scenario 6: Highly Nonlinear Parameter Function
#
# theta_6(z) = 1.25 + 1.5*z + 0.75*sin(2*pi*z)
#
# Nadaraya-Watson Estimation with Epanechnikov Kernel
# Lindley Distribution with Covariate-Dependent theta(z)

rm(list = ls())

set.seed(12345)

# 1. SIMULATION SETTINGS
B <- 1000

sample_sizes <- c(50, 100, 200, 500)

evaluation_grid <- c(
  0.10, 0.30, 0.50, 0.70, 0.90
)

bandwidth_grid <- c(
  0.05, 0.075, 0.10, 0.125,
  0.15, 0.175, 0.20, 0.25, 0.30
)

# 2. TRUE THETA FUNCTION
theta_fun <- function(z) {
  
  1.25 +
    1.5 * z +
    0.75 * sin(2 * pi * z)
  
}

# 3. LINDELY RANDOM NUMBER GENERATION
rlindley <- function(n, theta) {
  
  if (length(theta) == 1) {
    
    theta <- rep(theta, n)
  }
  
  if (length(theta) != n) {
    
    stop(
      "Length of theta must be 1 or equal to n."
    )
  }
  
  if (any(theta <= 0)) {
    
    stop(
      "Theta must be positive."
    )
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

# 4. EPANECHNIKOV KERNEL
epanechnikov_kernel <- function(u) {
  
  K <-
    ifelse(
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
  
  
  u <-
    (z0 - z) / h
  
  
  weights <-
    epanechnikov_kernel(u)
  
  
  denominator <-
    sum(weights)
  
  
  if (denominator <= 0) {
    
    return(mean(x))
  }
  
  
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
  
  
  m_hat <-
    numeric(
      length(grid)
    )
  
  
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
theta_from_m <- function(m) {
  
  m <-
    pmax(
      m,
      1e-10
    )
  
  
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
    pmax(
      theta_hat,
      1e-10
    )
  
  
  return(theta_hat)
}

# 8. ESTIMATE THETA(z)
estimate_theta_grid <- function(
    x,
    z,
    grid,
    h) {
  
  
  m_hat <-
    estimate_m_grid(
      x = x,
      z = z,
      grid = grid,
      h = h
    )
  
  
  theta_hat <-
    theta_from_m(
      m_hat
    )
  
  
  return(theta_hat)
}

# 9. LEAVE-ONE-OUT CROSS-VALIDATION
cv_nw <- function(
    x,
    z,
    h) {
  
  
  n <-
    length(x)
  
  
  squared_errors <-
    numeric(n)
  
  
  for (i in seq_len(n)) {
    
    
    x_train <-
      x[-i]
    
    z_train <-
      z[-i]
    
    
    prediction <-
      nw_estimator(
        x = x_train,
        z = z_train,
        z0 = z[i],
        h = h
      )
    
    
    squared_errors[i] <-
      (
        x[i] -
          prediction
      )^2
  }
  
  
  return(
    sum(squared_errors)
  )
}

# 10. OPTIMAL BANDWIDTH
select_bandwidth <- function(
    x,
    z,
    bandwidths) {
  
  
  cv_values <-
    numeric(
      length(bandwidths)
    )
  
  
  for (j in seq_along(bandwidths)) {
    
    cv_values[j] <-
      cv_nw(
        x = x,
        z = z,
        h = bandwidths[j]
      )
  }
  
  
  best_index <-
    which.min(
      cv_values
    )
  
  
  best_h <-
    bandwidths[
      best_index
    ]
  
  
  return(best_h)
}


# 11. TRAPEZOIDAL INTEGRATION
trapz <- function(
    x,
    y) {
  
  
  if (
    length(x) !=
    length(y)
  ) {
    
    stop(
      "x and y must have same length."
    )
  }
  
  
  if (
    length(x) < 2
  ) {
    
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
  
  
  mean_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = mean
    )
  
  
  names(mean_estimate)[2] <-
    "mean_theta_hat"
  
  
  variance_estimate <-
    aggregate(
      theta_hat ~ z,
      data = simulation_data,
      FUN = var
    )
  
  
  names(variance_estimate)[2] <-
    "variance"
  
  
  true_values <-
    unique(
      simulation_data[
        ,
        c(
          "z",
          "theta_true"
        )
      ]
    )
  
  
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
  
  
  result$bias <-
    result$mean_theta_hat -
    result$theta_true
  
  
  result$mse <-
    result$bias^2 +
    result$variance
  
  
  result$rmse <-
    sqrt(
      result$mse
    )
  
  
  return(result)
}

# 13. INTEGRATED PERFORMANCE
calculate_integrated_performance <- function(
    pointwise_result) {
  
  
  ISE <-
    trapz(
      pointwise_result$z,
      pointwise_result$mse
    )
  
  
  MISE <-
    ISE
  
  
  RMISE <-
    sqrt(
      MISE
    )
  
  
  Mean_Bias <-
    mean(
      pointwise_result$bias
    )
  
  
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
  "Scenario_6_Highly_Nonlinear_Epanechnikov_NW_Simulation"


if (
  !dir.exists(
    output_folder
  )
) {
  
  dir.create(
    output_folder
  )
}

# 15. STORAGE OBJECTS
pointwise_results <-
  list()


integrated_results <-
  list()


bandwidth_results <-
  list()


counter <-
  1

# 16. MAIN SIMULATION LOOP
cat("\n")
cat("\n")

cat(
  "SCENARIO 6: HIGHLY NONLINEAR THETA(z)\n"
)

cat(
  "theta(z) = 1.25 + 1.5*z + 0.75*sin(2*pi*z)\n"
)

cat(
  "NADARAYA-WATSON WITH EPANECHNIKOV KERNEL\n"
)

cat("\n")


for (
  n in sample_sizes
) {
  
  
  cat("\n")
  cat("\n")
  
  cat(
    "Sample size:",
    n,
    "\n"
  )
  
  cat("\n")
  
  
  simulation_results <-
    vector(
      "list",
      B
    )
  
  
  selected_bandwidths <-
    numeric(B)
  
  
  for (
    b in seq_len(B)
  ) {
    
    
    # Generate Z
    
    Z <-
      runif(
        n,
        min = 0,
        max = 1
      )
    
    
    # Generate theta(Z)
    
    theta_data <-
      theta_fun(
        Z
      )
    
    
    # Generate Lindley observations
    
    X <-
      rlindley(
        n = n,
        theta = theta_data
      )
    
    
    # Select bandwidth
    
    h_opt <-
      select_bandwidth(
        x = X,
        z = Z,
        bandwidths = bandwidth_grid
      )
    
    
    selected_bandwidths[b] <-
      h_opt
    
    
    # Estimate theta
    
    theta_hat <-
      estimate_theta_grid(
        x = X,
        z = Z,
        grid = evaluation_grid,
        h = h_opt
      )
    
    
    # True theta
    
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
    
    
    if (
      b %% 100 == 0
    ) {
      
      cat(
        "Replication:",
        b,
        "of",
        B,
        "\n"
      )
    }
  }
  
  
  # Combine replications
  
  scenario_data <-
    do.call(
      rbind,
      simulation_results
    )
  
  
  # Pointwise performance
  
  pointwise_table_temp <-
    calculate_pointwise_performance(
      scenario_data
    )
  
  
  pointwise_table_temp$scenario <-
    "Scenario_6"
  
  
  pointwise_table_temp$n <-
    n
  
  
  pointwise_table_temp$kernel <-
    "Epanechnikov"
  
  
  pointwise_results[[counter]] <-
    pointwise_table_temp
  
  
  # Integrated performance
  
  integrated_table_temp <-
    calculate_integrated_performance(
      pointwise_table_temp
    )
  
  
  integrated_table_temp$scenario <-
    "Scenario_6"
  
  
  integrated_table_temp$n <-
    n
  
  
  integrated_table_temp$kernel <-
    "Epanechnikov"
  
  
  integrated_results[[counter]] <-
    integrated_table_temp
  
  
  # Bandwidth results
  
  bandwidth_table_temp <-
    data.frame(
      replication = 1:B,
      bandwidth = selected_bandwidths,
      scenario = "Scenario_6",
      n = n,
      kernel = "Epanechnikov"
    )
  
  
  bandwidth_results[[counter]] <-
    bandwidth_table_temp
  
  
  counter <-
    counter + 1
  
  
  cat(
    "Completed Scenario_6 | n =",
    n,
    "| Epanechnikov kernel\n"
  )
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

# 20. AVERAGE SELECTED BANDWIDTH
average_bandwidth <-
  aggregate(
    bandwidth ~ scenario + n,
    data = bandwidth_table,
    FUN = mean
  )

names(
  average_bandwidth
)[3] <-
  "average_bandwidth"

# 21. FINAL MISE TABLE
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

# 22. DISPLAY POINTWISE RESULTS
cat("\n")
cat("\n")

cat(
  "POINTWISE RESULTS\n"
)

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

# 23. DISPLAY INTEGRATED RESULTS
cat("\n")
cat("\n")

cat(
  "INTEGRATED RESULTS\n"
)

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

# 24. DISPLAY AVERAGE BANDWIDTH
cat("\n")
cat("\n")

cat(
  "AVERAGE SELECTED BANDWIDTH\n"
)

cat("\n")

print(
  average_bandwidth
)

# 25. DISPLAY FINAL MISE TABLE
cat("\n")
cat("\n")

cat(
  "FINAL MISE TABLE\n"
)

cat("\n")

print(
  MISE_table
)

# 26. SAVE RESULTS
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

# 27. PLOT: TRUE VS ESTIMATED THETA
#n=50,100,200,500
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
  lty = c(
    1,
    2
  ),
  lwd = 2,
  bty = "n"
)


# 28. PLOT: POINTWISE RMSE

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

# 29. PLOT: MISE VS SAMPLE SIZE
mise_plot_data <-
  subset(
    integrated_table,
    scenario == "Scenario_6"
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



# 30. FINAL MESSAGE


cat("\n")
cat("\n")

cat(
  "SIMULATION COMPLETED SUCCESSFULLY\n"
)

cat("\n")

cat(
  "Objects available:\n"
)

cat(
  "pointwise_table\n"
)

cat(
  "integrated_table\n"
)

cat(
  "bandwidth_table\n"
)

cat(
  "MISE_table\n"
)

cat(
  "average_bandwidth\n"
)

cat("\n")

cat(
  "Results saved in:\n"
)

cat(
  output_folder,
  "\n"
)

