# Nonparametric Regression Evaluation: Comparison Plots & Paper Metrics
# Dataset: Inbuilt R dataset 'cars' (Z = speed, X = dist)
# 1. Load Data
data(cars)
Z <- cars$speed
X <- cars$dist
n <- length(X)

# 2. Kernel, Estimation, and Integration Functions
gaussian_kernel <- function(u) {
  (1 / sqrt(2 * pi)) * exp(-0.5 * u^2)
}

# Nadaraya-Watson Kernel Mean Estimator m_h(z)
nw_mean <- function(z_eval, Z_data, X_data, h) {
  sapply(z_eval, function(z) {
    weights <- gaussian_kernel((z - Z_data) / h)
    if (sum(weights) == 0) return(mean(X_data))
    sum(weights * X_data) / sum(weights)
  })
}

# Recover Stochastic Parameter Theta(z) via Lindley Closed-Form Solution
recover_theta <- function(m_hat) {
  (- (m_hat - 1) + sqrt((m_hat - 1)^2 + 8 * m_hat)) / (2 * m_hat)
}

# Trapezoidal Rule for Integrated Squared Error (ISE) over normalized z in [0, 1]
trapz_ise <- function(z_cov, errors_sq) {
  z_norm <- (z_cov - min(z_cov)) / (max(z_cov) - min(z_cov))
  ord <- order(z_norm)
  z_s <- z_norm[ord]
  e_s <- errors_sq[ord]
  sum(diff(z_s) * (e_s[-1] + e_s[-length(e_s)]) / 2)
}

# Helper function to compute paper performance metrics (LogLik Omitted)
compute_metrics <- function(m_hat, X_obs, Z_cov) {
  residuals <- m_hat - X_obs
  
  mean_bias <- mean(residuals)
  variance  <- var(m_hat)
  mse       <- mean(residuals^2)
  rmse      <- sqrt(mse)
  mae       <- mean(abs(residuals))
  mcse      <- sd(m_hat) / sqrt(length(m_hat))
  ise       <- trapz_ise(Z_cov, residuals^2)
  rmise     <- sqrt(ise)
  
  return(data.frame(
    Mean_Bias = round(mean_bias, 4),
    Variance  = round(variance, 4),
    MSE       = round(mse, 4),
    RMSE      = round(rmse, 4),
    MAE       = round(mae, 4),
    MCSE      = round(mcse, 4),
    ISE       = round(ise, 4),
    RMISE     = round(rmise, 4)
  ))
}

# 3. Bandwidth Selection via Leave-One-Out Cross-Validation (LOOCV)
bandwidth_grid <- seq(0.5, 6.0, by = 0.1)
cv_scores <- numeric(length(bandwidth_grid))

for (k in seq_along(bandwidth_grid)) {
  h_cand <- bandwidth_grid[k]
  loocv_err <- numeric(n)
  for (i in 1:n) {
    m_loo <- nw_mean(Z[i], Z[-i], X[-i], h_cand)
    loocv_err[i] <- (X[i] - m_loo)^2
  }
  cv_scores[k] <- sum(loocv_err)
}
optimal_h <- bandwidth_grid[which.min(cv_scores)]

# 4. Model Estimations & Performance Metric Calculations

# Model 1: Proposed Nonparametric Lindley 
m_hat_prop <- nw_mean(Z, Z, X, optimal_h)
theta_prop <- recover_theta(m_hat_prop)
res_prop   <- compute_metrics(m_hat_prop, X, Z)

# Model 2: Conventional Fixed-Parameter Lindley
m_bar_fixed  <- rep(mean(X), n)
theta_fixed  <- recover_theta(mean(X))
res_fixed    <- compute_metrics(m_bar_fixed, X, Z)

#  Model 3: Standard Linear Regression (OLS) 
fit_lm   <- lm(X ~ Z)
m_hat_lm <- predict(fit_lm)
res_lm   <- compute_metrics(m_hat_lm, X, Z)

# Model 4: Gamma GLM (Log Link) 
fit_glm   <- glm(X ~ Z, family = Gamma(link = "log"))
m_hat_glm <- predict(fit_glm, type = "response")
res_glm   <- compute_metrics(m_hat_glm, X, Z)

# 5. Paper Metric & Summary Tables

# Table 1: Master Model Comparison Table
comparison_table <- rbind(
  data.frame(Model = "Proposed Nonparametric Lindley", res_prop),
  data.frame(Model = "Fixed-Parameter Lindley", res_fixed),
  data.frame(Model = "Linear Regression (OLS)", res_lm),
  data.frame(Model = "Gamma GLM (Log link)", res_glm)
)

cat("\n=\n")
cat("TABLE 1: Master Model Performance Comparison Metrics\n")
cat("=\n")
print(comparison_table, row.names = FALSE)

# Table 2: Parameter Summary Statistics for theta(z) Trajectory
theta_summary_table <- data.frame(
  Parameter_Model = c("Proposed Nonparametric hat(theta)(z)", "Fixed-Parameter hat(theta)"),
  Min             = c(round(min(theta_prop), 4), round(theta_fixed, 4)),
  Q1              = c(round(quantile(theta_prop, 0.25), 4), round(theta_fixed, 4)),
  Median          = c(round(median(theta_prop), 4), round(theta_fixed, 4)),
  Mean            = c(round(mean(theta_prop), 4), round(theta_fixed, 4)),
  Q3              = c(round(quantile(theta_prop, 0.75), 4), round(theta_fixed, 4)),
  Max             = c(round(max(theta_prop), 4), round(theta_fixed, 4)),
  SD              = c(round(sd(theta_prop), 4), 0.0000)
)

cat("\n=\n")
cat("TABLE 2: Summary Statistics of Estimated Scale Parameter theta(z)\n")
cat("=\n")
print(theta_summary_table, row.names = FALSE)

# Table 3: Detailed Residual Distribution Analysis
calc_res_summary <- function(m_hat, X_obs) {
  e <- X_obs - m_hat
  c(Min = min(e), Q1 = quantile(e, 0.25), Median = median(e),
    Mean = mean(e), Q3 = quantile(e, 0.75), Max = max(e), SD = sd(e))
}

residual_table <- rbind(
  data.frame(Model = "Proposed Nonparametric Lindley", t(round(calc_res_summary(m_hat_prop, X), 4))),
  data.frame(Model = "Fixed-Parameter Lindley",       t(round(calc_res_summary(m_bar_fixed, X), 4))),
  data.frame(Model = "Linear Regression (OLS)",       t(round(calc_res_summary(m_hat_lm, X), 4))),
  data.frame(Model = "Gamma GLM (Log link)",         t(round(calc_res_summary(m_hat_glm, X), 4)))
)
colnames(residual_table) <- c("Model", "Min", "Q1", "Median", "Mean", "Q3", "Max", "SD")

cat("\n=\n")
cat("TABLE 3: Residual Distribution Analysis Across Models\n")
cat("=\n")
print(residual_table, row.names = FALSE)

# 6. Advanced 6-Panel Visualizations
# Palette colors matching academic figures
c_prop  <- "#E41A1C" # Crimson Red (Proposed)
c_fixed <- "#377EB8" # Royal Blue (Fixed)
c_lm    <- "#4DAF4A" # Forest Green (OLS)
c_glm   <- "#984EA3" # Deep Purple (Gamma GLM)
c_grid  <- "gray88"

z_grid <- seq(min(Z), max(Z), length.out = 200)
m_grid_prop <- nw_mean(z_grid, Z, X, optimal_h)
m_grid_lm   <- predict(fit_lm, newdata = data.frame(Z = z_grid))
m_grid_glm  <- predict(fit_glm, newdata = data.frame(Z = z_grid), type = "response")
theta_grid_prop <- recover_theta(m_grid_prop)

dev.new(width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1), bg = "white")

# Panel 1: LOOCV Bandwidth Selection Curve 
plot(bandwidth_grid, cv_scores, type = "l", col = "#FF7F00", lwd = 2.5,
     main = "", xlab = "Bandwidth (h)",
     ylab = "Cross-Validation Score", panel.first = grid(col = c_grid))
points(optimal_h, min(cv_scores), col = c_prop, pch = 19, cex = 1.6)
abline(v = optimal_h, col = c_prop, lty = 2, lwd = 1.5)
text(optimal_h + 0.35, min(cv_scores) + (max(cv_scores) - min(cv_scores)) * 0.1,
     labels = paste0("Optimal h = ", optimal_h), col = c_prop, font = 2, cex = 0.9)

# Panel 2: Fitted Mean Regression Curves 
plot(Z, X, pch = 21, bg = "gray70", col = "black", cex = 1.1,
     main = "",
     xlab = "Covariate Z (Speed)", ylab = "Response X (Stopping Distance)",
     panel.first = grid(col = c_grid))
lines(z_grid, m_grid_prop, col = c_prop,  lwd = 3, lty = 1)
lines(z_grid, rep(mean(X), 200), col = c_fixed, lwd = 2, lty = 2)
lines(z_grid, m_grid_lm,   col = c_lm,    lwd = 2, lty = 4)
lines(z_grid, m_grid_glm,  col = c_glm,   lwd = 2, lty = 3)
legend("topleft", legend = c("Observed Data", "Proposed Nonparametric", "Fixed Lindley", "OLS", "Gamma GLM"),
       col = c("black", c_prop, c_fixed, c_lm, c_glm),
       pch = c(21, NA, NA, NA, NA), pt.bg = "gray70", lty = c(NA, 1, 2, 4, 3),
       lwd = c(NA, 3, 2, 2, 2), bty = "n", cex = 0.75)

# Panel 3: Trajectory of Scale Parameter Theta(z)
plot(z_grid, theta_grid_prop, type = "n",
     main = expression(bold(paste("", hat(theta)(z)))),
     xlab = "Covariate Z (Speed)", ylab = expression(hat(theta)(z)),
     ylim = c(0, max(theta_grid_prop) * 1.15), panel.first = grid(col = c_grid))

polygon(c(z_grid, rev(z_grid)), c(theta_grid_prop, rep(0, length(z_grid))),
        col = adjustcolor(c_prop, alpha.f = 0.15), border = NA)

lines(z_grid, theta_grid_prop, col = c_prop, lwd = 3)
abline(h = theta_fixed, col = c_fixed, lwd = 2, lty = 2)

legend("topright", 
       legend = c(expression(paste("Proposed ", hat(theta)(z))),
                  expression(paste("Fixed Parameter ", hat(theta)))),
       col = c(c_prop, c_fixed), lty = c(1, 2), lwd = c(3, 2), bty = "n", cex = 0.8)

# Panel 4: Residuals vs. Fitted Values Diagnostic 
plot(m_hat_prop, X - m_hat_prop, pch = 19, col = adjustcolor(c_prop, 0.7),
     main = "", xlab = "Fitted Values m_hat(Z)",
     ylab = "Residuals (X - m_hat)", ylim = range(X - m_bar_fixed),
     panel.first = grid(col = c_grid))
points(m_hat_lm, X - m_hat_lm, pch = 17, col = adjustcolor(c_lm, 0.6))
points(m_hat_glm, X - m_hat_glm, pch = 18, col = adjustcolor(c_glm, 0.6))
abline(h = 0, lty = 2, col = "black", lwd = 1.5)
legend("topleft", legend = c("Proposed", "OLS", "Gamma GLM"),
       col = c(c_prop, c_lm, c_glm), pch = c(19, 17, 18), bty = "n", cex = 0.8)

# Panel 5: Empirical Residual Density Comparison 
dens_prop <- density(X - m_hat_prop)
dens_lm   <- density(X - m_hat_lm)
dens_glm  <- density(X - m_hat_glm)

plot(dens_prop, col = c_prop, lwd = 2.5, main = "",
     xlab = "Residuals", ylab = "Density",
     xlim = range(c(dens_prop$x, dens_lm$x, dens_glm$x)),
     ylim = c(0, max(c(dens_prop$y, dens_lm$y, dens_glm$y)) * 1.15),
     panel.first = grid(col = c_grid))
polygon(dens_prop, col = adjustcolor(c_prop, 0.2), border = c_prop)
lines(dens_lm, col = c_lm, lwd = 2, lty = 4)
lines(dens_glm, col = c_glm, lwd = 2, lty = 3)
legend("topright", legend = c("Proposed Nonparametric", "OLS", "Gamma GLM"),
       col = c(c_prop, c_lm, c_glm), lty = c(1, 4, 3), lwd = c(2.5, 2, 2), bty = "n", cex = 0.8)

#  Panel 6: Pointwise Absolute Errors across Covariate Z 
plot(Z, abs(X - m_hat_prop), type = "b", pch = 19, col = c_prop, lwd = 1.5,
     main = "",
     xlab = "Covariate Z (Speed)", ylab = "|Residual|", panel.first = grid(col = c_grid))
lines(Z, abs(X - m_hat_lm), type = "b", pch = 17, col = c_lm, lwd = 1.5, lty = 4)
lines(Z, abs(X - m_hat_glm), type = "b", pch = 18, col = c_glm, lwd = 1.5, lty = 3)
legend("topleft", legend = c("Proposed", "OLS", "Gamma GLM"),
       col = c(c_prop, c_lm, c_glm), pch = c(19, 17, 18), lty = c(1, 4, 3), bty = "n", cex = 0.8)

par(mfrow = c(1, 1))

