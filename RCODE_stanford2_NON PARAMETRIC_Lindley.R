# Fast Nonparametric Regression Evaluation on Lifetime Data
# Dataset: survival::stanford2 (Stanford Heart Transplant Study)
# Response X = time (Lifetime / Survival time in days)
# Covariate Z = age (Age of patient at transplant)

# 1. Load Package & Built-in Lifetime Data
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
library(survival)
data("stanford2")

# Clean complete cases
dat <- stanford2[complete.cases(stanford2$time, stanford2$age), ]
Z <- dat$age   # Covariate: Age (Thermal/Biological Stress)
X <- dat$time  # Response: Failure/Survival time (Days)
n <- length(X)

cat("Dataset: survival::stanford2\n")
cat("Analysis sample size:", n, "\n")
cat("Response (X): Survival Time (days) | Covariate (Z): Patient Age\n\n")

# 2. Vectorized Kernel & Estimation Functions
gaussian_kernel <- function(u) {
  (1 / sqrt(2 * pi)) * exp(-0.5 * u^2)
}

# Fast Vectorized Nadaraya-Watson Kernel Mean Estimator
nw_mean <- function(z_eval, Z_data, X_data, h) {
  diff_mat <- outer(z_eval, Z_data, "-") / h
  w_mat <- gaussian_kernel(diff_mat)
  w_sums <- rowSums(w_mat)
  w_sums[w_sums == 0] <- 1e-10
  (w_mat %*% X_data) / w_sums
}

# Integrated Squared Error via Trapezoidal Rule
trapz_ise <- function(z_cov, errors_sq) {
  z_range <- range(z_cov)
  if (diff(z_range) == 0) return(NA_real_)
  z_norm <- (z_cov - z_range[1]) / diff(z_range)
  ord <- order(z_norm)
  z_s <- z_norm[ord]
  e_s <- errors_sq[ord]
  sum(diff(z_s) * (e_s[-1] + e_s[-length(e_s)]) / 2)
}

# Helper function to compute paper metrics
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
  
  data.frame(
    Mean_Bias = round(mean_bias, 4),
    Variance  = round(variance, 4),
    MSE       = round(mse, 4),
    RMSE      = round(rmse, 4),
    MAE       = round(mae, 4),
    MCSE      = round(mcse, 4),
    ISE       = round(ise, 4),
    RMISE     = round(rmise, 4)
  )
}

# 3. Fast Vectorized LOOCV Bandwidth Selection
bandwidth_grid <- seq(max(diff(range(Z)) / 100, 0.5),
                      diff(range(Z)) / 2.5, length.out = 60)
cv_scores <- numeric(length(bandwidth_grid))

cat("Running fast LOOCV bandwidth optimization...\n")
dist_full <- outer(Z, Z, "-")

for (k in seq_along(bandwidth_grid)) {
  h_cand <- bandwidth_grid[k]
  w_mat <- gaussian_kernel(dist_full / h_cand)
  diag(w_mat) <- 0 # LOOCV condition
  w_sums <- rowSums(w_mat)
  w_sums[w_sums == 0] <- 1e-10
  m_loo <- (w_mat %*% X) / w_sums
  cv_scores[k] <- mean((X - m_loo)^2)
}
optimal_h <- bandwidth_grid[which.min(cv_scores)]
cat("Optimal bandwidth (h) found:", round(optimal_h, 4), "\n\n")

# 4. Model Estimations & Metric Calculations
m_hat_prop  <- nw_mean(Z, Z, X, optimal_h)
res_prop    <- compute_metrics(m_hat_prop, X, Z)

m_bar_fixed <- rep(mean(X), n)
res_fixed   <- compute_metrics(m_bar_fixed, X, Z)

fit_lm      <- lm(X ~ Z)
m_hat_lm    <- predict(fit_lm)
res_lm      <- compute_metrics(m_hat_lm, X, Z)

fit_quad    <- lm(X ~ Z + I(Z^2))
m_hat_quad  <- predict(fit_quad)
res_quad    <- compute_metrics(m_hat_quad, X, Z)

# 5. Summary & Comparison Tables
comparison_table <- rbind(
  data.frame(Model = "Proposed Gaussian Kernel (NW)", res_prop),
  data.frame(Model = "Constant-Mean Benchmark", res_fixed),
  data.frame(Model = "Linear Regression (OLS)", res_lm),
  data.frame(Model = "Quadratic Polynomial Regression", res_quad)
)

cat("==================================================\n")
cat("TABLE 1: Model Performance Comparison Metrics\n")
cat("==================================================\n")
print(comparison_table, row.names = FALSE)

# Table 2: Summary Statistics of Fitted Mean Curves
z_grid <- seq(min(Z), max(Z), length.out = 200)
m_grid_prop <- nw_mean(z_grid, Z, X, optimal_h)
m_grid_lm   <- predict(fit_lm, newdata = data.frame(Z = z_grid))
m_grid_quad <- predict(fit_quad, newdata = data.frame(Z = z_grid))

fit_summary <- function(y) {
  c(Min = min(y), Q1 = unname(quantile(y, 0.25)),
    Median = median(y), Mean = mean(y),
    Q3 = unname(quantile(y, 0.75)), Max = max(y), SD = sd(y))
}

fit_summary_table <- rbind(
  data.frame(Model = "Proposed Gaussian Kernel (NW)",
             as.data.frame(t(round(fit_summary(m_grid_prop), 4)))),
  data.frame(Model = "Constant-Mean Benchmark",
             as.data.frame(t(round(fit_summary(rep(mean(X), length(z_grid))), 4)))),
  data.frame(Model = "Linear Regression (OLS)",
             as.data.frame(t(round(fit_summary(m_grid_lm), 4)))),
  data.frame(Model = "Quadratic Polynomial Regression",
             as.data.frame(t(round(fit_summary(m_grid_quad), 4))))
)

cat("\n==================================================\n")
cat("TABLE 2: Summary Statistics of Fitted Mean Curves\n")
cat("==================================================\n")
print(fit_summary_table, row.names = FALSE)

# Table 3: Detailed Residual Distribution Analysis
calc_res_summary <- function(m_hat, X_obs) {
  e <- X_obs - m_hat
  c(Min = min(e), Q1 = unname(quantile(e, 0.25)),
    Median = median(e), Mean = mean(e),
    Q3 = unname(quantile(e, 0.75)), Max = max(e), SD = sd(e))
}

residual_table <- rbind(
  data.frame(Model = "Proposed Gaussian Kernel (NW)",
             as.data.frame(t(round(calc_res_summary(m_hat_prop, X), 4)))),
  data.frame(Model = "Constant-Mean Benchmark",
             as.data.frame(t(round(calc_res_summary(m_bar_fixed, X), 4)))),
  data.frame(Model = "Linear Regression (OLS)",
             as.data.frame(t(round(calc_res_summary(m_hat_lm, X), 4)))),
  data.frame(Model = "Quadratic Polynomial Regression",
             as.data.frame(t(round(calc_res_summary(m_hat_quad, X), 4))))
)
colnames(residual_table) <- c("Model", "Min", "Q1", "Median", "Mean", "Q3", "Max", "SD")

cat("\n==================================================\n")
cat("TABLE 3: Residual Distribution Analysis Across Models\n")
cat("==================================================\n")
print(residual_table, row.names = FALSE)

# 6. Advanced 6-Panel Visualizations
c_prop  <- "#E41A1C" 
c_fixed <- "#377EB8" 
c_lm    <- "#4DAF4A" 
c_quad  <- "#984EA3" 
c_grid  <- "gray88"

m_grid_fixed <- rep(mean(X), length(z_grid))

if (interactive()) {
  dev.new(width = 12, height = 8)
}
old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1), bg = "white")

# Panel 1: LOOCV Bandwidth Selection Curve
plot(bandwidth_grid, cv_scores, type = "l", col = "#FF7F00", lwd = 2.5,
     main = "LOOCV Bandwidth Selection", xlab = "Bandwidth (h)",
     ylab = "Mean CV Squared Error", panel.first = grid(col = c_grid))
points(optimal_h, min(cv_scores), col = c_prop, pch = 19, cex = 1.4)
abline(v = optimal_h, col = c_prop, lty = 2, lwd = 1.5)
legend("bottomright", legend = paste0("Optimal h = ", round(optimal_h, 3)),
       col = c_prop, lty = 2, bty = "n", cex = 0.8)

# Panel 2: Fitted Mean Regression Curves
plot(Z, X, pch = 21, bg = "gray70", col = "black", cex = 1,
     main = "Regression Fit Comparison",
     xlab = "Age at Stress/Transplant (Z)",
     ylab = "Lifetime / Survival Time (X)",
     panel.first = grid(col = c_grid))
lines(z_grid, m_grid_prop, col = c_prop, lwd = 3)
lines(z_grid, m_grid_fixed, col = c_fixed, lwd = 2, lty = 2)
lines(z_grid, m_grid_lm, col = c_lm, lwd = 2, lty = 4)
lines(z_grid, m_grid_quad, col = c_quad, lwd = 2, lty = 3)
legend("topright", legend = c("Observed Data", "Gaussian Kernel (NW)",
                              "Constant Mean", "OLS", "Quadratic"),
       col = c("black", c_prop, c_fixed, c_lm, c_quad),
       pch = c(21, NA, NA, NA, NA), pt.bg = "gray70",
       lty = c(NA, 1, 2, 4, 3), lwd = c(NA, 3, 2, 2, 2),
       bty = "n", cex = 0.65)

# Panel 3: Estimated Regression Function
plot(z_grid, m_grid_prop, type = "n",
     main = "Estimated Mean Function",
     xlab = "Age (Z)", ylab = expression(hat(m)(z)),
     ylim = range(c(m_grid_prop, m_grid_lm, m_grid_quad, m_grid_fixed)),
     panel.first = grid(col = c_grid))
polygon(c(z_grid, rev(z_grid)),
        c(m_grid_prop, rep(min(m_grid_prop), length(z_grid))),
        col = adjustcolor(c_prop, alpha.f = 0.12), border = NA)
lines(z_grid, m_grid_prop, col = c_prop, lwd = 3)
lines(z_grid, m_grid_fixed, col = c_fixed, lwd = 2, lty = 2)
lines(z_grid, m_grid_lm, col = c_lm, lwd = 2, lty = 4)
lines(z_grid, m_grid_quad, col = c_quad, lwd = 2, lty = 3)
legend("topright", legend = c("Gaussian Kernel (NW)", "Constant Mean", "OLS", "Quadratic"),
       col = c(c_prop, c_fixed, c_lm, c_quad),
       lty = c(1, 2, 4, 3), lwd = c(3, 2, 2, 2), bty = "n", cex = 0.7)

# Panel 4: Residuals vs Fitted Values
plot(m_hat_prop, X - m_hat_prop, pch = 19,
     col = adjustcolor(c_prop, 0.7), main = "Residuals vs Fitted",
     xlab = "Fitted Values", ylab = "Residuals",
     ylim = range(c(X - m_hat_prop, X - m_hat_lm, X - m_hat_quad)),
     panel.first = grid(col = c_grid))
abline(h = 0, lty = 2, col = "black", lwd = 1.5)

# Panel 5: Empirical Residual Density Comparison
dens_prop <- density(X - m_hat_prop)
dens_lm   <- density(X - m_hat_lm)
dens_quad <- density(X - m_hat_quad)
plot(dens_prop, col = c_prop, lwd = 2.5,
     main = "Residual Density",
     xlab = "Residuals", ylab = "Density",
     xlim = range(c(dens_prop$x, dens_lm$x, dens_quad$x)),
     ylim = c(0, max(c(dens_prop$y, dens_lm$y, dens_quad$y)) * 1.15),
     panel.first = grid(col = c_grid))
lines(dens_lm, col = c_lm, lwd = 2, lty = 4)
lines(dens_quad, col = c_quad, lwd = 2, lty = 3)
legend("topright", legend = c("Gaussian Kernel (NW)", "OLS", "Quadratic"),
       col = c(c_prop, c_lm, c_quad), lty = c(1, 4, 3),
       lwd = c(2.5, 2, 2), bty = "n", cex = 0.8)

# Panel 6: Pointwise Absolute Errors across Covariate Z
ord <- order(Z)
plot(Z[ord], abs(X[ord] - m_hat_prop[ord]), type = "b", pch = 19,
     col = c_prop, lwd = 1.5, main = "Pointwise Absolute Errors",
     xlab = "Age (Z)", ylab = "|Residual|",
     panel.first = grid(col = c_grid))
lines(Z[ord], abs(X[ord] - m_hat_lm[ord]), type = "b", pch = 17,
      col = c_lm, lwd = 1.5, lty = 4)
lines(Z[ord], abs(X[ord] - m_hat_quad[ord]), type = "b", pch = 18,
      col = c_quad, lwd = 1.5, lty = 3)
legend("topright", legend = c("Gaussian Kernel (NW)", "OLS", "Quadratic"),
       col = c(c_prop, c_lm, c_quad), pch = c(19, 17, 18),
       lty = c(1, 4, 3), bty = "n", cex = 0.8)

par(old_par)
