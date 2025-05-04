# inst/demo.R

library(NotSBS)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magrittr)

# Pre-defined Parameters
Time <- 400
p <- 100
theta_coef <- c(max(0.5-10*log(Time)/Time,0.42),0.5,0.5+6*sqrt(Time)/Time)

# Data Generating Process
data_sim <- duan_dgp(Time = Time, p = p, type = c(4,1,2), dep = TRUE, seed = 900, theta_coef = theta_coef)
trim <- round(2 * log(Time))
threshold <- 7.5

# Generate Seeded Intervals
lbd <- round(dim(X_t)[2]^(max(2/5, 1 - min(1, log(dim(X_t)[1])/log(dim(X_t)[2])))) * log(dim(X_t)[2])^1.1)
intervals <- seeded_intervals(Time, minl = lbd)

# Candidates - via Seeded Intervals
results <- cusum.fts(data_sim$x, V.diag = TRUE, lrv = TRUE, trim = trim, lbd = lbd)
results <- results[results$est.cp > trim & results$est.cp < (Time - trim), ]
results <- results[order(results$ed - results$st, decreasing = FALSE), ]

# Single Plot (oracle - iterations; fixed - one-shot)
plots_oracle <- plot_NotSBS_iterations(results, Time, m = 3, trim = trim, threshold = NULL, method = "oracle", c(max(0.5-10*log(Time)/Time,0.42),0.5,0.5+6*sqrt(Time)/Time))
for (p in plots_oracle) print(p)

# Use the derived (NH,ST] to choose threshold
plots_fixed <- plot_NotSBS_iterations(results, Time, m = 3, trim = trim, threshold = threshold, method = "fixed", c(max(0.5-10*log(Time)/Time,0.42),0.5,0.5+6*sqrt(Time)/Time))
for (p in plots_fixed) print(p)


# NotSBS to derive the estimated change points
res <- NotSBS_not(F_hat = F_hat, m = 3, trim = round(2 * log(Time)), threshold = NULL, method = "oracle", V_shap = "diag", lbd = lbd)
res

