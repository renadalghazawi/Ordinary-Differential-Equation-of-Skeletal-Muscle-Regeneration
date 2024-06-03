# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)
library(gridExtra)

# Data preparation
data <- data.frame(
  time = c(0, 1, 2, 3.5, 5, 7),
  M = c(1239, 9797, 2466, 1960, 1323, 888),
  M1 = c(2589, 8237, 8000, 10628, 6399, 2432),
  M2 = c(930, 3341, 3707, 10350, 7000, 3000),
  N = c(509, 4500, 272, 129, 98, 77),
  QSC = c(2700, 136, 19, 616, 1595, 1507),
  ASC = c(664, 1134, 123, 1309, 1260, 950),
  Mc = c(66, 40, 16, 900, 611, 741)
)

params <- c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")

# ODE function
ode_func <- function(time, state, params) {
  with(as.list(c(state, params)), {
    QSC_feedback <- ifelse(QSC > 2700, c19 * QSC * (QSC - 2700), 0)
    Mn_feedback <- ifelse(Mn > 30000, c21 * Mn * (Mn - 30000), 0)
    
    dN <- c1 * Md - c2 * N - c15 * N * Md
    dNd <- c15 * N * Md - c3 * M1 * Nd
    dM <- c4 * N - c5 * M * (Nd + Md) - c16 * M
    dM1 <- c5 * M * (Nd + Md) - c6 * M1 - (c17 / (c7 + Nd + Md)) * M1
    dM2 <- M1 * c17 / (c7 + Nd + Md) - c8 * M2
    
    dQSC <- -c9 * QSC * N - c18 * QSC * Md + c10 * M2 * ASC - QSC_feedback
    dASC <- c9 * QSC * N + c18 * QSC * Md - c10 * ASC * M2 + c11 * M1 * ASC - M2 * ASC * c12 - c22 * ASC
    dMc <- c12 * ASC * M2 - c13 * Mc - c20 * Mc
    
    dMn <- c13 * Mc - Mn_feedback
    dMd <- -c14 * Md * M1 - c15 * N * Md
    
    return(list(c(dN, dNd, dM, dM1, dM2, dQSC, dASC, dMc, dMn, dMd)))
  })
}

time_seq <- c(0, 1, 2, 3.5, 5, 7)

# Initial values
init <- c(N = 0, Nd = 0, M = 0, M1 = 0, M2 = 0, QSC = 2700, ASC = 0, Mc = 0, Mn = 0, Md = 30000)

# Initial guesses
initial_guesses <- list(
  c(
    c1 = runif(1, 0.8, 1),
    c2 = runif(1, 2, 3),
    c3 = runif(1, 8.1e-05, 3e-04),
    c4 = runif(1, 15, 20),
    c5 = runif(1, 3.5e-05, 5.5e-05),
    c6 = runif(1, 0.05, 0.1),
    c7 = runif(1, 9e-17, 1.1e-16),
    c8 = runif(1, 0.1, 0.2),
    c9 = runif(1, 2e-06, 7e-06),
    c10 = runif(1, 8e-05, 3.5e-04),
    c11 = runif(1, 6e-05, 0.000095),
    c12 = runif(1, 0.0005, 0.0008),
    c13 = runif(1, 2.5, 3.5),
    c14 = runif(1, 2e-04, 3e-04),
    c15 = runif(1, 1e-05, 2e-05),
    c16 = runif(1, 5, 10),
    c17 = runif(1, 500, 600),
    c18 = runif(1, 0.0001, 0.0003),
    c19 = runif(1, 0.01, 0.03),
    c20 = runif(1, 0.000015, 0.000055),
    c21 = runif(1, 1e-05, 9e-05),
    c22 = runif(1, 1e-04, 2e-04)
  )
)

# Objective function to minimize
objective_func <- function(parms) {
  c7 <- 1e-16
  
  modelfit <- ode(y = init, times = time_seq, func = ode_func, parms = parms, method = 'lsoda')
  modelfit <- as.data.frame(modelfit)
  
  errorN <- sum(((data$N - modelfit$N) / max(data$N))^2)
  errorM <- sum(((data$M - modelfit$M) / max(data$M))^2)
  errorM1 <- sum(((data$M1 - modelfit$M1) / max(data$M1))^2)
  errorM2 <- sum(((data$M2 - modelfit$M2) / max(data$M2))^2)
  errorQSC <- sum(((data$QSC - modelfit$QSC) / max(data$QSC))^2)
  errorASC <- sum(((data$ASC - modelfit$ASC) / max(data$ASC))^2)
  errorMc <- sum(((data$Mc - modelfit$Mc) / max(data$Mc))^2)
  
  total_error <- errorN + errorM + errorM1 + errorM2 + errorQSC + errorASC + errorMc
  return(total_error)
}

# Bounds for optimization
lower_bounds <- rep(0, length(params))
upper_bounds <- rep('inf', length(params))

# Best result initialization
best_result <- NULL
best_total_error <- Inf
for (guess in initial_guesses) {
  result <- optim(par = guess, fn = objective_func, method = "L-BFGS-B",
                  lower = lower_bounds, upper = upper_bounds,
                  control = list(maxit = 10000, factr = 1e-7, pgtol = 1e-7))
  total_error <- objective_func(result$par)
  
  if (total_error < best_total_error) {
    best_result <- result
    best_total_error <- total_error
  }
}

# Final results
best_optimized_parms <- best_result$par
print(best_optimized_parms)
print(best_total_error)

# Run simulation with best-fit parameters
time_seq <- seq(from = 0, to = 7, by = 1/365)
sol <- ode(y = init, times = time_seq, func = ode_func, parms = best_optimized_parms, method = 'lsoda')
sim_df <- as.data.frame(sol)

# Data transformation
data_long <- data %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "experimental_value")

sim_long <- sim_df %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "simulated_value")

# Join data frames and calculate residuals
joined_df <- left_join(sim_long, data_long, by = c("time", "variable")) %>%
  mutate(residual = simulated_value - experimental_value)

# Calculate standard deviation of residuals for each variable
std_devs <- joined_df %>%
  group_by(variable) %>%
  summarise(std_dev = sd(residual, na.rm = TRUE))

# Create plots for each variable
plots <- lapply(unique(data_long$variable), function(var) {
  var_std_dev <- std_devs %>% filter(variable == var) %>% pull(std_dev)
  
  ggplot(filter(joined_df, variable == var), aes(x = time)) +
    geom_line(aes(y = simulated_value), color = "blue") +
    geom_point(aes(y = experimental_value), color = "purple") +
    geom_ribbon(aes(ymin = simulated_value - var_std_dev, ymax = simulated_value + var_std_dev), alpha = 0.2, fill = "blue") +
    labs(title = paste("", var), x = "Time", y = "cells/mm3") +
    theme_minimal()
})

# Combine all the plots into one grid
plot_grid <- do.call(grid.arrange, c(plots, ncol = 2, nrow = 5))

# Print the combined plot grid
print(plot_grid)

# To see the plots one by one
for (plot in plots) {
  print(plot)
  readline(prompt = "Press [enter] to see the next plot")
}
