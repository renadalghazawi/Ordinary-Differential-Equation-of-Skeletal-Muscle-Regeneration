
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

data <- data.frame(
  time = c(0, 1, 2, 3.5, 5, 7),
  M = c(1239, 9797, 2466, 1960, 1323, 888),
  M1 = c(2589, 8237, 8000, 10628, 6399, 2432),
  M2 = c(930, 3341, 3707, 10350, 7000, 3000),
  N = c(509, 4500, 272, 129, 98, 77),
  QSC = c(2700, 136, 19, 616, 1595, 1507),
  ASC = c(664, 1134, 123, 1309, 1260, 950),
  Mc = c(66, 40, 16, 900, 611, 741),
  Mn = c(NA, NA, NA, NA, NA, 30000)  # Include Mn data point
)

params <- c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20", "c21", "c22")

# Your ODE function here
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
    c1 = runif(1, 0.8, 1),    # Unchanged
    c2 = runif(1, 2, 3),  # Unchanged
    c3 = runif(1, 8.1e-05, 3e-04), # Unchanged
    c4 = runif(1, 15, 20),      # May need to increase to boost M production
    c5 = runif(1, 3.5e-05, 5.5e-05), # Increased for faster M to M1 conversion
    c6 = runif(1, 0.05, 0.1),  # Significantly increased to deplete M1 faster
    c7 = runif(1, 9e-17, 1.1e-16),  # Unchanged
    c8 = runif(1, 0.1, 0.2),  # Significantly increased to deplete M2 faster
    c9 = runif(1, 2e-06, 7e-06), # Adjust if needed for ASC dynamics
    c10 = runif(1, 9e-05, 4.2e-04),  # Adjust if needed for ASC dynamics
    c11 = runif(1, 0.0000905, 0.0004),  # Adjust if needed for ASC dynamics
    c12 = runif(1, 0.0005, 0.0008),  # Adjust if needed for ASC dynamics
    c13 = runif(1, 2.5, 3.5),          # Adjust if needed for Mc dynamics
    c14 = runif(1, 2e-04, 3e-04),  # Adjust if needed for Md dynamics
    c15 = runif(1, 1e-05, 2e-05),  # Adjust if needed for N dynamics
    c16 = runif(1, 5, 10),      # Adjust if needed for M dynamics
    c17 = runif(1, 500, 600),      # Unchanged
    c18 = runif(1, 0.0001, 0.0003), # Adjust if needed for ASC dynamics
    c19 = runif(1, 0.01, 0.03),     # Unchanged
    c20 = runif(1, 0.000015, 0.000055), # Unchanged (fixed at 0)
    c21 = runif(1, 1e-05, 9e-05),  # Adjust if needed for Mn dynamics
    c22 = runif(1, 1e-04, 2e-04)   # Adjust if needed for ASC dynamics
  )
)

# Objective function to minimize
objective_func <- function(parms) {
  c7 <- 1e-16
  
  # Model fitting code
  modelfit <- ode(y = init, times = time_seq, func = ode_func, parms = parms, method = 'lsoda')
  modelfit <- as.data.frame(modelfit)
  
  # Calculate SSE for each variable
  errorN <- sum(((data$N - modelfit$N) / max(data$N, na.rm = TRUE))^2)
  errorM <- sum(((data$M - modelfit$M) / max(data$M, na.rm = TRUE))^2)
  errorM1 <- sum(((data$M1 - modelfit$M1) / max(data$M1, na.rm = TRUE))^2)
  errorM2 <- sum(((data$M2 - modelfit$M2) / max(data$M2, na.rm = TRUE))^2)
  errorQSC <- sum(((data$QSC - modelfit$QSC) / max(data$QSC, na.rm = TRUE))^2)
  errorASC <- sum(((data$ASC - modelfit$ASC) / max(data$ASC, na.rm = TRUE))^2)
  errorMc <- sum(((data$Mc - modelfit$Mc) / max(data$Mc, na.rm = TRUE))^2)
  errorMn <- ifelse(!is.na(data$Mn[6]), ((data$Mn[6] - modelfit$Mn[6]) / max(data$Mn, na.rm = TRUE))^2, 0)
  
  total_error <- errorN + errorM + errorM1 + errorM2 + errorQSC + errorASC + errorMc + errorMn
  return(total_error)
}

# Bounds for optimization
lower_bounds <- rep(0, length(params))
upper_bounds <- rep('inf', length(params))

# Best result initialization
best_result <- NULL
best_total_error <- Inf
for (guess in initial_guesses) {
  # Consider using a global optimization method here, e.g., Genetic Algorithms
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
best_total_error <- best_total_error
print(best_optimized_parms)
print(best_total_error)

time_seq <- seq(from = 0, to = 7, by = 1/500)

sol <- ode(y = init, times = time_seq, func = ode_func, parms = best_optimized_parms, method = 'lsoda')
sim_df <- as.data.frame(sol)

# Ensure both data frames have the same structure and time points
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
plot_grid <- do.call(grid.arrange, c(plots, ncol = 5, nrow = 2))

# Print the combined plot grid
print(plot_grid)


library(tidyr)
library(dplyr)

deMicheli_long <- pivot_longer(micheli, cols = -time, names_to = "variable", values_to = "micheli_value")
mckellar_long <- pivot_longer(data, cols = -time, names_to = "variable", values_to = "mckellar_value")
sim_long <- pivot_longer(sim_df, cols = -time, names_to = "variable", values_to = "simulated_value")

# Joining De Micheli and McKellar data to simulated data
joined_data <- left_join(sim_long, deMicheli_long, by = c("time", "variable"))
joined_data <- left_join(joined_data, mckellar_long, by = c("time", "variable"))



plots <- lapply(unique(sim_long$variable), function(var) {
  plot_data <- filter(joined_data, variable == var)
  
  ggplot(plot_data, aes(x = time)) +
    geom_line(aes(y = simulated_value), color = "black", size = 0.5) +  # Increased line thickness
    geom_point(aes(y = mckellar_value), color = "blue", size = 1) +
    geom_point(aes(y = micheli_value), color = "red", size = 1) +
    geom_segment(aes(xend = time, y = simulated_value, yend = mckellar_value), 
                 color = "blue", linetype = "dashed", size = 0.5) +  # Increased line thickness
    geom_segment(aes(xend = time, y = simulated_value, yend = micheli_value), 
                 color = "red", linetype = "dashed", size = 0.5) +  # Increased line thickness
    labs(title = paste("", var), x = "DPI", y = "cells/mm3") +
    scale_color_manual(values = c("McKellar" = "blue", "De Micheli" = "red"), name = "Dataset") +
    scale_shape_manual(values = c("McKellar" = 1, "De Micheli" = 2), name = "Dataset") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),  # Increased title font size
      axis.title.x = element_text(size = 10),  # Increased x-axis title font size
      axis.title.y = element_text(size = 10),  # Increased y-axis title font size
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Increased x-axis text font size
      axis.text.y = element_text(size = 10),  # Increased y-axis text font size
      legend.position = "bottom"
    )
})

# Assume 'plots' is a list of ggplot objects
plot_grid <- grid.arrange(grobs = plots, ncol = 2, nrow = 5)

