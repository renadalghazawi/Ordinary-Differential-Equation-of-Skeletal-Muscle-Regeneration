# Load libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Sample dataset structure (replace with scaled single-cell data)
## Replace NA values with calibrated scRNA-seq time course data 
data <- data.frame(
  time = c(0, 1, 2, 3.5, 5, 7),
  M1 = c(NA, NA, NA, NA, NA, NA),
  M2 = c(NA, NA, NA, NA, NA, NA),
  M = c(NA, NA, NA, NA, NA, NA),
  N = c(NA, NA, NA, NA, NA, NA),
  QSC = c(NA, NA, NA, NA, NA, NA),
  ASC = c(NA, NA, NA, NA, NA, NA),
  Mc = c(NA, NA, NA, NA, NA, NA)
)


# Define parameters from ODE system 
params <- c("cMd_out", "cNMd", "cN_in", "cN_out", "cM1Nd", "cM_in", "cMM1", "cM_out", "cM1_out", 
                        "cM2_inhib", "cM2_out", "cQSCN", "cQSCMd", "cASCM2", "cQSC_self", "cASC_pro", "cASC_out",
                        "cSC_diff", "cMc_fusion", "cMc_out", "cM1M2", "cMn_self")

# Define the ODE model function
ode_func <- function(time, state, params) {
  with(as.list(c(state, params)), { 
    
    ## Feedback mechanisms
    QSC_feedback <- ifelse(QSC > 2700, cQSC_self * QSC * (QSC - 2700), 0)
    Mn_feedback <- ifelse(Mn > 30000, cMn_self * Mn * (Mn - 30000), 0)
    
    ## Immune cells
    dN <- cN_in * Md - cN_out * N - cNMd * N * Md 
    dNd <- cNMd * N * Md - cM1Nd * M1 * Nd
    dM <- cM_in * N - cMM1 * M * (Nd + Md) - cM_out * M
    dM1 <- cMM1 * M * (Nd + Md) - cM1_out * M1 - (cM1M2 / (cM2_inhib + Nd + Md)) * M1
    dM2 <- M1 * cM1M2 / (cM2_inhib + Nd + Md) - cM2_out * M2
    
    ## MuSCs
    dQSC <- -cQSCN * QSC * N - cQSCMd * QSC * Md + cASCM2 * M2 * ASC - QSC_feedback
    dASC <- cQSCN * QSC * N + cQSCMd * QSC * Md - cASCM2 * ASC * M2 + cASC_pro * M1 * ASC - M2 * ASC * cSC_diff - cASC_out * ASC
    dMc <- cSC_diff * ASC * M2 - cMc_fusion * Mc - cMc_out * Mc
    
    ## Myonuclei
    dMn <- cMc_fusion * Mc - Mn_feedback
    dMd <- -cMd_out * Md * M1 - cNMd * N * Md
    
    return(list(c(dN, dNd, dM, dM1, dM2, dQSC, dASC, dMc, dMn, dMd)))
  })
}


# Initial values and time sequence for the model
init <- c(N = 0, Nd = 0, M = 0, M1 = 0, M2 = 0, QSC = 2700, ASC = 0, Mc = 0, Mn =0, Md = 30000)
time_seq <- c(0, 1, 2, 3.5, 5, 7)  # Time points of the experimental dataset

# Define the objective function for optimization
objective_func <- function(parms) {
  
  # Run the model and calculate the error between model output and data
  cM2_inhib <- 1e-16 # Fixed Param
  
  # Model fitting 
  modelfit <- ode(y = init, times = time_seq, func = ode_func, parms = params, method = 'lsoda')
  modelfit <- as.data.frame(modelfit)
  
  # Calculate SSE for each variable
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


# Initial parameter guesses (use realistic ranges)
initial_guesses <- runif(22, min = 0, max = 1)  # Random initial guesses

# Define bounds for optimization
lower_bounds <- rep(0, length(params))
upper_bounds <- rep(Inf, length(params))

# Perform optimization using the L-BFGS-B method
optim_result <- optim(par = initial_guesses, fn = objective_func, method = "L-BFGS-B",
                      lower = lower_bounds, upper = upper_bounds,
                      control = list(maxit = 10000, factr = 1e-7, pgtol = 1e-7))

# Output the results
best_parameters <- setNames(optim_result$par, params)
best_total_error <- optim_result$value
print(best_parameters)
print(best_total_error)

# Plotting the model results against the data
sim_df <- as.data.frame(ode(y = init, times = time_seq, func = ode_func, parms = best_parameters, method = 'lsoda'))
data_long <- pivot_longer(data, cols = -time, names_to = "variable", values_to = "experimental_value")
sim_long <- pivot_longer(sim_df, cols = -time, names_to = "variable", values_to = "simulated_value")

# Calculating and plotting residuals
joined_df <- left_join(sim_long, data_long, by = c("time", "variable"))
joined_df <- joined_df %>% mutate(residual = simulated_value - experimental_value)
std_devs <- joined_df %>% group_by(variable) %>% summarise(std_dev = sd(residual, na.rm = TRUE))

# Generate and display plots
plots <- lapply(unique(joined_df$variable), function(var) {
  var_data <- filter(joined_df, variable == var)
  ggplot(var_data, aes(x = time)) +
    geom_line(aes(y = simulated_value), color = "blue") +
    geom_point(aes(y = experimental_value), color = "purple") +
    geom_ribbon(aes(ymin = simulated_value - std_dev, ymax = simulated_value + std_dev), alpha = 0.2, fill = "blue") +
    labs(title = var, x = "Time", y = "Value") +
    theme_minimal()
})

# Arrange and display plots in a grid
plot_grid <- do.call(grid.arrange, c(plots, ncol = 2, nrow = ceiling(length(plots) / 2)))
print(plot_grid)
