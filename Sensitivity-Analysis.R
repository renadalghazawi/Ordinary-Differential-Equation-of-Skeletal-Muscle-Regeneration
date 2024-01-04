# Load required libraries
library(ggplot2)
library(reshape2)
library(gridExtra)

# Initialize an empty data frame to store sensitivity analysis results
# This data frame will contain parameters, state variables, and their relative changes
sensitivity_data <- data.frame(
  Parameter = character(), 
  StateVariable = character(),
  RelativeChange = numeric()
)

# Run simulations for each parameter and state variable
for (param in names(parameters)) {
  # Run the baseline model with current parameters
  baseline_output <- run_model(parameters)
  
  # Apply perturbations: 10% increase (1.1) and decrease (0.9)
  for (perturb in c(0.9, 1.1)) {
    # Copy and perturb the current parameter
    perturbed_params <- parameters
    perturbed_params[param] <- perturbed_params[param] * perturb
    # Run the model with the perturbed parameter
    perturbed_output <- run_model(perturbed_params)
    
    # Calculate the relative change for each state variable
    for (var in names(initial_state)) {
      relative_change <- (perturbed_output[, var] - baseline_output[, var]) / baseline_output[, var]
      # Store the results in the sensitivity data frame
      sensitivity_data <- rbind(sensitivity_data, data.frame(
        Parameter = param, 
        StateVariable = var,
        RelativeChange = relative_change
      ))
    }
  }
}

# Create boxplots for each state variable to visualize the sensitivity
# Assuming 'sensitivity_data' data frame is already created

# Extract a list of unique state variables
state_vars <- unique(sensitivity_data$StateVariable)

# Loop through each state variable to create individual boxplots
for (state_var in state_vars) {
  # Filter the data for the current state variable
  filtered_data <- subset(sensitivity_data, StateVariable == state_var)
  
  # Create a boxplot for the current state variable
  p <- ggplot(filtered_data, aes(x = Parameter, y = RelativeChange, fill = Parameter)) +
    geom_boxplot() +
    labs(title = paste("Boxplot of Relative Changes in", state_var, "due to Parameter Perturbations"),
         x = "Parameter", y = "Relative Change") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Display the boxplot
  print(p)
}

# Create a comprehensive boxplot grid for all state variables
plots <- ggplot(sensitivity_data, aes(x = Parameter, y = RelativeChange, fill = Parameter)) +
  geom_boxplot() +
  facet_wrap(~ StateVariable, scales = "free", ncol = 2) +  # Organize into two columns
  labs(title = "Boxplot of Relative Changes in State Variables due to Parameter Perturbations",
       x = "Parameter", y = "Relative Change") +
  theme_minimal() +
  theme(legend.position = "none")

# Display the grid of boxplots
print(plots)


# Reshape the data for plotting
sensitivity_long <- melt(sensitivity_data, id.vars = c("Time", "Parameter", "PerturbationLevel"))

#Temporal analysis
## Create a heatmap for each state variable
state_variables <- names(initial_state)
plots <- lapply(state_variables, function(var) {
  heatmap_data <- subset(sensitivity_long, variable == var)
  
  ggplot(heatmap_data, aes(x = Time, y = Parameter, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0) +
    labs(title = paste("Sensitivity of", var, "to Parameters Over Time"), x = "Time", y = "Parameter") +
    theme_minimal()
})

# Display the heatmaps
do.call(grid.arrange, plots)
