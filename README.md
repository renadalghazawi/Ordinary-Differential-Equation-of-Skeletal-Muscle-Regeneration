# Ordinary Differential Equation of Skeletal Muscle Regeneration
These are scripts developed to construct a mathematical model of muscle regeneration informed and calibrated by scRNA-seq data. The specific focus is to develop a non-linear ODE model that accurately captures the behavior of critical cell populations involved in muscle repair, primarily myogenic lineage cells and immune cells. The main goal was to optimally fit the model to a single, normalized scRNA-seq dataset, ensuring a detailed portrayal of cellular dynamics during muscle regeneration.

# Mathematical Derivation of Equations
Each cell type population is represented by a variable and their change over time is expressed using an ODE. Rates of cell activities, such as their infiltration or coming into existence (influx), transforming into another type (differentiation/polarization), apoptosis or leaving (outflux), and affecting each other (signaling), are determined by rate constants represented by 'c', as well as numbers of cells performing or regulating the activity.
![image](https://github.com/renadalghazawi/Ordinary-Differential-Equation-of-Skeletal-Muscle-Regeneration/assets/57923969/d608f9e4-217d-489c-8992-76aa9d069ed4)

# Optimization
Parameters representing regulatory interactions were initially estimated through manual iterative refinement, informed by expected outcomes guided from published experimental data. To refine these estimates, the L-BFGS-B algorithm was implemented within a NLLS framework. The objective is to minimize the sum of squared errors (SSE) between observed data and model predictions normalized by the maximum observed value. 

The optimization step was an iterative process implemented through a loop structure. This loop executed the following steps for each set of initial parameter guesses: 

1. Parameter Feeding: Initial parameter guesses are inputted into the optimization function.
2. Model Fitting: The `optim` function calls `objective_func`, which integrates the ODE   model over the specified time sequence using the current set of parameters. This step generates model predictions for each state variable.
3. Error computation: `objective_func` then calculates the total error as the SSE between the model predictions and the observed data for all state variables.
4. Result Evaluation: Upon completion of each iteration, the total error is evaluated. If this error is lower than that from any previous iteration, the current set of parameters is considered the best fit so far.
Mathematical Framework of Cellular Dynamics
5. Iteration and output: The loop tests new parameter guesses until convergence, finalizing the set with the lowest total error.
