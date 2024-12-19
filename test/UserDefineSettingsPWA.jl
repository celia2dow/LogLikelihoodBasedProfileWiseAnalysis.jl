# Dowling (2024)
# Adapted from Murphy et al. (2023) 

# Calls MakeSyntheticData.jl to generate synthetic data and calls 
# PWAFunction.jl given a user-defined error and mechanistice model,
# and user-defined guesses and bounds for the parameters.

#######################################################################################
## Package installation
#######################################################################################
# Top level business
using Pkg
include("envActivate.jl")

# Define the path to the environment - edit
env_path = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth " * 
        "dynamics/Code/Murphy2023ErrorModels-main/LogLikelihoodBasedProfileWiseAnalysis" 

# Define the environment's required packages - edit
required_packages = ["Plots", "NLopt", "Interpolations", "Distributions", 
        "Roots", "LaTeXStrings", "CSV", "DataFrames", 
        "DifferentialEquations", "Random", "Measures",
        "StatsPlots", "Colors"]

# Activate the environment with the appropriate packages via the function with inputs:
# (activOption, envPath, requiredPackages)
envActivate(1, env_path, required_packages)

# Use the package
using LogLikelihoodBasedProfileWiseAnalysis

#######################################################################################
## User defined settings
#######################################################################################
## Seed the simulation if desired
seed_num = 1234

## Number of points between guess and bounds for loglikelihood exploration
npts = 40;

## Experiment name - edit
experi_name = "TwoCellPopulationsNormalNoise"

## Parameters
# Create the struct for information storage
store_variables = defineStruct() 
# store_variables stores parameter information in the following form: 
#       store_variables([name1, name 2, ...], [true_value1, true_value2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [initial_value1, initial_value2, ...])
# Or if true values are unknown:
#       store_variables([name1, name 2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [guess_value1, initial_value2, ...])
# Or if true values and bounds are unnecessary:
#       store_variables([name1, name 2, ...], [guess_value1, initial_value2, ...])

# Mechanistic model - edit
r1=1.0;
r2=0.5;
guess_r1 = 1.2;
guess_r2 = 0.8;
r1_lb = 0.5; r1_ub = 1.5;
r2_lb = 0.1; r2_ub = 0.9;
model_params = store_variables(     # Store the true model parameter values and their names
    ["R_{d}","ζ"], 
    [r1, r2], 
    [r1_lb, r2_lb], 
    [r1_ub, r2_ub], 
    [guess_r1,guess_r2]) 

# Noise model - edit
error_type = "Normal"               # Normal, Lognormal
Sd = 5;
guess_Sd = 8;
Sd_lb = 1.0; Sd_ub = 10.0;
noise_params = store_variables(     # Store the true noise parameter valuess and their names
    "σ_{N}", 
    Sd, 
    Sd_lb, 
    Sd_ub, 
    guess_Sd) 

# Initial conditions - edit
C10=100.0;
C20=25.0;
X_array = store_variables(          # Store the initial conditions and the names of each state variable X_i(t)
    ["C_{1}(t)", "C_{2}(t)"], 
    [C10, C20]) 

# Data time-points - edit
t_max = 2.0; 
times = LinRange(0.0,t_max, 16); # synthetic data measurement times at which to generate data points

## Mathematical model: system of differential equations (ODE or PDE) - edit # NEED TO ADD THE CAPABILITY TO WORK WITH PDEs TO THE WORKFLOW
# 1 state variable example:
# f(u,p,t) = (p[1]*u/3) * (1 - max(0, (1+p[2]/p[1]) * (u-p[3])^3/(u^3))) # p[1]=λ, p[2]=zeta, and p[3]=R_d
# More than 1 state variable example:
function DE!(du,u,p,t)  # ! is "bang" convention indicating that the function modifies its inputs
    r1,r2=p
    du[1]=-r1*u[1];
    du[2]=r1*u[1] - r2*u[2];
end
   
#######################################################################################
## Data loading or generation
#######################################################################################

# Generate seed number if not given
if !isdefined(Main, :seed_num)
    seed_num = rand(1000:9999)
end

# Is synthetic data allowed to be negative? true/false
noNegatives = true

# Generate data synthetically
Random.seed!(seed_num)
data = makeSyntheticData(times, X_array.initial_values, DE!, model_params.true_values, error_type, noise_params.true_values, noNegatives)
df_data = DataFrame(data, :auto)  # Use :auto to automatically name columns
CSV.write( "synthetic_data" * experi_name * string(seed_num) * ".csv", df_data)

# Or read a CSV file of data
# data = DataFrame(CSV.File("path/to/your/file.csv"))  # Update with the correct path

#######################################################################################
## Path Wise Analysis
#######################################################################################

pwaFunction(experi_name, times, X_array, DE!, model_params, error_type, noise_params, seed_num, data, npts)