# Dowling (2024)
# Adapted from Murphy et al. (2023) 

# Calls MakeSyntheticData.jl to generate synthetic data and calls 
# PWAFunction.jl given a user-defined error and mechanistice model,
# and user-defined guesses and bounds for the parameters.

#######################################################################################
## Package installation
#######################################################################################

# ## PACKAGE INSTALLATION OPTION 1
# # Activate the package environment
# using Pkg
# env_path = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth " * 
#            "dynamics/Code/Murphy2023ErrorModels-main/LogLikelihoodBasedProfileWiseAnalysis" 
#            # Change this depending on the location of the package in your directory
# Pkg.activate(env_path)  # Activate the environment
# # Define the environment's required packages
# required_packages = ["Plots", "NLopt", "Interpolations", "Distributions", 
#                      "Roots", "LaTeXStrings", "CSV", "DataFrames", 
#                      "DifferentialEquations", "Random", "Measures",
#                      "StatsPlots", "StructuralIdentifiability", "Colors"]
# # Get the current environment's package names
# installed_packages = keys(Pkg.installed())
# # Flag to track if any package is missing
# missing_packages = false
# # Check for each required package and add if not already in Project.toml
# for pkg in required_packages
#     if !(pkg in installed_packages)
#         global missing_packages = true
#         println("Missing package: $pkg. Installing...")  # Print the missing package name
#         if pkg == "StructuralIdentifiability"
#             Pkg.add(Pkg.PackageSpec(;name="StructuralIdentifiability", version="v0.5.1")) # need to use older version
#         else
#             Pkg.add(pkg)
#         end
#     else
#         println("$pkg is already installed.")  # Print the already installed package
#     end
# end

## PACKAGE INSTALLATION OPTION 2
import Pkg
env_path = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth " * 
           "dynamics/Code/Murphy2023ErrorModels-main/LogLikelihoodBasedProfileWiseAnalysis" 
           # Change this depending on the location of the package in your directory
Pkg.activate(env_path)  # Activate the environment
# Define the environment's required packages
required_packages = ["Plots", "NLopt", "Interpolations", "Distributions", 
                     "Roots", "LaTeXStrings", "CSV", "DataFrames", 
                     "DifferentialEquations", "Random", "Measures",
                     "StatsPlots", "StructuralIdentifiability", "Colors"]
# Ensure that the required packages are installed
Pkg.add(required_packages)

# Check that Manifest.toml is consistent with Project.tomls of the current project and all its
# dependencies, updating it if necessary, then instantiate. 
Pkg.resolve()
# Use the package
using LogLikelihoodBasedProfileWiseAnalysis

#######################################################################################
## User defined settings
#######################################################################################

## Parameters
struct store_variables # DO NOT EDIT
    names::Union{Vector{String}, String}                    # String for names
    true_values::Union{Vector{Float64}, Float64, Int64, Nothing}   # True values if known
    lower_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}  # Lower bounds
    upper_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}  # Upper bounds
    initial_values::Union{Vector{Float64}, Float64, Int64}         # Initial guesses or initial conditions
    # Constructor function for all five fields
    function store_variables(names::Union{Vector{String}, String}, 
        true_values::Union{Vector{Float64}, Float64, Int64, Nothing},
        lower_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}, 
        upper_bounds::Union{Vector{Float64}, Float64, Int64, Nothing},  
        initial_values::Union{Vector{Float64}, Float64, Int64})
        return new(names, true_values, lower_bounds, upper_bounds, initial_values)
    end
    # Constructor function for four fields
    function store_variables(names::Union{Vector{String}, String}, 
        lower_bounds::Union{Vector{Float64}, Float64, Int64, Nothing}, 
        upper_bounds::Union{Vector{Float64}, Float64, Int64, Nothing},  
        initial_values::Union{Vector{Float64}, Float64, Int64})
        return new(names, nothing, lower_bounds, upper_bounds, initial_values) # Assign nothing to true_values
    end
    # Constructor function for two fields (default true_values to nothing)
    function store_variables(names::Union{Vector{String}, String}, 
        initial_values::Union{Vector{Float64}, Float64, Int64})
        return new(names, nothing, nothing, nothing, initial_values)  # Assign nothing to true_values and bounds
    end
end
# Store in the following form: 
#       store_variables([name1, name 2, ...], [true_value1, true_value2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [initial_value1, initial_value2, ...])
# Or if true values are unknown:
#       store_variables([name1, name 2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [guess_value1, initial_value2, ...])
# Or if true values and bounds are unnecessary:
#       store_variables([name1, name 2, ...], [guess_value1, initial_value2, ...])

#######################################################################################
## Example 1: Linear Normal Noise ODE 
## (2 state variables, 2 model paramters, synthetic data, true paramters known)
#######################################################################################

## Seed the simulation if desired
seed_num = 1234

## Number of points between guess and bounds for loglikelihood exploration
npts = 40;

## Experiment name - edit
experi_name = "TwoCellPopulationsNormalNoise"

## Parameters
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

## Mathematical model: system of differential equations (ODE or PDE) - edit
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

# Generate data synthetically
Random.seed!(seed_num)
data = makeSyntheticData(times, X_array.initial_values, DE!, model_params.true_values, error_type, noise_params.true_values)
df_data = DataFrame(data, :auto)  # Use :auto to automatically name columns
CSV.write( "synthetic_data" * experi_name * string(seed_num) * ".csv", df_data)

# Or read a CSV file of data
# data = DataFrame(CSV.File("path/to/your/file.csv"))  # Update with the correct path

#######################################################################################
## Path Wise Analysis
#######################################################################################

pwaFunction(experi_name, times, X_array, DE!, model_params, error_type, noise_params, seed_num, data, npts)

#######################################################################################
## Example 2: Linear LogNormal Noise ODE 
## (2 state variables, 2 model paramters, synthetic data, true paramters known)
#######################################################################################

## Seed the simulation if desired
seed_num = 1

## Number of points between guess and bounds for loglikelihood exploration
npts = 40;

## Experiment name - edit
experi_name = "TwoCellPopulationsLognormalNoise"

## Parameters
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
error_type = "Lognormal"            # Normal, Lognormal
Sd = 0.4;
guess_Sd = 0.35;
Sd_lb = 0.2; Sd_ub = 0.8;
noise_params = store_variables(     # Store the true noise parameter valuess and their names
    "σ_{L}", 
    Sd, 
    Sd_lb, 
    Sd_ub, 
    guess_Sd) 

# Initial conditions - edit
C10=100.0;
C20=10.0;
X_array = store_variables(          # Store the initial conditions and the names of each state variable X_i(t)
    ["C_{1}(t)", "C_{2}(t)"], 
    [C10, C20]) 

# Data time-points - edit
t_max = 5.0; 
times = LinRange(0.0,t_max, 31); # synthetic data measurement times at which to generate data points

## Mathematical model: system of differential equations (ODE or PDE) - edit
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

# Generate data synthetically
Random.seed!(seed_num)
data = makeSyntheticData(times, X_array.initial_values, DE!, model_params.true_values, error_type, noise_params.true_values)
df_data = DataFrame(data, :auto)  # Use :auto to automatically name columns
CSV.write( "synthetic_data" * experi_name * string(seed_num) * ".csv", df_data)

# Or read a CSV file of data
# data = DataFrame(CSV.File("path/to/your/file.csv"))  # Update with the correct path

#######################################################################################
## Path Wise Analysis
#######################################################################################

pwaFunction(experi_name, times, X_array, DE!, model_params, error_type, noise_params, seed_num, data, npts)

#######################################################################################
## Example 3: Linear Normal Noise ODE 
## (1 state variable, 2 model paramters, synthetic data, true paramters known)
#######################################################################################

## Seed the simulation if desired
seed_num = 1234

## Number of points between guess and bounds for loglikelihood exploration
npts = 40;

## Experiment name - edit
experi_name = "RadialDeathModelNormalNoise"

## Parameters
# Mechanistic model - edit
r1=1.0;
r2=1.25;
r3=75.0;
guess_r1 = 1.05;
guess_r2 = 1.1;
guess_r3 = 70.0;
r1_lb = 0.9; r1_ub = 1.1;
r2_lb = 1.0; r2_ub = 2.0;
r3_lb = 60.0; r3_ub = 90.0;
model_params = store_variables(     # Store the true model parameter values and their names
    ["λ","ζ","R_{d}"], 
    [r1, r2, r3], 
    [r1_lb, r2_lb, r3_lb], 
    [r1_ub, r2_ub, r3_ub], 
    [guess_r1, guess_r2, guess_r3]) 

# Noise model - edit
error_type = "Normal"               # Normal, Lognormal
Sd = 5.0;
guess_Sd = 8.0;
Sd_lb = 1.0; Sd_ub = 10.0;
noise_params = store_variables(     # Store the true noise parameter valuess and their names
    "σ_{N}", 
    Sd, 
    Sd_lb, 
    Sd_ub, 
    guess_Sd) 

# Initial conditions - edit
C10=10.0;
X_array = store_variables(          # Store the initial conditions and the names of each state variable X_i(t)
    "R(t)", 
    C10) 

# Data time-points - edit
t_max = 30.0; 
times = LinRange(0.0,t_max, 20); # synthetic data measurement times at which to generate data points

## Mathematical model: system of differential equations (ODE or PDE) - edit
f(u,p,t) = (p[1]*u/3) * (1 - max(0, (1+p[2]/p[1]) * (u-p[3])^3/(u^3))) # p[1]=λ, p[2]=zeta, and p[3]=R_d
   
#######################################################################################
## Data loading or generation
#######################################################################################

# Generate seed number if not given
if !isdefined(Main, :seed_num)
    seed_num = rand(1000:9999)
end

# Generate data synthetically
Random.seed!(seed_num)
data = makeSyntheticData(times, X_array.initial_values, f, model_params.true_values, error_type, noise_params.true_values)
df_data = DataFrame(data, :auto)  # Use :auto to automatically name columns
CSV.write( "synthetic_data" * experi_name * string(seed_num) * ".csv", df_data)

# Or read a CSV file of data
# data = DataFrame(CSV.File("path/to/your/file.csv"))  # Update with the correct path

#######################################################################################
## Path Wise Analysis
#######################################################################################

pwaFunction(experi_name, times, X_array, f, model_params, error_type, noise_params, seed_num, data, npts)
