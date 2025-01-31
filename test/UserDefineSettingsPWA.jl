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
        "StatsPlots", "StatsBase", "Colors", "Suppressor"]

# Activate the environment with the appropriate packages via the function with inputs:
# (activOption, envPath, requiredPackages)
envActivate(1, env_path, required_packages)

# Use the package
using LogLikelihoodBasedProfileWiseAnalysis

#######################################################################################
## User defined settings
#######################################################################################
## Seed the simulation if desired - edit
seed_num = 1234

## Number of points between guess and bounds for loglikelihood exploration - edit
npts = 60;

## Parameters
# Create the struct for information storage
store_variables = defineStruct()
# Store parameter and variable information in the following form:
#       store_variables([name1, name 2, ...], [true_value1, true_value2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [guess1, guess2, ...])
# Or if true values are unknown:
#       store_variables([name1, name 2, ...], [lower_bound1, lower_bound2, ...], [upper_bound1, upper_bound2,...], [guess1, guess2, ...])
# Or if true values and bounds are unnecessary:
#       store_variables([name1, name 2, ...], [guess1, guess2, ...])
# Or if true values, bounds and guesses are unnecessary:
#       store_variables([name1, name 2, ...])

# Mechanistic model - edit
#r1=1.0;
#r2=0.5;
guess_r1 = 0.15;
guess_r2 = 0.025;
guess_r3 = 200
r1_lb = 0.01; r1_ub = 0.2;
r2_lb = 0.001; r2_ub = 0.5;
r3_lb = 180; r3_ub = 225;
model_params = store_variables(     # Store the true model parameter values and their names
    ["λ","ζ", "R_d"],                #   Labels
    #[r1, r2],                      #   True values if known
    [r1_lb, r2_lb, r3_lb],                 #   Lower bounds to guesses
    [r1_ub, r2_ub, r3_ub],                 #   Upper bounds to guesses
    [guess_r1,guess_r2,guess_r3])            #   Initial guesses

# Noise model - edit
error_type = "Normal"               # Normal, Lognormal
#Sd = 5;
guess_Sd = 10                        # Lognormal: 0.05                        # Normal: 8;
Sd_lb = 8; Sd_ub = 12;          # Lognoamral: Sd_lb = 0.01; Sd_ub = 0.1; # Normal: Sd_lb = 1.0; Sd_ub = 15.0;
noise_params = store_variables(     # Store the true noise parameter valuess and their names
    "σ_{N}",                        #   Label
    #Sd,                            #   True value if known
    Sd_lb,                          #   Lower bound to guesses
    Sd_ub,                          #   Upper bound to guesses
    guess_Sd)                       #   Initial guess

# State variables - edit
X_variables = store_variables(      # Store the names of each state variable X_i(t)
    ["R(t)", "R_n(t)"])                       #   Label - even if only one state variable, must be in an array []

# Initial conditions - edit
#R0=140.0;
guess_R0 = 175
guess_Rn0 = 0
R0_lb = 160; R0_ub = 190
Rn0_lb = 0; Rn0_ub = 10
X_init_condits = store_variables(   # Store the initial conditions of each state variable X_0 and their names
    ["R(0)", "R_n(0)"],                         #   Label
    #R0,                             #   True value if known
    [R0_lb, Rn0_lb],                          #   Lower bound to guesses
    [R0_ub, Rn0_ub],                         #   Upper bound to guesses
    [guess_R0, guess_Rn0])                       #   Initial guess

## Mathematical model: system of differential equations (ODE or PDE) - edit # NEED TO ADD THE CAPABILITY TO WORK WITH PDEs TO THE WORKFLOW
# 1 state variable example:
mechanistic_model_name = "RadialDeathOuterNecrotic"
# dif_equ(u,p,t) = (p[1] .* u ./ 3) .* (1 .- max.(0, (1 + p[2] / p[1]) .* (u .- p[3]) .^3 ./(u .^3))) # Radial death (outer only): p[1]=λ, p[2]=ζ, and p[3]=R_d
# dif_equ(u,p,t) = (p[1].*u/3) .* (1 .- u/p[2]) # Logistic growth: p[1]=λ and p[2]=R_max

# More than 1 state variable example:
function dif_equ!(du,u,p,t) # ! is "bang" convention indicating that the function modifies its inputs
    λ,ζ,R_d=p
    du[1]=(λ .* u[1] ./ 3) .* (1 .- max.(0, (1 + ζ / λ) .* (u[1] .- R_d) .^3 ./(u[1] .^3)));
    if u[1] >= R_d          # R_d is the radial width between the outer radius and the necrotic core radius
        du[2]=du[1];        # If R(t)>R_d, the necrotic core begins growing
    else
        du[2]=0;            # Otherwise the necrotic core doesn't exist
    end
end

## Specifications required if ... - edit
# ... loading data from file - edit
desired_time_points = [1 3 8 12 17]; # Data measurement times at which to use data points (Resolution A, Murphy et al. 2022)

# OUTER RADIUS
#file_name_all = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth dynamics/Code/Murphy2023ErrorModels-main/WM793b_IncuCyte_confocal_outer.csv"  # Raw data
file_name_all_outer = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth dynamics/Code/Murphy2023ErrorModels-main/MurphyData_logisticGrowthModel_NormalErr_initCellNum-5000_RNG1234_tumourGrowthData.csv"  # Tidy data for tumour growth
file_name_specific_outer = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth dynamics/Code/Murphy2023ErrorModels-main/MurphyData_logisticGrowthModel_NormalErr_initCellNum-5000_RNG1234_tumourGrowthData_specificTimePoints_rawNumMsrmntsPerTp.csv" # Tidy data with extracted data at specific time-points; no restrictions on number of measurements
#file_name_specific_outer = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth dynamics/Code/Murphy2023ErrorModels-main/MurphyData_logisticGrowthModel_NormalErr_initCellNum-5000_RNG1234_tumourGrowthData_specificTimePoints_restrictNumMsrmntsPerTp.csv" # Tidy data with extracted data at specific time-points; restrictions on number of measurements

# NECROTIC RADIUS
file_name_all_necrotic = "/Users/aliladearie/Documents/RA_work/Heterogeneity in tumour spheroid growth dynamics/Code/Murphy2023ErrorModels-main/WM793b_confocal_necrotic.csv"  # Raw data

day_growth_begins = 4           # The number of days after tumour formation is initialised determined as the first day of tumour growth
restrict_measures_per_tp = 1    # Determines if the data will be extracted to have the same number of measurements per timepoint:
                                #   0: the data will be extracted at the desired time-points with no restrictions
                                #   1: the number of measurements at each time-point will be restricted by the minimum amount of measurements made for any one state variable at any one time-point
                                #   If this parameter isn't defined, it is assumed that the loaded data will be used without adjustment
                                #   COMMENT OUT THIS PARAMETER IF YOU WISH TO USE THE DATA AS IS, AFTER TIDYING AND TIME-POINT EXTRACTION IF RELEVANT
restricted_num = 50             # Attempt to enforce this as the restricted number of measurements per time-point (will be overruled if one particular time-point had fewer measurements)
initial_cell_number = 5000      # The number of cells with which the tumour is initialised to grow
type_of_data =                  # Determines how data will be loaded per state variable:
    ["extracted_tps", "raw"]    #   "raw": data is in original format and needs to be extracted for the desired cell seeding number, and the desired time-points
                                #   "tidy": data_all and times_all have been extracted for the desired cell seeding number, but need to be extracted for desired timepoints
                                #   "extracted_tps": data_all, times_all, data_specific and times_specific have all been extracted for the appropriate time-points
variables_included =            # The measurements taken in the dataset
    ["outerR", "necroticR"]  

    # ... enerating data synthetically
#t_max = 2.0;                       # Data time-points maximum



#######################################################################################
## Data loading or generation
#######################################################################################

# Define experiment name
core_name = "MurphyData_" * mechanistic_model_name * "Model_" * error_type * "Err_initCellNum-" * string(initial_cell_number) * "_RNG" * string(seed_num)

if @isdefined(restrict_measures_per_tp)
    if restrict_measures_per_tp == 1
        experi_name = core_name * "_" * join(type_of_data,"-") * "_restrictNumMsrmntsPerTp"
    elseif restrict_measures_per_tp == 0
        experi_name = core_name * "_" * join(type_of_data,"-") * "_rawNumMsrmntsPerTp"
    end
else
    experi_name = core_name * "_" * join(type_of_data,"-") * "_asIsPerTp"
end

# Generate seed number if not given
if !isdefined(Main, :seed_num)
    seed_num = rand(1000:9999)
end

# Set the random seed
Random.seed!(seed_num)

# Extract the data/times at all valid time-points (_all), and at the desired time-points (_specific)
# Outer radius

data_all_necrotic, times_all_necrotic, data_specific_necrotic, times_specific_necrotic, restricted_num =
    loadExtractData(file_name_all_necrotic,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data[2],core_name,variables_included[2],seed_num,restricted_num)
#    loadExtractData(file_name_all,file_name_specific,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data[2],core_name,variables_included[2],seed_num,restricted_num)# Only for "extracted_tps" data-type (loading tidy data that has been extracted at time points)
#    loadExtractData(file_name_all,file_name_specific,initial_cell_number,day_growth_begins,desired_time_points,type_of_data[2],core_name,variables_included[2],seed_num,restricted_num)                          # Only for "extracted_tps" data-type (loading tidy data that has been extracted at time points); appropriate for data that you wish to use as is
println("Smallest number of outer radius measurements per timepoint: " * string(restricted_num))

data_all_outer, times_all_outer, data_specific_outer, times_specific_outer, restricted_num =
#    loadExtractData(file_name_all_outer,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data[1],core_name,variables_included[1],seed_num,restricted_num)
    loadExtractData(file_name_all_outer,file_name_specific_outer,restrict_measures_per_tp,initial_cell_number,day_growth_begins,desired_time_points,type_of_data[1],core_name,variables_included[1],seed_num,restricted_num) # Only for "extracted_tps" data-type (loading tidy data that has been extracted at time points)
#    loadExtractData(file_name_all_outer,file_name_specific_outer,initial_cell_number,day_growth_begins,desired_time_points,type_of_data[1],core_name,variables_included[1],seed_num,restricted_num)                          # Only for "extracted_tps" data-type (loading tidy data that has been extracted at time points); appropriate for data that you wish to use as is
println("Smallest number of necrotic radius measurements per timepoint: " * string(restricted_num))


# Store the data/times in the appropriate format if more than one state variables
data_all = storeMultiVars(data_all_outer,data_all_necrotic)
times_all = storeMultiVars(times_all_outer,times_all_necrotic);
data_specific = storeMultiVars(data_specific_outer,data_specific_necrotic)
times_specific = storeMultiVars(times_specific_outer,times_specific_necrotic)

#=
# Generate data synthetically - edit
t_max = 2.0;                       # Data time-points maximum
times_specific = LinRange(0.0,t_max, 16); # synthetic data measurement times at which to generate data points
noNegatives = true                # Is synthetic data allowed to be negative? true/false
data = makeSyntheticData(times_specific, X_array.initial_values, DE!, model_params.true_values, error_type, noise_params.true_values, noNegatives)
df_data = DataFrame(data, :auto)    # Use :auto to automatically name columns
CSV.write( "synthetic_data" * experi_name * string(seed_num) * ".csv", df_data)
=#

#######################################################################################
## Check for potential errors
#######################################################################################
if !isnothing(model_params.guess)
    for i in eachindex(model_params.guess)
        if model_params.guess[i]<model_params.lower_bound[i] || model_params.guess[i]>model_params.upper_bound[i]
            error("Initial guesses for " * model_params.name[i] * " are not within the bounds given.\n")
        end
    end
end
if !isnothing(noise_params.guess)
    for i in eachindex(noise_params.guess)
        if noise_params.guess[i]<noise_params.lower_bound[i] || noise_params.guess[i]>noise_params.upper_bound[i]
            error("Initial guesses for " * noise_params.name[i] * " are not within the bounds given.\n")
        end
    end
end
#######################################################################################
## Path Wise Analysis
#######################################################################################

pwaFunction(experi_name, X_variables, X_init_condits, model_params, error_type, noise_params, dif_equ!, times_all, times_specific, data_all, data_specific, npts, seed_num)
#pwaFunction(experi_name, X_variables, X_init_condits, model_params, error_type, noise_params, dif_equ, times_all_outer, times_specific_outer, data_all_outer, data_specific_outer, npts, seed_num)