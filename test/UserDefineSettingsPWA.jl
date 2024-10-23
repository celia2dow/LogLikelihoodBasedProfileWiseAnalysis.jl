# Dowling (2024)
# Adapted from Murphy et al. (2023) 

# Calls MakeSyntheticData.jl to generate synthetic data and calls 
# PWAFunction.jl given a user-defined error and mechanistic model,
# and user-defined guesses and bounds for the parameters.

#module LogLikelihoodBasedProfileWiseAnalysis
#export MakeSyntheticData
#export PWAfunction

include("Functions.jl")
include("PWAFunction.jl")

#######################################################################################
## USER DEFINED SETTINGS
#######################################################################################

## Parameters
# Model and noise parameters
r1=1.0;
r2=0.5;
Sd = 0.4;
model_params = [r1, r2]
noise_params = Sd
# Initial conditions
C10=100.0;
C20=10.0;
ICs = [C10, C20]
# Data time-points
t_max = 5.0;
times = LinRange(0.0,t_max, 20); # synthetic data measurement time

## Mathematical model: system of differential equations (ODE or PDE)
# 1 state variable
# f(u,p,t) = (p[1]*u/3) * (1 - max(0, (1+p[2]/p[1]) * (u-p[3])^3/(u^3))) # p[1]=λ, p[2]=zeta, and p[3]=R_d
# More than 1 state variable
function DE!(du,u,p,t)  # ! is "bang" convention indicating that the function modifies its inputs
    r1,r2=p
    du[1]=-r1*u[1];
    du[2]=r1*u[1] - r2*u[2];
end

## Data: can be authentic or synthetically generated
MakeSyntheticData(model_params, noise_params, ICs, times, DE!)
    



