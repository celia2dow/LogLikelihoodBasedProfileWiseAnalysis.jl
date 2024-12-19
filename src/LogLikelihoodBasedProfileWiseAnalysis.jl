## Create a module defining all of the packages and functions required for this package
module LogLikelihoodBasedProfileWiseAnalysis

## Install package dependencies
using Plots
using LinearAlgebra
using NLopt
using Interpolations
using Distributions
using Roots
using LaTeXStrings
using DataFrames
using CSV
using DifferentialEquations
using Random
using StatsPlots
using Colors
using Measures

## Export functions to be used
# Plots
export plot
export plot!
export scatter!
export savefig
# Interpolations
export LinearInterpolation
# Distributions
export Normal
export LogNormal
# LaTeXStrings
export @L_str
# DataFrames
export DataFrame            
# CSV
export CSV
# Random
export Random
# Using StatsPlots
export quantile
export qqplot
# Colors
export palette
export RGB
# Measures               
export mm
# makeSyntheticData
export makeSyntheticData    # Function for generating synthetic data from a given model and paramters
# pwaFunction
export pwaFunction          # Function for running the profile-wise analysis and generating summary figures
# odesolver
export odesolver            # Function for solving an ODE system with a given set of paramters, initial conditions, and times
# optimisingFunctions
export optimise             # Function for finding the paramter values θ that maximise the loglikelihood of the data 
export sum_loglikelihood    # Function for finding the sum of loglikelihood values
# funcInterpCI
export funcInterpCI         # Function for interpolating a confidence interval
# defineStruct
export defineStruct         # Function for defining the struct to store parameter information

## Include the files containing the functions
include("makeSyntheticData.jl")
include("pwaFunction.jl")
include("odesolver.jl")
include("optimisingFunctions.jl")
include("funcInterpCI.jl")
include("defineStruct.jl")

end