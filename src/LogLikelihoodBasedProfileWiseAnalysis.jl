## Create a module defining all of the packages and functions required for this package
module LogLikelihoodBasedProfileWiseAnalysis

## Install package dependencies
using Plots
using LinearAlgebra
using NLopt
using .Threads 
using Interpolations
using Distributions
using Roots
using LaTeXStrings
using DataFrames
using CSV
using DifferentialEquations
using Random
using StatsPlots
using StructuralIdentifiability
using Colors

## Export functions to be used
export makeSyntheticData
export pwaFunction
export odesolver
export DataFrame
export CSV
export sum_loglikelihood
export optimise

## Include the files containing the functions
include("makeSyntheticData.jl")
include("pwaFunction.jl")
include("odesolver.jl")
inxlusw("optimisingFunctions.jl")

end