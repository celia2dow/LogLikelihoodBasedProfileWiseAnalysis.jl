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
using CSV
using DataFrames
using DifferentialEquations
using Random
using StatsPlots
using StructuralIdentifiability

## Export functions to be used
export makeSyntheticData
export pwaFunction
export odesolver

## Include the files containing the functions
include("makeSyntheticData.jl")
include("pwaFunction.jl")
include("odesolver.jl")

end