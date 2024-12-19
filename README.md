# LogLikelihoodBasedProfileWiseAnalysis

[![Build Status](https://github.com/celia2dow/LogLikelihoodBasedProfileWiseAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/celia2dow/LogLikelihoodBasedProfileWiseAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package performs profile-wise analysis on a dataset given a provided mathematical ordinary differential equation model and an accompanying error model. 

The dataset may be:
- generated synthetically, requiring true parameter values to be provided, or
- provided in .csv format where each row of data corresponds to a different state variable

All negative values generated in a synthetic dataset are set to 0 - need to update `makeSyntheticData.jl` if negative values are desired.
