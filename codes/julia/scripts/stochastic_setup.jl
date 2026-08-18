# rm(sourceDir*"Manifest.toml", force=true)

# Recreate environment
using Pkg
Pkg.activate(sourceDir)
# Pkg.instantiate()
# Pkg.resolve()


# using Plots
using DataFrames
using CSV
using Statistics
using Random
using ProgressMeter

# Disable plot display, only save
# ENV["GKSwstype"] = "nul"  # For headless operation
# gr(show=false)  # Don't display plots

# Ensure directories exist

threads = Threads.nthreads()
# sourceDir = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/"
include(joinpath(sourceDir,"src/RACIPEdata.jl"))
using .RACIPEdata

include(joinpath(sourceDir,"src/RACIPErunner.jl"))
using .RACIPErunner

include(joinpath(sourceDir,"src/ODEgenerator.jl"))
using .ODEgenerator


include(joinpath(sourceDir, "src/StochasticSimulations.jl"))
using .StochasticSimulations
