module SolarStrategy

using DataFrames # also consider using JuliaDB and Query
using CSV
using Plots # default
using PlotlyBase
using Dates
using Optim
using Distributions
using Random
using StatsBase
using Distributed
using ProgressMeter
using JSON

# dev dependencies, not used as of right now
# using Revise
# using Alert
# selecting a Plots backend
plotly(ticks=:native)
# consider using Gadfly http://gadflyjl.org/stable/
# using Plotly # not as fast, but interactive
# using PlotlyJS # a lot of dependencies, slow loading

include("energy_draw.jl")
include("time.jl")
include("solar_radiation.jl")
include("track.jl")
include("utils.jl")
include("weather.jl")
include("strategy_calculation.jl")

export get_track_and_segments, keep_extremum_only_peaks_segments, get_segments_for_track
export generate_clouds, calculate_weather_weights_for_segments, read_weather_json, write_weather_json
export generate_density_data
export iterative_optimization
export plots_for_results
export simulate_run_finish_time
export get_mean_data
export parametrized_track_simplification, get_track_and_segments_for_selected_points_modified
export SolarCar, Environment
# #### concept
# #=
# overall concept of modelling:

# speeds as an input
# at first compute the energy loss, since it can be done in vectorized way
# energy loss and time on each sector is calculated

# now, since we know the times, calculate the energy income
# this also can now be done in vectorized way

# there are several possible models for energy income

# start with stub, develop proper models later
# =#

# #### future use
# # for optimization (overall list: https://www.juliaopt.org/packages/ ):
# # https://github.com/robertfeldt/BlackBoxOptim.jl - looks like the case
# # https://github.com/JuliaNLSolvers/Optim.jl - not sure if it handles derivative-free optimi
# # https://github.com/JuliaOpt/Pajarito.jl - integer linear programming
# # https://github.com/anriseth/MultiJuMP.jl - multicriterial optimiztion
# # also http://julia.cookbook.tips/doku.php?id=optim , Nelder-Mead Simplex	NM
# # https://github.com/SciML/Optimization.jl - wraper for optimizers
# # https://github.com/JuliaFirstOrder/ProximalOperators.jl - proximal operators for semi-continuous functions

# # https://juliahub.com/ui/Packages/LightGraphs/Xm08G/1.3.5?t=0 for graphs

# # also see:
# # https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl - from https://github.com/JuliaSmoothOptimizers 
# # https://jump.dev - for constrained problems
# # https://github.com/JuliaOpt/NLopt.jl - an interface to NLopt, later grew to
# # https://github.com/SciML/Optimization.jl - invoke NLopt from here

end # module
