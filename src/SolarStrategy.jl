module SolarStrategy

using DataFrames
# also consider using JuliaDB and Query
using CSV
using Plots # default
using PlotlyBase
# using PlotlySave
import Pluto
using TimeZones
using Dates
using Optim
using LineSearches
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
include("strategy_calculation.jl")


#### concept
#=
overall concept of modelling:

speeds as an input
at first compute the energy loss, since it can be done in vectorized way
energy loss and time on each sector is calculated

now, since we know the times, calculate the energy income
this also can now be done in vectorized way

there are several possible models for energy income

start with stub, develop proper models later
=#

# preparing the track data
track = get_track_data("data/data_australia.csv")

# TODO: track preprocessing
plot(track.distance, track.altitude, title="Track raw data")
# track_test = keep_extremum_only(track)
# plot(track_test.distance, track_test.altitude, title="track extremum only data")
track_test_peaks = keep_extremum_only_peaks(track)
plot(track_test_peaks.distance, track_test_peaks.altitude, title="Track extremum only data built w/ Peaks.jl")

# TODO: slope angle preprocessing

# # few big chunks, LBFGS, 645 seconds, negative speeds, consider revising or constraints, 3886 seconds
# @time result_chunks_lbfgs = optimize(x -> solar_trip_chunks(x, track), [41.0, 42.0, 43.0, 44.0, 45.0], LBFGS())
# minimized_inputs_chunks_lbfgs = Optim.minimizer(result_chunks_lbfgs)
# minimized_result_chunks_lbfgs = Optim.minimum(result_chunks_lbfgs)
# lower = [0.0]
# upper = [100.]
# initial_x = [43.0]
# @time result_chunks_lbfgs_2 = optimize(x -> solar_trip_chunks(x, track), lower, upper, initial_x, Fminbox(LBFGS()))
# minimized_inputs_chunks_lbfgs_2 = Optim.minimizer(result_chunks_lbfgs_2)
# minimized_result_chunks_lbfgs_2 = Optim.minimum(result_chunks_lbfgs_2)
# result = optimize(f, [30.0])
# result = optimize(f, [30.0], LBFGS())
# result = optimize(f, [30.0], GradientDescent())

line_search = LineSearches.BackTracking();
iterations_num = 10000
number_of_chunks = 5
start_speed = 40.


function f(x)
	return solar_trip_chunks(abs.(x), track)
end

td = TwiceDifferentiable(f, fill(start_speed, number_of_chunks), autodiff = :forward)

lower_bound = fill(0.0, number_of_chunks)
upper_bound = fill(100.0, number_of_chunks)

tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)

@time result_chunks = optimize(td, tdc, fill(start_speed, number_of_chunks),
    IPNewton(),
    Optim.Options(
        x_tol = 1e-6,
        f_tol = 1e-8,
        g_tol = 1e-6
    )
)

minimized_inputs_chunks = Optim.minimizer(result_chunks)

inputs_ms = convert_kmh_to_ms(minimized_inputs_chunks)
speed_vector = propagate_speeds(inputs_ms, track)
power_use_chunks, solar_power_chunks, energy_in_system_chunks, time_chunks, time_s_chunks = solar_trip_calculation(speed_vector, track);
last(time_s_chunks)

plot(track.distance, track.altitude, label="altitude", ylabel="altitude", 
    title="Speed (km/h) vs distance", right_margin = 15Plots.mm
)
plot!(twinx(), track.distance, speed_vector * 3.6, color=:red, ylabel="speed (km/h)",
    ylim=[0, 60], label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance"
)


plot(track.distance,
    [power_use_chunks solar_power_chunks energy_in_system_chunks zeros(size(track,1))],
	label=["Energy use" "Energy income" "Energy in system" "Failure threshold"],
    title="Energy graph (distance)",
	xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
	color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)


plot(time_chunks.utc_time,
    [power_use_chunks solar_power_chunks energy_in_system_chunks zeros(size(track,1))],
	label=["Energy use" "Energy income" "Energy in system" "Failure threshold"],
	xlabel="Time", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
	color=[:blue :green :cyan :red]
	# ,ylims=[-10000, 40000]
	, legend = :topleft 
	, right_margin = 15Plots.mm
	, title = "Energy graph (time)"
)
plot!(twinx(), time_chunks.utc_time, speed_vector * 3.6, color=:red, ylabel="speed (km/h)",
    ylim=[0, 60], label="speed (km/h)", ymirror = true,
    title = "Energy graph (time)"
)

### HIERARCHICAL

track_size = 100
# size(keep_extremum_only_peaks(track),1)
# short_track = keep_extremum_only_peaks(track)[1:track_size, :]
short_track = keep_extremum_only_peaks(track)
track_size = size(short_track.distance, 1)
distance_perc = last(short_track).distance / last(track.distance)
proposed_start_energy = distance_perc * 5100
start_energy_short = proposed_start_energy
initial_speed = 40.0
chunks_amount_hierarchical = 10
start_datetime_hierarchical = DateTime(2022, 7, 1, 0, 0, 0)


@time result_hierarchical = hierarchical_optimization(
    initial_speed, short_track, chunks_amount_hierarchical,
    start_energy_short, 0., start_datetime_hierarchical, 1, track_size
)

inputs_ms_hier = convert_kmh_to_ms(result_hierarchical)
power_use_hier, solar_power_hier, energy_in_system_hier, time_hier, time_s_hier = solar_trip_calculation(inputs_ms_hier, short_track, start_energy_short)
println(last(time_s_hier))

plot(
    short_track.distance, short_track.altitude,
    label="altitude", ylabel="altitude",
    title="Speed (km/h) vs distance", right_margin = 15Plots.mm,
    size=(1200, 500)
    )
plot!(
    twinx(), short_track.distance, result_hierarchical,
    color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true,
    title="Speed (km/h) vs distance",
    size=(1200, 500)
    )

plot(short_track.distance, [power_use_hier solar_power_hier energy_in_system_hier zeros(track_size)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Hierarchical",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)

# TODO: penalty for huge speed change?
# OR see why that happens

# TODO: run line-by-line in debug


# TODO: refactor code, calculations in separate function, wrapper for optimization
# make an optimization with full vector
# make different wrappers for calculation: only cost and all data

#### future use
# for optimization (overall list: https://www.juliaopt.org/packages/ ):
# https://github.com/robertfeldt/BlackBoxOptim.jl - looks like the case
# https://github.com/JuliaNLSolvers/Optim.jl - not sure if it handles derivative-free optimi
# https://github.com/JuliaOpt/Pajarito.jl - integer linear programming
# https://github.com/anriseth/MultiJuMP.jl - multicriterial optimiztion
# also http://julia.cookbook.tips/doku.php?id=optim , Nelder-Mead Simplex	NM
# https://github.com/SciML/Optimization.jl - wraper for optimizers
# https://github.com/JuliaFirstOrder/ProximalOperators.jl - proximal operators for semi-continuous functions

# https://juliahub.com/ui/Packages/LightGraphs/Xm08G/1.3.5?t=0 for graphs


end # module
