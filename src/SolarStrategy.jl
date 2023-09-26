module SolarStrategy

using DataFrames
# also consider using JuliaDB and Query
using CSV
using Plots # default
using PlotlyBase
# using PlotlySave
using TimeZones
using Dates
using Optim
using LineSearches
using ProtoStructs
using OhMyREPL
using Distributions
using Random
using StatsBase
# using FLoops
using Distributed
using ProgressMeter
using BenchmarkTools
using JSON

# using Revise
# includet("energy_draw.jl")
# includet("time.jl")
# includet("solar_radiation.jl")
# includet("track.jl")
# includet("utils.jl")
# includet("weather.jl")
# includet("strategy_calculation.jl")
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
# track = get_track_data("data/data_australia.csv")
track, segments = get_track_and_segments("data/data_australia.csv")
plot(
    track.distance / 1000., track.altitude, title="Гоночный маршрут",
    xlabel="Дистанция (км)",
    ylabel="Высота (м)",
    legend=false
)
track_peaks, segments_peaks = keep_extremum_only_peaks_segments(track)
plot(track_peaks.distance, track_peaks.altitude, title="Track extremum only data built w/ Peaks.jl")
# TODO: track preprocessing
track_high = copy(track)
track_high.altitude = track_high.altitude .* 10;
segments_high = get_segments_for_track(track_high);
track_peaks_high, segments_peaks_high = keep_extremum_only_peaks_segments(track_high)
plot(track_peaks_high.distance, track_peaks_high.altitude, title="Track extremum only data built w/ Peaks.jl alt * 10")

plot(
    track_peaks.distance / 1000., track_peaks.altitude, title="Сокращённый гоночный маршрут",
    xlabel="Дистанция (км)",
    ylabel="Высота (м)",
    legend=false
)


from_index = 1
to_index = 3000

from_index = 3000
to_index = 7000

from_index = 7000
to_index = size(track_peaks.distance, 1)

plot(
    track_peaks.distance[from_index:to_index] / 1000., track_peaks.altitude[from_index:to_index],
    title="Сокращённый гоночный маршрут, подзадача 3",
    xlabel="Дистанция (км)",
    ylabel="Высота (м)",
    legend=false
)

track_aus, segments_aus = get_track_and_segments("data/data_australia_random.csv");
track_aus.altitude = track_aus.altitude * 10;
segments_aus = get_segments_for_track(track_aus);
plot(track_aus.distance, track_aus.altitude, title="Track aus short")


dimensions=15;
w, elat, elon = generate_clouds(
	-10,
	130,
	-18,
	135,
	-12.5,
	131.2,
	0.5,
	0.5,
	dimensions,
	0.75
);

weather_coeff = calculate_weather_weights_for_segments(
    w,
    elat,
    elon,
    segments_peaks
);
segments_peaks.weather_coeff = weather_coeff

# @time res = iterative_optimization(track, segments, 5, 5100., DateTime(2022,7,1,0,0,0));
# 41 sec for two iterations
# regular track
@time res_peaks = iterative_optimization(
    track_peaks, segments_peaks,
    5,
    5,
    5100.,
    DateTime(2023,1,1,10,0,0)
);
plots_for_results(res_peaks, track_peaks, segments_peaks)

simulate_run_finish_time(
    fill(37.90363171504429, size(segments_peaks.diff_distance, 1)),
    track_peaks,
    segments_peaks,
    5100.,
    DateTime(2023,1,1,10,0,0)
)

# high track
segments_peaks_high.weather_coeff = weather_coeff
@time res_peaks_high = iterative_optimization(
    track_peaks_high, segments_peaks_high,
    5,
    5,
    5100.,
    DateTime(2023,1,1,10,0,0)
);
plots_for_results(res_peaks_high, track_peaks_high, segments_peaks_high)

@time single_speed_high = minimize_single_speed(
    track_peaks_high,
    segments_peaks_high,
    5100.,
    DateTime(2023,1,1,10,0,0),
    31.
)

simulate_run_finish_time(
    fill(28., size(segments_peaks_high.diff_distance, 1)),
    track_peaks_high,
    segments_peaks_high,
    5100.,
    DateTime(2023,1,1,10,0,0)
)

# full track
weather_coeff_full = calculate_weather_weights_for_segments(
    w,
    elat,
    elon,
    segments
);
segments.weather_coeff = weather_coeff_full
@time res_full = iterative_optimization(
    track, segments,
    5,
    5,
    5100.,
    DateTime(2023,1,1,10,0,0)
);
plots_for_results(res_full, track, segments)

simulate_run_finish_time(
    fill(37.90363171504429, size(segments.diff_distance, 1)),
    track,
    segments,
    5100.,
    DateTime(2023,1,1,10,0,0)
)

# full high track
segments_high.weather_coeff = weather_coeff_full
@time res_full_high = iterative_optimization(
    track_high, segments_high,
    5,
    5,
    5100.,
    DateTime(2023,1,1,10,0,0)
);
plots_for_results(res_full_high, track_high, segments_high)



@time res_aus = iterative_optimization(
    track_aus, segments_aus,
    5,
    5,
    5100.,
    DateTime(2023,1,1,0,0,0)
);
plots_for_results(res_aus, track_aus, segments_aus)



# 13 secs for 2 iterations
# 126 seconds full

# iter_prev=res[end-1];
# iter_last=res[end];
# subtasks_prev = iter_prev.subtasks;
# subtasks_last = iter_last.subtasks;

# segments_size_prev = map((x) -> x.subtask_boundaries.size, subtasks_prev);
# segments_size_last = map((x) -> x.subtask_boundaries.size, subtasks_last);

# linetypes_test = [ :stepmid]

# plot(segments_size_prev, line=linetypes_test)
# plot(segments_size_last, line=linetypes_test)

# # graphs looks ok, we can continue to optimization





single_speed = minimize_single_speed(
    track_peaks,
    segments_peaks,
    5100.,
    DateTime(2023,1,1,10,0,0),
    31.
)

simulate_run_finish_time(
    fill(first(single_speed), size(segments_peaks.diff_distance, 1)),
    track_peaks,
    segments_peaks,
    5100.,
    DateTime(2023,1,1,10,0,0)
)





n_speeds = minimize_n_speeds(
    track_peaks,
    segments_peaks,
    2,
    5100.,
    DateTime(2023,1,1,10,0,0),
    first(single_speed)
    # 35.
)

boundaries = calculate_boundaries(1, size(track_peaks, 1), 2)

out_speeds = set_speeds_boundaries(n_speeds, boundaries)

simulate_run(
    out_speeds,
    track_peaks,
    segments_peaks,
    5100.,
    DateTime(2023,1,1,10,0,0)
)


@time n_speeds_high = minimize_n_speeds(
    track_peaks_high,
    segments_peaks_high,
    11,
    5100.,
    DateTime(2023,1,1,10,0,0),
    first(single_speed_high)
    # 35.
)

boundaries = calculate_boundaries(1, size(track_peaks_high, 1), 11)

out_speeds_high = set_speeds_boundaries(n_speeds_high, boundaries)

simulate_run_finish_time(
    out_speeds_high,
    track_peaks_high,
    segments_peaks_high,
    5100.,
    DateTime(2023,1,1,10,0,0)
)

# TODO: somehow ensure that we are NOT running out of energy


############################
# making an optimization without hierarchy/iterations

variables_num = 100
variable_boundaries = calculate_boundaries(
    1,
    size(track_peaks, 1),
    variables_num
)

function f_wrap(input_speeds)
    return solar_partial_trip_wrapper_iter(
        input_speeds, segments_peaks, variable_boundaries,
        5100., 0.,
        DateTime(2022,7,1,0,0,0)
    )
end

init_speeds = fill(44., variables_num)

td = TwiceDifferentiable(f_wrap, init_speeds; autodiff = :forward)
lower_bound = fill(0.0, variables_num)
upper_bound = fill(100.0, variables_num)
tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
# line_search = LineSearches.BackTracking();
# result = optimize(td, fill(speed, vars_amount),
    #Newton(; linesearch = line_search),
@time result = optimize(td, tdc, init_speeds 
# .+ rand(vars_amount) .- 0.5
    ,
    IPNewton(),
    Optim.Options(
        x_tol = 1e-12,
        f_tol = 1e-12,
        g_tol = 1e-12
    )
)
minimized_speeds = Optim.minimizer(result)
# TODO: calculate distance at split points
# calculcate time
# calculate energy

plot(
    minimized_speeds,
    line=:stepmid,
    title=string(
        "Single optim results "
    ) 
)

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


# @time result_hierarchical = hierarchical_optimization(
#     initial_speed, short_track, chunks_amount_hierarchical,
#     start_energy_short, 0., start_datetime_hierarchical, 1, track_size
# )

@time result_hierarchical = hierarchical_optimization_alloc(
    initial_speed, short_track, chunks_amount_hierarchical,
    start_energy_short, 0., start_datetime_hierarchical, 1, track_size
)

iter_speeds = []
push!(iter_speeds, fill(initial_speed, size(short_track.distance, 1)))

@time result_hierarchical = hierarchical_optimization_alloc!(iter_speeds,
    initial_speed, short_track, chunks_amount_hierarchical,
    start_energy_short, 0., start_datetime_hierarchical, 2, track_size
)

# for i=1:size(iter_speeds,1)
#     plot(short_track.distance, iter_speeds[i], linetype = :steppost,
#     title="Speed (km/h) on iteration $(i)")
# end

plot(short_track.distance, iter_speeds[1], linetype = :steppost,
    title="Speed (km/h) on iteration 1")
plot(short_track.distance, iter_speeds[2], linetype = :steppost,
    title="Speed (km/h) on iteration $(2)")
plot(short_track.distance, iter_speeds[3], linetype = :steppost,
    title="Speed (km/h) on iteration $(3)")
plot(short_track.distance, iter_speeds[4], linetype = :steppost,
    title="Speed (km/h) on iteration $(4)")
plot(short_track.distance, iter_speeds[5], linetype = :steppost,
    title="Speed (km/h) on iteration 5")
plot(short_track.distance, iter_speeds[6], linetype = :steppost,
    title="Speed (km/h) on iteration 5")

inputs_ms_hier = convert_kmh_to_ms(result_hierarchical)
power_use_hier, solar_power_hier, energy_in_system_hier, time_hier, time_s_hier = solar_trip_calculation(inputs_ms_hier, short_track, start_energy_short)
println(last(time_s_hier)) #247523.356

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
    size=(1200, 500),
    linetype = :steppost
    )

linetypes = [:path :steppost]
plot(
    short_track.distance, [short_track.altitude, result_hierarchical],
    layout=(2,1),
    labels=["altitude", "speed (km/h)"], ylabel=["altitude","2"],
    line=(linetypes,2), lab=map(string, linetypes),
    # lines=linetypes,
    title="Speed (km/h) vs distance", right_margin = 15Plots.mm,
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

# also see:
# https://github.com/JuliaSmoothOptimizers/ADNLPModels.jl - from https://github.com/JuliaSmoothOptimizers 
# https://jump.dev - for constrained problems
# https://github.com/JuliaOpt/NLopt.jl - an interface to NLopt, later grew to
# https://github.com/SciML/Optimization.jl - invoke NLopt from here
end # module
