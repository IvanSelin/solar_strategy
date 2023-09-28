# using DataFrames
# # also consider using JuliaDB and Query
# using CSV
# using Plots # default
# using PlotlyBase
# # using PlotlySave
# using TimeZones
# using Dates
# using Optim
# using LineSearches
# using ProtoStructs
# using OhMyREPL
# using Distributions
# using Random
# using StatsBase
# # using FLoops
# using Distributed
# using ProgressMeter
# using BenchmarkTools
# using JSON

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
include("SolarStrategy.jl")
using .SolarStrategy

using Plots
using Dates
# using CSV
plotly(ticks=:native)


# include("energy_draw.jl")
# include("time.jl")
# include("solar_radiation.jl")
# include("track.jl")
# include("utils.jl")
# include("weather.jl")
# include("strategy_calculation.jl")


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
track_peaks_file, segments_peaks_file = get_track_and_segments("data/data_australia_peaks.csv")
# peaks_temp = copy(track_peaks)
# peaks_temp.distance = peaks_temp.distance / 1000.
# CSV.write("data/data_australia_peaks.csv", peaks_temp)
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


# from_index = 1
# to_index = 3000

# from_index = 3000
# to_index = 7000

# from_index = 7000
# to_index = size(track_peaks.distance, 1)

# plot(
#     track_peaks.distance[from_index:to_index] / 1000., track_peaks.altitude[from_index:to_index],
#     title="Сокращённый гоночный маршрут, подзадача 3",
#     xlabel="Дистанция (км)",
#     ylabel="Высота (м)",
#     legend=false
# )

# track_aus, segments_aus = get_track_and_segments("data/data_australia_random.csv");
# track_aus.altitude = track_aus.altitude * 10;
# segments_aus = get_segments_for_track(track_aus);
# plot(track_aus.distance, track_aus.altitude, title="Track aus short")


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

w2, elat2, elon2 = generate_clouds(
    -10,
    130,
    -35,
    140,
    -25,
    133,
    0.5,
    0.5,
    dimensions,
    0.75
)
write_weather_json(w2,elat2,elon2,"weather_coeff_full_2")

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
