include("SolarStrategy.jl")
using .SolarStrategy

using Plots
using Dates
plotly(ticks=:native)

# preparing the track data

track, segments = get_track_and_segments("data/data_australia.csv")
track_peaks_file, segments_peaks_file = get_track_and_segments("data/data_australia_peaks.csv")
plot(
    track.distance / 1000., track.altitude, title="Гоночный маршрут",
    xlabel="Дистанция (км)",
    ylabel="Высота (м)",
    legend=false
)
track_peaks, segments_peaks = keep_extremum_only_peaks_segments(track)
plot(track_peaks.distance, track_peaks.altitude, title="Track extremum only data built w/ Peaks.jl")

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

# preparing the solar car data
# using the SOL data as an example
solar_car = SolarCar(
    390.0,
    0.18,
    1.0,
    0.0023,
    0.000041,
    0.87,
    40.0,
    0.86,
    0.228,
    4.0
)

# and general environment data
env = Environment(
    9.8019,
    1.18
)

# parametric track simplification
k=1.75
_, points_k = parametrized_track_simplification(track, k)
track_k, segments_k = get_track_and_segments_for_selected_points_modified(track, points_k)

track_high_k, segments_high_k = get_track_and_segments_for_selected_points_modified(track_high, points_k)

# this way reduction is very low
k2 = 7.25
_, points_k2 = parametrized_track_simplification(track_high, k2)
track_high_k2, segments_high_k2 = get_track_and_segments_for_selected_points_modified(track_high, points_k2)

# generating weather conditions
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

# high track
segments_peaks_high.weather_coeff = weather_coeff
@time res_peaks_high = iterative_optimization(
    track_peaks_high, segments_peaks_high,
    5,
    5,
    5100.,
    solar_car,
    env,
    DateTime(2023,1,1,10,0,0)
);
plots_for_results(res_peaks_high, track_peaks_high, segments_peaks_high)

simulate_run_finish_time(
    fill(28., size(segments_peaks_high.diff_distance, 1)),
    track_peaks_high,
    segments_peaks_high,
    5100.,
    DateTime(2023,1,1,10,0,0),
    solar_car,
    env
)

# k track high
weather_coeff_k = calculate_weather_weights_for_segments(
    w,
    elat,
    elon,
    segments_k
);
segments_high_k.weather_coeff = weather_coeff_k
@time res_full_high_k = iterative_optimization(
    track_high_k, segments_high_k,
    5,
    5,
    5100.,
    solar_car,
    env,
    DateTime(2023,1,1,10,0,0)
);
plots_for_results(res_full_high_k, track_high_k, segments_high_k)

function k_speeds_to_regular(speeds, points)
	# k_speeds = res_full_high_k[end].solution.speeds

	new_speed_vector_for_regular = []
	for i=2:size(points,1)
		append!(
			new_speed_vector_for_regular,
			fill(
				speeds[i-1],
				points[i] - points[i-1]
			)
		)
	end
	# println("new size $(size(new_speed_vector_for_regular,1))")
	# println("old size $(poi_df.intersect[index]-1)")
	# @assert size(new_speed_vector_for_regular,1)==poi_df.intersect-1
	return new_speed_vector_for_regular
end

speeds_for_regular_track = k_speeds_to_regular(
    res_full_high_k[end].solution.speeds,
    points_k
)

weather_coeff_full = calculate_weather_weights_for_segments(
    w,
    elat,
    elon,
    segments
);
segments_high.weather_coeff = weather_coeff_full
simulate_run_finish_time(
    Float64.(speeds_for_regular_track),
    track_high,
    segments_high,
    5100.,
    DateTime(2023,1,1,10,0,0),
    solar_car,
    env
)


# high k2 track
weather_coeff_k2 = calculate_weather_weights_for_segments(
    w,
    elat,
    elon,
    segments_high_k2
);
segments_high_k2.weather_coeff = weather_coeff_k2
@time res_full_high_k2 = iterative_optimization(
    track_high_k2, segments_high_k2,
    5,
    5,
    5100.,
    solar_car,
    env,
    DateTime(2023,1,1,10,0,0)
);
plots_for_results(res_full_high_k2, track_high_k2, segments_high_k2)

speeds_for_regular_track2 = Float64.(
    k_speeds_to_regular(
        res_full_high_k2[end].solution.speeds,
        points_k2
    )
)

simulate_run_finish_time(
    speeds_for_regular_track2,
    track_high,
    segments_high,
    5100.,
    DateTime(2023,1,1,10,0,0),
    solar_car,
    env
)
