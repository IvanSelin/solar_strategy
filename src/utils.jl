function convert_kmh_to_ms(speed_kmh)
    return speed_kmh / 3.6;
end

function convert_kmh_to_ms_typed(speed_kmh :: Vector{<: Real}) :: Vector{<: Real}
    return speed_kmh ./ 3.6;
end

function convert_kmh_to_ms!(speed_kmh :: Vector{<: Real})
    return speed_kmh .= speed_kmh ./ 3.6;
end

function plots_for_results(res, track, segments)
    for r in res
        display(
            plot_resulted_plan(
                track,
                segments,
                r.solution.speeds,
                r.solution.energies,
                r.solution.times,
                "iteration $(r.number)"
            )
        )
    end
end

function simulate_run_finish_time(
	speeds :: Vector{<: Real}, track :: DataFrame, segments :: DataFrame, start_energy :: Real, start_datetime :: DateTime, solar_car :: SolarCar, env :: Environment)
	minimized_speeds_ms = speeds / 3.6
	
	power_use, solar_power, time_s = solar_trip_boundaries_typed(
		minimized_speeds_ms, segments, start_datetime, solar_car, env
	)
	# track points, not segments, that's why it is size is +1 
	energy_in_system_new = start_energy .+ solar_power .- power_use
	lowest_energy = minimum(energy_in_system_new)
	last_energy = last(energy_in_system_new)
	pushfirst!(energy_in_system_new, start_energy)
	finish_time = start_datetime + Dates.Millisecond(round(last(time_s * 1000)))

	track_plot = plot(track.distance, track.altitude, title="Track",
		color=:green,
		ylabel="altitude(m)", label="track")

	speed_plot = plot(
		get_mean_data(track.distance),
		speeds,
		# seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed, proj. finish at $(finish_time)",
		ylabel="speed(kmh)",
		label="speed"
	)

	low_energy_red = fill(0., size(track.distance, 1))
	
	energy_plot = plot(
		track.distance,
		[energy_in_system_new low_energy_red],
		linewidth=[1 3],
		title="Energy with manually set speeds, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)",
		label="energy"
	)
	
	plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(750,700), legend=false)
end

function plot_resulted_plan(track, segments, minimized_speeds, energy_in_system, time_s, title_text)

    lowest_energy = minimum(energy_in_system)
	last_energy = last(energy_in_system)
    projected_finish_time = last(time_s)
    var_num = size(segments.diff_distance, 1)

    track_plot = plot(
        track.distance,
        track.altitude,
        title="Track, $(var_num) segments, $title_text",
		color=:green,
		ylabel="altitude(m)",
		xlimits=(-20, last(track.distance)+1000)
	)

	speed_plot = plot(
		get_mean_data(track.distance),
		minimized_speeds,
		# seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed, proj. finish at $(projected_finish_time)",
		ylabel="speed(kmh)",
		xlimits=(-20, last(track.distance)+1000)
	)

	low_energy_red = fill(0., size(track.distance, 1))

	energy_plot = plot(
		track.distance,
		[energy_in_system low_energy_red],
		linewidth=[1 3],
		title="Energy, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)",
		xlimits=(-20, last(track.distance)+1000)
	)
	
	return plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(700,700), legend=false)
end