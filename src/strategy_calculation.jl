using DataFrames

include("utils.jl")

function solar_trip_calculation(input_speed, track, 
    start_energy::Float64=5100.)
    # input speed in m/s

    # calculating time needed to spend to travel across distance
    time_df = calculate_travel_time_datetime(input_speed, track)

    #### calculcations
    # mechanical calculations are now in separate file
    mechanical_power = mechanical_power_calculation(input_speed, track.slope, track.diff_distance)

    # electical losses
    electrical_power = electrical_power_calculation(track.diff_distance, input_speed)
    # converting mechanical work to elecctrical power and then power use
    # power_use = calculate_power_use(mechanical_power, electrical_power)
    power_use_accumulated_wt_h = calculate_power_use_accumulated(mechanical_power, electrical_power)

    # get solar energy income
    solar_power = solar_power_income(time_df, track, input_speed)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power)
    # TODO: night charging with additional solar panels

    # #### plotting
    # plot(track.distance, power_use, title="Power spent on toute")
    # plot(track.distance, solar_power, title="Power gained on the route")

    # plot(track.distance, power_use_accumulated_wt_h, title="Power spent on the route, accumulated")
    # plot(track.distance, solar_power_accumulated, title="Power gained on the route, accumulated")

    # plot(track.distance, solar_power_accumulated - power_use_accumulated_wt_h, title="Power balance w/o battery")
    battery_capacity = 5100 # wt, to be used later for physical constraints
    energy_in_system = start_energy .+ solar_power_accumulated .- power_use_accumulated_wt_h
    # plot(track.distance, energy_in_system, title="Power balance with battery")

    # TODO: calculate night charging - do it later since it is not critical as of right now
    # TODO: block overcharging - cost function?
    # at first do the black-box optimization, then gradient one
    # will start with Optim
    # TODO: find an optimal single speed - make a loss function and start optimization process
    time_seconds = calculate_travel_time_seconds(input_speed, track)
    # TODO: find an optimal speed vector
    return power_use_accumulated_wt_h, solar_power_accumulated, energy_in_system, time_df, time_seconds
end

function solar_trip_cost(input_speed, track)
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(input_speed, track)
    cost = last(time_s) + 10 * abs(minimum(energy_in_system)) + 100 * sum(input_speed[input_speed .< 0.0])
    return cost
end

function solar_trip_target_cost(input_speed::Vector{Float64}, target_energy::Float64, track,
    start_energy::Float64, finish_energy::Float64)
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(input_speed, track)
    cost = last(time_s) + 10 * abs(last(energy_in_system) - target) + 100 * sum(input_speed[input_speed .< 0.0])
    return cost
end

# consider using Union for types for multiple dispatch
# https://stackoverflow.com/questions/65094714/efficient-way-to-implement-multiple-dispatch-for-many-similar-functions 
# also consider rewriting core functions so they will work in .func mode

function solar_trip_wrapper(speed::Float64, track)
    # input speeds in km/h
    # as of right now one speed for all track parts
    input_speed = convert_single_kmh_speed_to_ms_vector(first(speed), length(track.distance))
    return solar_trip_cost(input_speed, track)
end

function solar_trip_wrapper(speed::Vector{Float64}, track)
    # input in km/h
    return solar_trip_cost(convert_kmh_to_ms(speed), track)
end

# function to test optimization with several big chunks to optimize
# everything under Number
function solar_trip_chunks(speeds, track)
# function solar_trip_chunks(speeds::Vector{Float64}, track)
    speed_ms = convert_kmh_to_ms(speeds)
    speed_vector = propagate_speeds(speed_ms, track)
    return solar_trip_cost(speed_vector, track)
end

function show_results_wrapper(input::Number, track::DataFrame)
    input_speed = convert_single_kmh_speed_to_ms_vector(first(input), length(track.distance))
    return show_result_graphs(input_speed, track)
end

function show_results_wrapper(inputs::Vector{Float64}, track::DataFrame)
    inputs_ms = convert_kmh_to_ms(inputs);
    speed_vector = propagate_speeds(inputs_ms, track);
    return show_result_graphs(speed_vector, track)
end

function show_result_graphs(inputs, track)
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(inputs, track);
    power_use_plot = plot(track.distance, power_use, title="Power spent on route");
    display(power_use_plot)
    power_income_plot = plot(track.distance, solar_power, title="Power gained on the route");
    display(power_income_plot)
    # plot(track.distance, power_use_accumulated_wt_h, title="Power spent on the route, accumulated")
    # plot(track.distance, solar_power_accumulated, title="Power gained on the route, accumulated")

    # plot(track.distance, solar_power - power_use, title="Power balance w/o battery")
    # battery_capacity = 5100 # wt
    # energy_in_system = battery_capacity .+ solar_power_accumulated .- power_use_accumulated_wt_h
    energy_plot = plot(track.distance, [energy_in_system zeros(size(track,1))],
    label=["Energy" "Failure threshold"], title="Energy in system", lw=3,
    xlabel="Distance (m)", ylabel="Energy (W*h)", size=(1000, 500),
    color=[:blue :red]);
    display(energy_plot)

    plot(track.distance, track.altitude, label="altitude", ylabel="altitude", title="Speed (m/s) vs distance")
    speed_distance_plot = plot!(twinx(), inputs, color=:red, ylabel="speed", label="speed (m/s)", ymirror = true, title="Speed (m/s) vs distance")
    display(speed_distance_plot)

    speed_time_plot = plot(time.utc_time, inputs, title="Speed (m/s) vs time")
    display(speed_time_plot)

    power_both_plot = plot(track.distance, [power_use solar_power energy_in_system zeros(size(track,1))],
    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph",
    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
    color=[:blue :green :cyan :red]
    # ,ylims=[-10000, 40000]
    )
    # save("energy.png", power_both_plot)
    display(power_both_plot)

    power_in_time_plot = plot(time.utc_time, [power_use solar_power energy_in_system zeros(size(track,1))],
    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph in time",
    xlabel="Time", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
    color=[:blue :green :cyan :red]
    # ,ylims=[-10000, 40000]
    )
    # save("energy.png", power_both_plot)
    display(power_in_time_plot)

    println(last(time_s));
end



function minimize_single_speed(speed::Float64, track)
    # regular optimization, Nelder-Mead, 9 seconds 
    @time result = optimize(x -> solar_trip_wrapper(x, track), [speed])
    # alert();
    minimized_inputs = Optim.minimizer(result)
    minimized_result = Optim.minimum(result)
	println("minimized speed input is $minimized_inputs kmh")
	println("minimized result time is $minimized_result seconds")
    show_results_wrapper(first(minimized_inputs), track);
	return minimized_inputs;
end

function minimize_5_chunks(track)
    # get few big chunks with the same speed, Nelder-Mead, ~210 seconds
    @time result_chunks = optimize(x -> solar_trip_chunks(x, track), [41.0, 42.0, 43.0, 44.0, 45.0])
    minimized_inputs_chunks = Optim.minimizer(result_chunks)
    minimized_result_chunks = Optim.minimum(result_chunks)
    show_results_wrapper(minimized_inputs_chunks, track);
end

function minimize(speed::Vector{Float64}, track)
    @time result = optimize(x -> solar_trip_wrapper(x, track), speed)
    minimized_inputs = Optim.minimizer(result)
    show_results_wrapper(minimized_inputs, track);
end

function split_track(track, chunks_amount, energy)
    track_len = length(track.distance)
    sector_size = track_len ÷ chunks_amount
    track_array = []
    start_energy = zeros(chunks_amount)
    finish_energy = zeros(chunks_amount)
    for i = 1:chunks_amount
        start_energy[i] = energy[(i-1)*sector_size + 1];
        # start_energy = [start_energy, energy[(i-1)*sector_size + 1]]
        # finish_energy = [finish_energy, energy[i * sector_size + 1]]
        if i==chunks_amount
            finish_energy[i] = last(energy);
            push!(track_array, track[(i-1)*sector_size + 1 : track_len, : ]);
        else
            finish_energy[i] = energy[i*sector_size];
            push!(track_array, track[(i-1)*sector_size + 1 : i * sector_size, : ]);
        end
    end
    # push!(track_array, track[track_len - last_sector_size + 1 : track_len, : ]);
    # start_energy[chunks_amount] = energy[track_len - last_sector_size + 1];
    # finish_energy[chunks_amount] = last(energy);
    return track_array, start_energy, finish_energy
end

function split_track(track, chunks_amount)
    # println("entered split_track with track array of size $(size(track, 1))")
    # show(track)
    track_len = size(track, 1)
    track_array = [];
    if track_len <= chunks_amount
        for i = 1:track_len
            push!(track_array, DataFrame(track[i,:]) )
        end
        return track_array
    end
    sector_size = track_len ÷ chunks_amount
    for i = 1:chunks_amount
        if i==chunks_amount
            push!(track_array, track[(i-1)*sector_size + 1 : track_len, : ]);
        else
            push!(track_array, track[(i-1)*sector_size + 1 : i * sector_size, : ]);
        end
    end
    return track_array;
end

function recursive_track_division(track_array, chunks_amount, total_chunks)
    # println("entering recursive_track_division with array of size $(size(track_array, 1))")
    if size(track_array, 1) == total_chunks
        return track_array;
    end
    track_parts = []
    for track in track_array
        if size(track, 1) > 1
            append!(track_parts, split_track(track, chunks_amount))
        else
            push!(track_parts, track)
        end
    end
    return recursive_track_division(track_parts, chunks_amount, total_chunks);
end

function recursive_track_optimization(track_array, chunks_amount::Int64, total_chunks::Int64,
    speed_array::Vector{Float64}, start_enegry::Float64, finish_energy::Float64)
    # println("entering recursive_track_optimization with array of size $(size(track_array, 1))")
    if size(track_array, 1) == total_chunks
        return track_array;
    end
    track_parts = []
    optimized_speeds = []
    for i in eachindex(track_array) #track in track_array
        # at first create speed inputs, for 1 track peice there is chunks_amount speeds
        speed_chunks = fill(speed_array[i], chunks_amount)
        # speed_chunks_ms = convert_kmh_to_ms(speed_chunks)
        @time result_chunks = optimize(x -> solar_trip_chunks(x, track_array[i]), speed_chunks)
        minimized_input_chunks = Optim.minimizer(result_chunks)
        minimized_result_chunks = Optim.minimum(result_chunks)
        # TODO: ensure that track splitting is identical to speed_propagation
        # OR rewrite the whole thing :)
        if size(track, 1) > 1
            append!(track_parts, split_track(track, chunks_amount))
            append!(optimized_speeds, minimized_input_chunks)
        else
            push!(track_parts, track)
            push!(optimized_speeds, minimized_input_chunks)
        end
    end
    # return additional info
    return recursive_track_division(track_parts, chunks_amount, total_chunks);
end

function minimize_hierarchical(chunks_amount::Int64, initial_speed::Number, track, 
    start_enegry::Float64, finish_energy::Float64)
    # small number of chunks - 2 to 5

    speed_chunks = fill(initial_speed, chunks_amount)
    speed_chunks_ms = convert_kmh_to_ms(speed_chunks)
    speed_vector = propagate_speeds(speed_chunks_ms, track)
    @time result_chunks = optimize(x ->solar_trip_chunks(x, track), speed_vector)
    minimized_input_chunks = Optim.minimizer(result_chunks)
    minimized_result_chunks = Optim.minimum(result_chunks)
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(minimized_input_chunks, track);

    # split track in chunks_amount
    # fix an energy amount at start and finish of chunk
    tracks, start_energies, finish_energies = split_track(track, energy_in_system, chunks_amount);
    for i in 1:chunks_amount
        minimize_hierarchical(chunks_amount, minimized_input_chunks[i], tracks[i])
        # energies?
    end
    # rewrite optimization function's cost to finish value - maybe even optimize for finish energy of 0?
    # repeat optimization process for each chunk
    # take initial speed as a speed that was used earlier

    # TODO: continue
    
end

function solar_trip_division(speed, track, divide_at)
	# input_speed = set_speeds(speed, track, divide_at)
	# return solar_trip_cost(convert_kmh_to_ms(input_speed), track)
    return solar_trip_cost(convert_kmh_to_ms(speed), track)
end

function set_speeds_293(speeds, track, divide_at)
	println("set speeds line 293")
	output_speeds = fill(last(speeds), size(track.distance, 1))
	for i=1:size(divide_at,1)-1
		if i==1
			output_speeds[1:divide_at[1]] .= speeds[1]
		else
			output_speeds[divide_at[i-1] + 1:divide_at[i]] .= speeds[i]
		end
	end
	return output_speeds
end

function recursive_optimization(speed, track, track_raw, divide_at, iteration)
    track = flatten_track(track_raw, divide_at);
	@time result = optimize(x -> solar_trip_division(x, track, divide_at), speed)

	minimized = Optim.minimizer(result)
	println(minimized)
	inputs = convert_kmh_to_ms(minimized)
	speed_vec = set_speeds(inputs, track_raw, divide_at)
	power_use_div, solar_power_div, energy_in_system_div, time_div, time_s_div = solar_trip_calculation(speed_vec, track_raw)

	println(last(time_s_div))

	zero_index_div = argmin(abs.(energy_in_system_div))

	new_divide = copy(divide_at)
	push!(new_divide, zero_index_div)
	sort!(new_divide)
	new_divide_index = findfirst(x -> x==zero_index_div, new_divide)
	println(new_divide)

	new_track = flatten_track(track_raw, new_divide)

	if iteration == run_until
		return power_use_div, solar_power_div, energy_in_system_div, time_div, time_s_div;
	end
	
	new_speed = copy(minimized)
	insert!(new_speed, new_divide_index, 40.)
	recursive_optimization(new_speed, new_track, track_raw, new_divide, iteration + 1, run_until)
	
	return power_use_div, solar_power_div, energy_in_system_div, time_div, time_s_div;
end

function flatten_track(track_raw, divide_indexes)
	if size(divide_indexes,1) == 0
		return track_raw
	end
	track = track_raw[divide_indexes,:]
	distance = copy(track.distance)
	altitude = copy(track.altitude)
	pushfirst!(distance,track_raw[1,:distance])
	pushfirst!(altitude,track_raw[1,:altitude])
	track.diff_distance = diff(distance)
	track.diff_altitude = diff(altitude)
	track.slope = atand.(track.diff_altitude ./ track.diff_distance)
	return track
end

function calculate_split_indexes(size_to_distribute, chunks_amount)
	if size_to_distribute <= chunks_amount
		return collect(1:size_to_distribute)
	end
	# TODO: split more evenly, if size%chunks < size/chunks/2 then use floor, otherwise - ceiling
	chunks_indexes = zeros(Int, chunks_amount)
	# if size_to_distribute % chunks_amount < size_to_distribute / chunks_amount / 2
		step_size = floor(Int, size_to_distribute / chunks_amount)
	# else
	# 	step_size = ceil(Int, size_to_distribute / chunks_amount)
	# end
	#div_ceiling = ceil(Int, size_to_distribute / chunks_amount)
	accumulator = 0
	for i=1:chunks_amount - 1
		accumulator += step_size
		chunks_indexes[i] = accumulator
	end
	chunks_indexes[chunks_amount] = size_to_distribute
	return chunks_indexes
end

function calculate_split_indexes(start_index, end_index, chunks_amount)
	return start_index .-1 .+ calculate_split_indexes(end_index - start_index + 1, chunks_amount)
end

function calculate_split_bounds_segments(start_point, end_point, chunks_amount)
	# участки - segments, sections 
	length_in_segments = end_point - start_point + 1;

	splitted_segments_length = calculate_split_indexes(length_in_segments, chunks_amount) .+ start_point .- 1
	segments_bounds = []
	# segments_bounds_array = zeros(size(splitted_segments_length, 1), 2)

	push!(segments_bounds, (start_point, splitted_segments_length[1]))
	for i=2:length(splitted_segments_length)
		push!(segments_bounds, (splitted_segments_length[i-1] + 1, splitted_segments_length[i]))
	end
	# return array of tuples?
	return segments_bounds
end

function calculate_boundaries(start_point, end_point, segments_amount)::Vector{Boundaries}
	# участки - segments, sections 
	length_in_segments = end_point - start_point;
	length_in_points = end_point - start_point + 1;

	# calculates split points (without start)
	split_points = calculate_split_indexes(length_in_segments, segments_amount) .+ start_point
	# so, we are adding start
	pushfirst!(split_points, start_point)

	boundaries_array::Vector{Boundaries} = []
	for i=2:length(split_points)
		segment = Boundaries(split_points[i-1], split_points[i], split_points[i] - split_points[i-1]);
		push!(boundaries_array, segment)
	end
	# return array of tuples?
	return boundaries_array
end

function calculate_split_bounds_points(start_point, end_point, chunks_amount)
	# участки - segments, sections 
	length_in_points = end_point - start_point + 1;

	splitted_segments_length = calculate_split_indexes(length_in_points, chunks_amount)
	points_bounds = []
	segments_bounds_array = zeros(size(splitted_segments_length, 1), 2)

	push!(points_bounds, (1, splitted_segments_length[1]))
	for i=2:length(splitted_segments_length)
		push!(points_bounds, (splitted_segments_length[i-1], splitted_segments_length[i]))
	end
	# return array of tuples?
	return points_bounds
end

function split_track_by_indexes(track, indexes)
	current_index = 1
	results = []
	for index in indexes
		# println("from $(current_index) to $(index)")
		push!(results, track[current_index:index, :])
		current_index = index + 1
	end
	return results
end

function solar_trip_calculation_bounds(input_speed, track, start_datetime,
    start_energy::Float64=5100.)
    # input speed in m/s
	# @debug "func solar_trip_calculation_bounds input_speed size is $(size(input_speed, 1)), track size is $(size(track.distance, 1)))"

    # calculating time needed to spend to travel across distance
    time_df = calculate_travel_time_datetime(input_speed, track, start_datetime)

    #### calculcations
    # mechanical calculations are now in separate file
    mechanical_power = mechanical_power_calculation(input_speed, track.slope, track.diff_distance)

    # electical losses
    electrical_power = electrical_power_calculation(track.diff_distance, input_speed)
    # converting mechanical work to elecctrical power and then power use
    # power_use = calculate_power_use(mechanical_power, electrical_power)
    power_use_accumulated_wt_h = calculate_power_use_accumulated(mechanical_power, electrical_power)

    # get solar energy income
	# @debug "track size is $(size(track.latitude, 1))"
    solar_power = solar_power_income(time_df, track, input_speed)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power)
    # TODO: night charging with additional solar panels

    # #### plotting
    # plot(track.distance, power_use, title="Power spent on toute")
    # plot(track.distance, solar_power, title="Power gained on the route")

    # plot(track.distance, power_use_accumulated_wt_h, title="Power spent on the route, accumulated")
    # plot(track.distance, solar_power_accumulated, title="Power gained on the route, accumulated")

    # plot(track.distance, solar_power_accumulated - power_use_accumulated_wt_h, title="Power balance w/o battery")
    # battery_capacity = 5100 # wt, to be used later for physical constraints
    energy_in_system = start_energy .+ solar_power_accumulated .- power_use_accumulated_wt_h
    # plot(track.distance, energy_in_system, title="Power balance with battery")

    # TODO: calculate night charging - do it later since it is not critical as of right now
    # TODO: block overcharging - cost function?
    # at first do the black-box optimization, then gradient one
    # will start with Optim
    # TODO: find an optimal single speed - make a loss function and start optimization process
    time_seconds = calculate_travel_time_seconds(input_speed, track)
    # TODO: find an optimal speed vector
    return power_use_accumulated_wt_h, solar_power_accumulated, energy_in_system, time_df, time_seconds
end

function solar_trip_calculation_bounds_alloc(input_speed, track, start_datetime,
    start_energy::Float64=5100.)
    # input speed in m/s
	# @debug "func solar_trip_calculation_bounds input_speed size is $(size(input_speed, 1)), track size is $(size(track.distance, 1)))"

    # calculating time needed to spend to travel across distance
    time_df = calculate_travel_time_datetime(input_speed, track, start_datetime)

    #### calculcations
    # mechanical calculations are now in separate file
    mechanical_power = mechanical_power_calculation_alloc.(input_speed, track.slope, track.diff_distance)

    # electical losses
    electrical_power = electrical_power_calculation(track.diff_distance, input_speed)
    # converting mechanical work to elecctrical power and then power use
    # power_use = calculate_power_use(mechanical_power, electrical_power)
    power_use_accumulated_wt_h = mechanical_power + electrical_power
	cumsum!(power_use_accumulated_wt_h, power_use_accumulated_wt_h)
	power_use_accumulated_wt_h = power_use_accumulated_wt_h / 3600.

    # get solar energy income
	# @debug "track size is $(size(track.latitude, 1))"
    solar_power = solar_power_income_alloc.(
		track.latitude,
		track.longitude, 
		track.altitude, 
		time_df.utc_time,
		track.diff_distance,
		input_speed
	)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power)
    energy_in_system = start_energy .+ solar_power_accumulated .- power_use_accumulated_wt_h

    # TODO: calculate night charging - do it later since it is not critical as of right now
    time_seconds = calculate_travel_time_seconds(input_speed, track)
    return power_use_accumulated_wt_h, solar_power_accumulated, energy_in_system, time_df, time_seconds
end

function set_speeds(speeds, track, divide_at)
	println("set speed line 507")
	output_speeds = fill(last(speeds), size(track.distance, 1))
	for i=1:size(divide_at,1)-1
		if i==1
			output_speeds[1:divide_at[1]] .= speeds[1]
		else
			output_speeds[divide_at[i-1]:divide_at[i]] .= speeds[i]
		end
	end
	return output_speeds
end

function set_speeds_segments(speeds, track_len, segments)
	output_speeds = fill(last(speeds), track_len)
	for i=1:size(segments, 1)
		output_speeds[segments[i][1]:segments[i][2]] .= speeds[i]
	end
	return output_speeds
end

function set_speeds_segments(speeds, segments)
	track_len = last(segments)[2] - first(segments)[1] + 1
	first_index = first(segments)[1]
	output_speeds = fill(last(speeds), track_len)
	for i=1:size(segments, 1)
		output_speeds[segments[i][1] - first_index + 1 : segments[i][2] - first_index + 1] .= speeds[i]
	end
	return output_speeds
end

# differentiable with ReverseDiff.jl
function set_speeds_grad(speeds, track, divide_at)
	# output_speeds = fill(last(speeds), size(track.distance, 1))
	output_speeds = fill(speeds[1], divide_at[1])
	for i=2:size(divide_at,1)
		output_speeds = vcat(output_speeds, fill(speeds[i], divide_at[i] - divide_at[i-1]))
	end
	return output_speeds
end

function solar_partial_trip_wrapper(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds(speeds_ms, track, indexes)
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) + 10 * abs(finish_energy - last(energy_in_system))
    return cost
end

function solar_partial_trip_wrapper_alloc(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds(speeds_ms, track, indexes)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds_alloc(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) + 10 * abs(finish_energy - last(energy_in_system))
	return cost
	# return solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
end

function solar_partial_trip_wrapper_iter(speeds, track, segments, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds_segments(speeds_ms, size(track.distance,1), segments)
	power_use, solar_power, energy_in_system, time, time_s = 
		solar_trip_calculation_bounds_alloc(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) +  (finish_energy - last(energy_in_system))^2
	return cost
	# return solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
end

function hierarchical_optimization(
    speed, track, chunks_amount, start_energy, finish_energy, start_datetime,
    iteration, end_index
    )
	# 0. if track is non-divisible on chunks_amount, then return (array of speeds)
	# 1. split the whole track in chunks (chunks division and speed propagation with same logic - divide at the same idexes)
	# 2. optimize it on chunks (initial speed = speed, use it for all chunks)
	# 3. save chunks_amount input speeds
	# 4. simulate it again to get energy levels at start and finish of each chunk
	# 5. go through resulting speeds and track chunks to optimize them (entering recursion)

	# 0 - exit condition, stub for now
	# if iteration == 5
	# 	return speed
	# end

	# @debug "func hierarchical_optimization speed is $(speed), track_size is $(size(track.distance, 1))"
	
	# track is non-divisible, if its size is <= 1, return speed
	if size(track.distance, 1) == 1
		return speed
	end

	# 1 - splitting the track
	# determine split indexes
	track_size = size(track.distance, 1)
    start_index = end_index - track_size + 1
	split_indexes = calculate_split_indexes(track_size, chunks_amount)
	# @debug "split indexes are $(split_indexes), chunks are $(chunks_amount)"
	# actually split the track
	tracks = split_track_by_indexes(track, split_indexes)
	# for the case when there are less indexes than chunks
	chunks_amount = size(split_indexes,1)

	# 2 - set up optimization itself
	function f(speed)
		return solar_partial_trip_wrapper(
            abs.(speed), track, split_indexes, start_energy, finish_energy,
            start_datetime
            )
	end
	td = TwiceDifferentiable(f, fill(speed, chunks_amount); autodiff = :forward)
	lower_bound = fill(0.0, chunks_amount)
	upper_bound = fill(100.0, chunks_amount)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
	line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, chunks_amount),
	    #Newton(; linesearch = line_search),
	result = optimize(
        td, tdc, fill(speed, chunks_amount)
        .+ (rand(chunks_amount) .* 0.5)
        ,
		IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10
	    )
	)

	# 3 - save optimized speeds
	minimized_speeds = abs.(Optim.minimizer(result))
	
	# 4 - sumulate again to obtain energies and times around split indexes
	minimized_speeds_ms = convert_kmh_to_ms(minimized_speeds)
	minimized_speed_vector = set_speeds(minimized_speeds_ms, track, split_indexes)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds(
        minimized_speed_vector, track, start_datetime, start_energy)
	# println("iteration $(iteration), speed is $(speed) planned finish energy is $(finish_energy)")
	# println("time is $(last(time_s)), cost is $(f(minimized_speeds))")
	# println("minimized speeds are: $(minimized_speeds)")
	# println("simulated finish energy is $(last(energy_in_system))")
	# # println("calculated cost is $( last(time_s) + 100 * abs(last(energy_in_system) - finish_energy) + 100 * sum(abs.(energy_in_system[energy_in_system .< 0.0])) + 100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .< 0.0])) + 100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .> 100.0])) )")
	# println("finish energy difference penalty is: $(10 * abs(last(energy_in_system) - finish_energy))")
	# # println("energy less than 0. penalty is: $(100 * sum(abs.(energy_in_system[energy_in_system .< 0.0])))")
	# # println("speed less than 0. penalty is: $(100 * sum(abs.(speed_vector[speed_vector .< 0.0])))")
	# # println("speed more than 100. penalty is: $(100 * sum(abs.(speed_vector[speed_vector .> 100.0 / 3.6])))")
	split_energies = energy_in_system[split_indexes]
	pushfirst!(split_energies, start_energy)
	split_times = time[split_indexes, :utc_time]
	pushfirst!(split_times, start_datetime)
	# println("split energies are $(split_energies)")
    println("")
	
	# 5 - go though each track piece and enter function again
	# hierarchical_optimization(minimized_speeds[1], tracks[1], chunks_amount, start_energy, split_energies[1], start_datetime, iteration + 1)
	# @debug "split_energies size is $(size(split_energies, 1)), chunks_amount is $(chunks_amount)"
	result_speeds = []
	for i=1:chunks_amount
		result_speeds_chunk = hierarchical_optimization(
            minimized_speeds[i], tracks[i], chunks_amount, split_energies[i],
            split_energies[i+1], split_times[i], iteration + 1,
            start_index + split_indexes[i] - 1
            )
		append!(result_speeds, result_speeds_chunk)
	end

	return result_speeds
end

function hierarchical_optimization_alloc(speed, track, chunks_amount, start_energy, finish_energy, start_datetime, iteration, end_index)
	# 0. if track is non-divisible on chunks_amount, then return (array of speeds)
	# 1. split the whole track in chunks (chunks division and speed propagation with same logic - divide at the same idexes)
	# 2. optimize it on chunks (initial speed = speed, use it for all chunks)
	# 3. save chunks_amount input speeds
	# 4. simulate it again to get energy levels at start and finish of each chunk
	# 5. go through resulting speeds and track chunks to optimize them (entering recursion)

	# 0 - exit condition, stub for now
	# if iteration == 5
	# 	return speed
	# end

	# @debug "func hierarchical_optimization speed is $(speed), track_size is $(size(track.distance, 1))"
	
	# track is non-divisible, if its size is <= 1, return speed
	if size(track.distance, 1) == 1
		return speed
	end

	# 1 - splitting the track
	# determine split indexes
	track_size = size(track.distance, 1)
	start_index = end_index - track_size + 1
	split_indexes = calculate_split_indexes(track_size, chunks_amount)
	# @debug "split indexes are $(split_indexes), chunks are $(chunks_amount)"
	# actually split the track
	tracks = split_track_by_indexes(track, split_indexes)
	# for the case when there are less indexes than chunks
	chunks_amount = size(split_indexes,1)

	# 2 - set up optimization itself
	function f(speed)
		return solar_partial_trip_wrapper_alloc(speed, track, split_indexes, start_energy, finish_energy, start_datetime)
	end
	td = TwiceDifferentiable(f, fill(speed, chunks_amount); autodiff = :forward)
	lower_bound = fill(0.0, chunks_amount)
	upper_bound = fill(100.0, chunks_amount)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
	# line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, chunks_amount),
	    #Newton(; linesearch = line_search),
	result = optimize(td, tdc, fill(speed, chunks_amount) 
	.+ (rand(chunks_amount) .* 0.5)
		,
		IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10
	    )
	)

	# 3 - save optimized speeds
	minimized_speeds = abs.(Optim.minimizer(result))
	
	# 4 - sumulate again to obtain energies and times around split indexes
	minimized_speeds_ms = convert_kmh_to_ms(minimized_speeds)
	minimized_speed_vector = set_speeds(minimized_speeds_ms, track, split_indexes)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds_alloc(minimized_speed_vector, track, start_datetime, start_energy)
	# println("iteration $(iteration), speed is $(speed) planned finish energy is $(finish_energy)")
	# println("track from $(start_index) to $(end_index)")
	# println("start datetime $(start_datetime)")
	# println("stare energy $(start_energy)")
	# println("split indexes are $(split_indexes)")
	# # println("distances are $(track[split_indexes, :])")
	# # println("minimized speeds are: $(minimized_speeds)")
	# println("simulated finish energy is $(last(energy_in_system))")
	# # println("finish energy difference penalty is: $(100 * abs(last(energy_in_system) - finish_energy))")
	split_energies = energy_in_system[split_indexes]
	pushfirst!(split_energies, start_energy)
	split_times = time[split_indexes, :utc_time]
	pushfirst!(split_times, start_datetime)
	# println("split energies are $(split_energies)")
	# println("")
	
	# 5 - go though each track piece and enter function again
	# @debug "split_energies size is $(size(split_energies, 1)), chunks_amount is $(chunks_amount)"
	result_speeds = []
	for i=1:chunks_amount
		result_speeds_chunk = hierarchical_optimization_alloc(minimized_speeds[i], tracks[i], chunks_amount, split_energies[i], split_energies[i+1], split_times[i], iteration + 1 , start_index + split_indexes[i] - 1)
		append!(result_speeds, result_speeds_chunk)
	end

	return result_speeds
end

function hierarchical_optimization_alloc!(speed_by_iter, speed, track, chunks_amount, start_energy, finish_energy, start_datetime, iteration, end_index)
	# 0. if track is non-divisible on chunks_amount, then return (array of speeds)
	# 1. split the whole track in chunks (chunks division and speed propagation with same logic - divide at the same idexes)
	# 2. optimize it on chunks (initial speed = speed, use it for all chunks)
	# 3. save chunks_amount input speeds
	# 4. simulate it again to get energy levels at start and finish of each chunk
	# 5. go through resulting speeds and track chunks to optimize them (entering recursion)

	# 0 - exit condition, stub for now
	# if iteration == 5
	# 	return speed
	# end

	# @debug "func hierarchical_optimization speed is $(speed), track_size is $(size(track.distance, 1))"
	
	# track is non-divisible, if its size is <= 1, return speed
	if size(track.distance, 1) == 1
		return speed
	end

	# 1 - splitting the track
	# determine split indexes
	track_size = size(track.distance, 1)
	start_index = end_index - track_size + 1
	split_indexes = calculate_split_indexes(track_size, chunks_amount)
	# @debug "split indexes are $(split_indexes), chunks are $(chunks_amount)"
	# actually split the track
	tracks = split_track_by_indexes(track, split_indexes)
	# for the case when there are less indexes than chunks
	chunks_amount = size(split_indexes,1)

	# 2 - set up optimization itself
	function f(speed)
		return solar_partial_trip_wrapper_alloc(speed, track, split_indexes, start_energy, finish_energy, start_datetime)
	end
	td = TwiceDifferentiable(f, fill(speed, chunks_amount); autodiff = :forward)
	lower_bound = fill(0.0, chunks_amount)
	upper_bound = fill(100.0, chunks_amount)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
	# line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, chunks_amount),
	    #Newton(; linesearch = line_search),
	result = optimize(td, tdc, fill(speed, chunks_amount) 
	.+ (rand(chunks_amount) .* 0.5)
		,
		IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10
	    )
	)

	# 3 - save optimized speeds
	minimized_speeds = abs.(Optim.minimizer(result))
	
	# 4 - sumulate again to obtain energies and times around split indexes
	minimized_speeds_ms = convert_kmh_to_ms(minimized_speeds)
	minimized_speed_vector = set_speeds(minimized_speeds_ms, track, split_indexes)
    if size(speed_by_iter, 1) < iteration
        push!(speed_by_iter, zeros( size(speed_by_iter[1], 1) ))
    end
    speed_iter_n = speed_by_iter[iteration]
    speed_iter_n[start_index:end_index] = copy(minimized_speed_vector)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds_alloc(minimized_speed_vector, track, start_datetime, start_energy)
	# println("iteration $(iteration), speed is $(speed) planned finish energy is $(finish_energy)")
	# println("track from $(start_index) to $(end_index)")
	# println("start datetime $(start_datetime)")
	# println("stare energy $(start_energy)")
	# println("split indexes are $(split_indexes)")
	# # println("distances are $(track[split_indexes, :])")
	# # println("minimized speeds are: $(minimized_speeds)")
	# println("simulated finish energy is $(last(energy_in_system))")
	# # println("finish energy difference penalty is: $(100 * abs(last(energy_in_system) - finish_energy))")
	split_energies = energy_in_system[split_indexes]
	pushfirst!(split_energies, start_energy)
	split_times = time[split_indexes, :utc_time]
	pushfirst!(split_times, start_datetime)
	# println("split energies are $(split_energies)")
	# println("")
	
	# 5 - go though each track piece and enter function again
	# @debug "split_energies size is $(size(split_energies, 1)), chunks_amount is $(chunks_amount)"
	result_speeds = []
	for i=1:chunks_amount
		result_speeds_chunk = hierarchical_optimization_alloc!(speed_by_iter,
        minimized_speeds[i], tracks[i], chunks_amount,
        split_energies[i], split_energies[i+1], split_times[i],
        iteration + 1 , start_index + split_indexes[i] - 1)
		append!(result_speeds, result_speeds_chunk)
	end

	return result_speeds
end


mutable struct Boundaries
	from::Integer
	to::Integer
	size::Integer
	# points::DataFrame # maybe better another sub-struct or arrays?
	# segments::DataFrame
	# there is no sense in placing points(track) and segments here, 
	# since they will need another indexing
	# it is better to just use indexes

	# overloaded constructor will only work without @proto macro
	Boundaries(from, to) = new(from, to, to - from)
	Boundaries(from, to, size) = new(from, to, size)
end


@proto mutable struct SubtaskProblem
	start_energy::AbstractFloat
	finish_energy::AbstractFloat
	initial_speed::AbstractFloat
	start_datetime::DateTime
end

mutable struct SubtaskVariable
	boundaries::Boundaries
	speed::AbstractFloat
end

@proto struct SubtaskSolution
	speeds::Vector{AbstractFloat} # n-1 speeds (for segments)
	energies::Vector{AbstractFloat} # n energies (for points)
	times::Vector{AbstractFloat} # n times (for points)
end

mutable struct Subtask
	subtask_boundaries::Boundaries
	variables_boundaries::Vector{Boundaries}
	problem::SubtaskProblem
	solution::SubtaskSolution
end

mutable struct Iteration
	subtasks::Vector{Subtask}
	number::Integer
end

function iterative_optimization_new(track, segments, scaling_coef, start_energy)
	# general algorithm:
	# 0. data setup
	# 1. exit loop check
	# 2. 	split the track into subtasks
	# 3. 		each subtask comes with its own chunks (variables)
	# 4. 		solve optimization problem for every subtask with its chunks (variables)
	# 5. 	tie everything together (collect speeds to one array)
	# 6. 	make full simulation with said speeds
	# 7. 	prepare data for next iteration. should be subtask's chunks as subtasks
	# 8. 	go to 2:
	# 9. final calculations?

	# 0. data setup
	track_size = size(track,1)
	# initial splitting not needed
	# array of arrays. why? 1st level for iterations, 2nd level for substasks on interation
	subtasks_splits_general = []
	push!(subtasks_splits_general, calculate_split_bounds_segments(1,track_size, 1) )
	# variables_splits_iteration = []
	# will be initialized in loop
	# subtasks_splits_segments =  [ calculate_split_bounds_segments(track_size, scaling_coef) ]

	# there should be as many speeds, as variables on last iteration
	speeds_general = []
	push!(speeds_general, [ [43.97116943001747] ] )
	# push!(speeds_general, minimize_single_speed(50., track) )


	start_energies_general = []
	push!(start_energies_general, [ [start_energy] ])

	zero_subtask = Subtask(
		Boundaries(1, track_size, track_size - 1),
		# calculate_boundaries(1, track_size, scaling_coef),
		[],
		SubtaskProblem(
			start_energy,
			0.,
			43.97116943001747,
			DateTime(2022,7,1,0,0,0)
		),
		SubtaskSolution(
			[],
			[],
			[]
		)
	)

	iteration_1 = Iteration(
		[ zero_subtask ],
		1
	);

	iterations::Vector{Iteration} = []
	push!(iterations, iteration_1)

	iteration_num = 1;
	# 1. exit loop check
	is_track_divisible_further = true
	# while iteration <= 2
	while is_track_divisible_further
		# amount_of_subtasks = scaling_coef^(iteration - 1) # also amount of variables from prev. iteration
		# amount_of_variables = scaling_coef^iteration

		iteration = iterations[iteration_num]
		println("Iteration $(iteration.number)")
		
		# 2. split the track into subtasks
		# or grab the result of previous division for subtasks
		# подумать где у меня точки, где участки, а где трассы по результатам деления
		# и ещё подумать как получаются подзадачи и их переменные на разных итерациях
		# № итерации, подзадачи, количество переменных
		# 1			1			scale
		# 2			scale		scale^2
		# ...
		# n			scale^(n-1)	scale^n - но не совсем, где как места хватит

		# по идее достаточно делить трассу, когда определяемся с участками для подзадач
		# потом в конце итерации после всех оптимизаций склеивать единый массив из точек разделения
		# и на следующей итерации использовать уже его

		# # generates separate track dataframes
		# tracks_iteration_chunk = split_track_by_indexes(
		# 		track_iteration_chunk,
		# 		split_indexes_chunk
		# 		)

		is_track_divisible_further = false
		# is_there_single_subtask_where_track_is_divisible = false
		for subtask_index in eachindex(iteration.subtasks)
		# for subtask_index in eachindex(subtasks_splits_iteration)
			# println("Subtask $subtask_index")
			subtask = iteration.subtasks[subtask_index]
			# println("Analyzing subtask from $(subtask.subtask_boundaries.from) to $(subtask.subtask_boundaries.to)")

			# split each task on parts
			subtask.variables_boundaries = calculate_boundaries(
				subtask.subtask_boundaries.from,
				subtask.subtask_boundaries.to,
				scaling_coef
			)

			if (subtask.subtask_boundaries.size) >= scaling_coef
				# at least one subtask can be divided, so, there should be next iteration
				is_track_divisible_further = true
			end

			# # check if track is divisible further
			# if (subtask.subtask_boundaries.size) >= scaling_coef
			# 	is_track_divisible_further = true

			# 	# split each task on parts
			# 	subtask.variables_boundaries = calculate_boundaries(
			# 		subtask.subtask_boundaries.from,
			# 		subtask.subtask_boundaries.to,
			# 		scaling_coef
			# 	)

			# 	# append!(subtask.variables_segments, subtask_variables_segments)
			# 	# fill it later straight with results

			# 	# TODO : check here if could be divided further
			# 	for variable_index in eachindex(subtask_variables_split)
			# 		if (subtask_variables_split[variable_index][2] - 
			# 			subtask_variables_split[variable_index][1] + 1) >= scaling_coef
			# 			# is_there_single_subtask_where_track_is_divisible = true
			# 			is_track_divisible_further = true
			# 		end
			# 	end

			# 	for variable_segment in subtask_variables_segments
			# 		if (variable_segment.to - variable_segment.from + 1) >= scaling_coef
			# 			# is_there_single_subtask_where_track_is_divisible = true
			# 			is_track_divisible_further = true
			# 		end
			# 	end


			# 	# collecting the variables split for the next iteration
			# 	append!(variables_split_iteration, subtask_variables_split)
			# 	# append!(variables_segments_iteration, subtask_variables_segments)
			# 	# exit codition here when we can no longer split?
			# 	# new splits only here
			# else
			# 	push!(variables_split_iteration, subtask_splits)
			# 	# push!(variables_segments_iteration, subtask.subtask_segment)
			# 	subtask_variables_split = subtask_splits

			# 	# use the same subtask on next iter?
			# end

			# perform optimization

			# optimize
			# here or only if segments were changed?
			# it is better to optimize everytime,
			# since we are improving accuracy every iteration
			
			# TODO: get needed parametes for optimization

			#######################

			
			# # TODO: how to calculate amount of speeds?
			# input_speeds = fill(subtask.problem.initial_speed, length(subtask.variables_boundaries))

			# track_subtask = [subtask.subtask_boundaries.from:subtask.subtask_boundaries.to,:]
			# # finish_energy_subtask = 
			# # start_datetime_subtask = 
		
			# # TODO: stopped here, check how variables_boundaries used inside the function
			# function f_iter(input_speeds)
			# 	return solar_partial_trip_wrapper_iter(
			# 		input_speeds, track_subtask, subtask.variables_boundaries,
			# 		subtask.problem.start_energy, subtask.problem.finish_energy,
			# 		subtask.problem.start_datetime
			# 	)
			# end

			# td = TwiceDifferentiable(f_iter, fill(speed, chunks_amount); autodiff = :forward)
			# lower_bound = fill(0.0, chunks_amount)
			# upper_bound = fill(100.0, chunks_amount)
			# tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
			# # line_search = LineSearches.BackTracking();
			# # result = optimize(td, fill(speed, chunks_amount),
			# 	#Newton(; linesearch = line_search),
			# result = optimize(td, tdc, fill(speed, chunks_amount) 
			# .+ (rand(chunks_amount) .* 0.5)
			# 	,
			# 	IPNewton(),
			# 	Optim.Options(
			# 		x_tol = 1e-10,
			# 		f_tol = 1e-10,
			# 		g_tol = 1e-10
			# 	)
			# )
			# minimized_speeds = Optim.minimizer(result)

			# TODO: check optimization procedure in compliance with article
			# TODO: save result somewhere - in subtask struct
			# OR, in subtaskResult struct

			# 3. each subtask comes with its own chunks (variables)
			# 4. solve optimization problem for every subtask with its chunks (variables)
		end
		# is_track_divisible_further = !is_there_single_subtask_where_track_is_divisible
		# push!(subtasks_splits_general, variables_split_iteration)

		println()
		# 5. tie everything together (collect speeds to one array)
		# 6. make full simulation with said speeds
		# 7. 	prepare data for next iteration. should be subtask's chunks as subtasks

		if is_track_divisible_further
			iteration_num += 1;
			next_iteration = Iteration(
				[],
				iteration_num
			)
			for subtask in iteration.subtasks
				for variable_boundaries in subtask.variables_boundaries
					new_subtask = Subtask(
						variable_boundaries,
						[],
						SubtaskProblem(
							0.,
							0.,
							0.,
							DateTime(2022,1,1,0,0,0)
						),
						SubtaskSolution(
							[],
							[],
							[]
						)
					)
					push!(next_iteration.subtasks, new_subtask)
				end
			end

			push!(iterations, next_iteration)
		end
	end # 8. 	go to 2:
	# 9. final calculations?
	println("Calc done")
	return iterations

end

function iterative_optimization(
	speed_by_iter, speed, track, chunks_amount, start_energy, finish_energy,
	start_datetime, end_index
	)
	# same approach as optimization_hierarchical, but we are going iter-by-iter

	# and also at the end of every iteration (when all the subroutines end),
	# perform a run with current inputs (result from current iteration depth),
	# so that we obtain more accurate time and energy

	# 0. if track is non-divisible on chunks_amount, then return (array of speeds)
	# 1. split the whole track in chunks (chunks division and speed propagation with same logic - divide at the same idexes)
	# 2. optimize it on chunks (initial speed = speed, use it for all chunks)
	# 3. save chunks_amount input speeds
	# 4. simulate it again to get energy levels at start and finish of each chunk
	# 5. go through resulting speeds and track chunks to optimize them (entering recursion)


	# data set-up for iterative variant
	split_indexes_array = []
	start_indexes = []
	end_indexes = []
	start_energies = []
	push!(start_energies, [start_energy])
	finish_energies = []
	push!(finish_energies, [finish_energy])
	start_datetimes = []
	push!(start_datetimes, [start_datetime])
	iteration = 1
	track_sizes = []
	tracks = []
	push!(tracks, [track])
	speeds_by_iter = []
	push!(speeds_by_iter, [speed])
	# should we store reduced data (only for variables/chunks), or full propagated array
	# say, we have total track length of 5, and on 1st iteration should we only store "speed",
	# or "speed-speed-speed-speed-speed"
	# i thnink it's better to store the latter
	# or maybe both?

	# general algorithm:
	# 1. exit loop check
	# 2. 	split the track into subtasks
	# 3. 		each subtask comes with its own chunks (variables)
	# 4. 		solve optimization problem for every subtask with its chunks (variables)
	# 5. 	tie everything together (collect speeds to one array)
	# 6. 	make full simulation with said speeds
	# 7. 	prepare data for next iteration. should be subtask's chunks as subtasks
	# 8. 	go to 2:
	# 9. final calculations?

	# on every iteration we will perform chunks_amount^(iteration-1) optimizations
	# TODO: formulate an exit condition
	while iteration <= 2
		println("Iteration $iteration")
		
		split_indexes = []
		iter_speeds = []
		speeds_vectors = []
		tracks_on_iter = []
		start_energies_iter = []
		finish_energies_iter = []
		chunks_to_divide = []
		start_indexes_iter = []
		end_indexes_iter = []
		# now on iteration X we have to perform optimization 
		# for every chunk on this iteration
		# TODO: it will fail on later stages. or not?
		for subtask in 1:chunks_amount^(iteration - 1)
			# optimization happens here

			println("Subtask $subtask")
			# chunk_number -> subtask
			# TODO: change variable names
			track_iteration_chunk = tracks[iteration][subtask]
			start_energy_chunk = start_energies[iteration][subtask]
			finish_energy_chunk = finish_energies[iteration][subtask]
			start_datetime_chunk = start_datetimes[iteration][subtask]
			speed_chunk = speeds_by_iter[iteration][subtask]
			track_size = size(track_iteration_chunk, 1)

			
			split_indexes_chunk = calculate_split_indexes(track_size, chunks_amount)
			push!(split_indexes, split_indexes_chunk)

			# why do i need this?
			# tracks on this exact iteration and chunk, will be used for future iterations
			tracks_iteration_chunk = split_track_by_indexes(
				track_iteration_chunk,
				split_indexes_chunk
				)
			# TODO: check if we need to wrap in array
			push!(tracks_on_iter, tracks_iteration_chunk)
			# in case if track size is less than chunks_amount 
			chunks_amount_current = size(split_indexes_chunk, 1)
			push!(chunks_to_divide, chunks_amount_current)

			function f_iter(input_speed)
				return solar_partial_trip_wrapper_alloc(
					input_speed, track_iteration_chunk, split_indexes_chunk,
					start_energy_chunk, finish_energy_chunk, start_datetime_chunk)
			end

			# td = TwiceDifferentiable(f_iter, fill(speed_chunk, chunks_amount_current); autodiff = :forward)
			# lower_bound = fill(0.0, chunks_amount_current)
			# upper_bound = fill(100.0, chunks_amount_current)
			# tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
			# # line_search = LineSearches.BackTracking();
			# # result = optimize(td, fill(speed, chunks_amount),
			# 	#Newton(; linesearch = line_search),
			# result_chunk = optimize(td, tdc, fill(speed_chunk, chunks_amount_current) 
			# .+ (rand(chunks_amount_current) .* 0.5)
			# 	,
			# 	IPNewton(),
			# 	Optim.Options(
			# 		x_tol = 1e-10,
			# 		f_tol = 1e-10,
			# 		g_tol = 1e-10
			# 	)
			# )

			# minimized_speeds_chunk = Optim.minimizer(result_chunk)

			# minimized_speeds_ms = convert_kmh_to_ms(minimized_speeds_chunk)
			# minimized_speed_vector = set_speeds(minimized_speeds_ms,
			# 	track_iteration_chunk, split_indexes_chunk
			# 	)
			# push!(speeds_vectors, minimized_speed_vector)
			# append!(iter_speeds, minimized_speed_vector)
			# now we have speeds for certain part of the track
			# gotta save 'em to later calculate whole trip for iteration
			# write resulting speeds
			# power_use, solar_power, energy_in_system, time,
			# time_s = solar_trip_calculation_bounds_alloc(minimized_speed_vector,
			# 	track_iteration_chunk, start_datetime_chunk, start_energy_chunk
			# 	)
			# save 

			# # no, we shouldn't really be doing that
			# split_energies = energy_in_system[split_indexes_chunk]
			# pushfirst!(split_energies, start_energy_chunk)
			# split_times = time[split_indexes_chunk, :utc_time]
			# pushfirst!(split_times, start_datetime_chunk)

			# somewhere here goes pushing tracks_iteration_chunk for future iterations
			

		end

		println()

		# calculating the run after the iteration
		power_use, solar_power, energy_in_system, time,
		time_s = solar_trip_calculation_bounds_alloc(iter_speeds,
			track, start_datetime, start_energy
			)

		# now we have to prepare data for next iteration
		# smaller tracks with different energies and start times
		# for each subtrack: start_e, finish_e, start_datetime
		for i in eachindex(split_indexes)
			# TODO: rework this code piece since we have 1 sim at the end of iteration
			# so, we will only need to add to start energy once?
			indexes_ith_chunk = split_indexes[i]
			# think of start and end index in order to split track

			# split indexes are local, not absolute
			# so an addition of start_index is needed


			split_energies = energy_in_system[split_indexes[i]]
			pushfirst!(split_energies, start_energy[iteration][i])
			split_times = time[split_indexes[i], :utc_time]
			pushfirst!(split_times, start_datetimes[iteration][i])
			# for each chunk we have chunks_amount tracks
			for j in chunks_to_divide

			end
		end

		# going the easy way
		# as simple as possible
		# loop through every big chunk we have 
		# we can loop through split_indexes since it is array of arrays
		for chunk_number in 1:chunks_amount^(iteration - 1)

			# and extract energy and time for every small chunk
			# and store it for next iteration
		end

		push!(tracks, tracks_on_iter)
		# code below just as an example
		# for chunk_number in 1:chunks_amount^(iteration - 1)
		# 	split_energies = energy_in_system[split_indexes_chunk]
		# 	pushfirst!(split_energies, start_energy_chunk)
		# 	split_times = time[split_indexes_chunk, :utc_time]
		# 	pushfirst!(split_times, start_datetime_chunk)
		# end

		push!(split_indexes_array, split_indexes)
		
		# TODO: adjust

		iteration += 1
	end

end