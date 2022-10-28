using DataFrames

function solar_trip_calculation(input_speed::Vector{Float64}, track, 
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

function solar_trip_cost(input_speed::Vector{Float64}, track)
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
function solar_trip_chunks(speeds::Vector{<:Number}, track)
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
    show_results_wrapper(first(minimized_inputs), track);
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
