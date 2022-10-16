module SolarStrategy

using DataFrames
# also consider using JuliaDB and Query
using CSV
using Plots # default
using PlotlyBase
using TimeZones
using Dates
using Optim
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
track_short = first(track, 5);

# TODO: track preprocessing


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
function solar_trip_test(speeds::Vector{<:Number}, track)
# function solar_trip_test(speeds::Vector{Float64}, track)
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
    energy_plot = plot(track.distance, energy_in_system, title="Power balance with battery");
    display(energy_plot)

    plot(track.distance, track.altitude, label="altitude", ylabel="altitude", title="Speed (m/s) vs distance")
    speed_distance_plot = plot!(twinx(), inputs, color=:red, ylabel="speed", label="speed (m/s)", ymirror = true, title="Speed (m/s) vs distance")
    display(speed_distance_plot)

    speed_time_plot = plot(time.utc_time, inputs, title="Speed (m/s) vs time")
    display(speed_time_plot)
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
    @time result_chunks = optimize(x -> solar_trip_test(x, track), [41.0, 42.0, 43.0, 44.0, 45.0])
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
    show(track)
    track_len = size(track, 1)
    sector_size = track_len ÷ chunks_amount
    # TODO: check how it works if size of track is smaller than chunks_amount
    track_array = [];
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
    # if size(track, 1) <= chunks_amount
    #     return [track]
    # end
    # track_parts = split_track(track, chunks_amount)
    track_parts = []
    for track in track_array
        if size(track, 1) > 1
            push!(track_parts, split_track(track, chunks_amount))
        else
            push!(track_parts, track)
        end
    end
    # return recursive_track_division(track_parts, chunks_amount)
    if size(track_parts, 1) == total_chunks
        return track_parts;
    else
        return recursive_track_division(track_parts, chunks_amount, total_chunks);
    end
end



function minimize_hierarchical(chunks_amount::Int64, initial_speed::Number, track, 
    start_enegry::Float64, finish_energy::Float64)
    # small number of chunks - 2 to 5

    speed_chunks = fill(initial_speed, chunks_amount)
    speed_chunks_ms = convert_kmh_to_ms(speed_chunks)
    speed_vector = propagate_speeds(speed_chunks_ms, track)
    @time result_chunks = optimize(x ->solar_trip_test(x, track), speed_vector)
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

# # few big chunks, LBFGS, 645 seconds, negative speeds, consider revising or constraints, 3886 seconds
# @time result_chunks_lbfgs = optimize(x -> solar_trip_test(x, track), [41.0, 42.0, 43.0, 44.0, 45.0], LBFGS())
# minimized_inputs_chunks_lbfgs = Optim.minimizer(result_chunks_lbfgs)
# minimized_result_chunks_lbfgs = Optim.minimum(result_chunks_lbfgs)
# lower = [0.0]
# upper = [100.]
# initial_x = [43.0]
# @time result_chunks_lbfgs_2 = optimize(x -> solar_trip_test(x, track), lower, upper, initial_x, Fminbox(LBFGS()))
# minimized_inputs_chunks_lbfgs_2 = Optim.minimizer(result_chunks_lbfgs_2)
# minimized_result_chunks_lbfgs_2 = Optim.minimum(result_chunks_lbfgs_2)
# result = optimize(f, [30.0])
# result = optimize(f, [30.0], LBFGS())
# result = optimize(f, [30.0], GradientDescent())

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

# https://juliahub.com/ui/Packages/LightGraphs/Xm08G/1.3.5?t=0 for graphs


end # module
