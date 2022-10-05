module SolarStrategy

using DataFrames
# also consider using JuliaDB and Query
using CSV
using Plots # default
using PlotlyBase
using TimeZones
using Dates
using Optim
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

# TODO: track preprocessing


function solar_trip_calculation(input_speed::Vector{Float64}, track)
    # input speed in m/s

    # calculating time needed to spend to travel across distance
    time_df = calculate_travel_time_datetime(input_speed, track)

    #### calculcations
    # mechanical calculations are now in separate file
    mechanical_power = mechanical_power_calculation(input_speed, track.slope, track.diff_distance)

    # electical losses
    electrical_power = electrical_power_calculation(track.diff_distance, input_speed)
    # converting mechanical work to elecctrical power and then power use
    power_use = calculate_power_use(mechanical_power, electrical_power)
    power_use_accumulated_wt_h = calculate_power_use_accumulated(mechanical_power, electrical_power)

    # get solar energy income
    intensity = solar_radiation_pvedication_time(time_df, track)
    plot(track.distance, intensity, title="Solar intensity")
    solar_power = solar_power_income(time_df, track, input_speed)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power)
    # TODO: night charging with additional solar panels

    #### plotting
    plot(track.distance, power_use, title="Power spent on toute")
    plot(track.distance, solar_power, title="Power gained on the route")

    plot(track.distance, power_use_accumulated_wt_h, title="Power spent on the route, accumulated")
    plot(track.distance, solar_power_accumulated, title="Power gained on the route, accumulated")

    plot(track.distance, solar_power_accumulated - power_use_accumulated_wt_h, title="Power balance w/o battery")
    battery_capacity = 5100 # wt
    energy_in_system = battery_capacity .+ solar_power_accumulated .- power_use_accumulated_wt_h
    plot(track.distance, energy_in_system, title="Power balance with battery")

    # TODO: calculate night charging - do it later since it is not critical as of right now
    # TODO: block overcharging - cost function?
    # at first do the black-box optimization, then gradient one
    # will start with Optim
    # TODO: find an optimal single speed - make a loss function and start optimization process
    time_seconds = calculate_travel_time_seconds(input_speed, track)
    cost = last(time_seconds) + 10 * abs(minimum(energy_in_system)) + sum(input_speed[input_speed .< 0.0])
    # TODO: find an optimal speed vector
    return cost
end

# consider using Union for types for multiple dispatch
# https://stackoverflow.com/questions/65094714/efficient-way-to-implement-multiple-dispatch-for-many-similar-functions 
# also consider rewriting core functions so they will work in .func mode

function solar_trip_wrapper(speed::Float64, track)
    # input speeds in km/h
    # as of right now one speed for all track parts
    input_speed = convert_single_kmh_speed_to_ms_vector(first(speed), length(track.distance))

    return solar_trip_calculation(input_speed, track)
end

function solar_trip_wrapper(speed::Vector{Float64}, track)
    # input in km/h
    return solar_trip_calculation(convert_kmh_to_ms(speed), track)
end

# function to test optimization with several big chunks to optimize
# everything under Number
function solar_trip_test(speeds::Vector{<:Number}, track)
# function solar_trip_test(speeds::Vector{Float64}, track)
    speed_ms = convert_kmh_to_ms(speeds)
    speed_vector = propagate_speeds(speed_ms, track)
    return solar_trip_calculation(speed_vector, track)
end


# regular optimization, Nelder-Mead, 9 seconds 
@time result = optimize(x -> solar_trip_wrapper(x, track), [30.0])
minimized_inputs = Optim.minimizer(result)
minimized_result = Optim.minimum(result)
# get few big chunks with the same speed, Nelder-Mead, ~210 seconds
@time result_chunks = optimize(x -> solar_trip_test(x, track), [41.0, 42.0, 43.0, 44.0, 45.0])
minimized_inputs_chunks = Optim.minimizer(result_chunks)
minimized_result_chunks = Optim.minimum(result_chunks)
# few big chunks, LBFGS, 645 seconds, negative speeds, consider revising or constraints
@time result_chunks_lbfgs = optimize(x -> solar_trip_test(x, track), [41.0, 42.0, 43.0, 44.0, 45.0], LBFGS())
minimized_inputs_chunks_lbfgs = Optim.minimizer(result_chunks_lbfgs)
minimized_result_chunks_lbfgs = Optim.minimum(result_chunks_lbfgs)
# result = optimize(f, [30.0])
# result = optimize(f, [30.0], LBFGS())
# result = optimize(f, [30.0], GradientDescent())


# TODO: refactor code, calculations in separate function, wrapper for optimization
# make an optimization with full vector

#### future use
# for optimization (overall list: https://www.juliaopt.org/packages/ ):
# https://github.com/robertfeldt/BlackBoxOptim.jl - looks like the case
# https://github.com/JuliaNLSolvers/Optim.jl - not sure if it handles derivative-free optimi
# https://github.com/JuliaOpt/Pajarito.jl - integer linear programming
# https://github.com/anriseth/MultiJuMP.jl - multicriterial optimiztion
# also http://julia.cookbook.tips/doku.php?id=optim , Nelder-Mead Simplex	NM

# https://juliahub.com/ui/Packages/LightGraphs/Xm08G/1.3.5?t=0 for graphs


end # module
