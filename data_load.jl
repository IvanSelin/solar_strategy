using DataFrames
# also consider using JuliaDB and Query
using CSV
# using Plots # default
using Plotly # not as fast, but interactive
# using PlotlyJS # a lot of dependencies, slow loading

include("mechanical.jl")
include("time.jl")


track_csv = CSV.read("data_australia.csv", DataFrame)
#plot(track_csv.distance, track_csv.elevation)


#### constants
# panels_efficiency = 0.228 # 910/4000
# electrics_efficiency = 0.86
# battery_efficiency = 0.98
power_onboard = 40 # Wt, 0.04kWt
# battery_capacity = 5.100 # kWt*h
# panels_area = 4 # m^2
# panels_area_charge = 6 # m^2

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


####### preprocessing of DataFrame
# meters from km
track_csv.distance = track_csv.distance * 1000
# create new dataframe with size of n-1 (and without spoiled sin column)
# sin column is ignored, since it is incomplete
# btw, it's better to store the slope angle
# track is for storing diff data, there are n-1 elements
track = select(track_csv[2:size(track_csv,1),:], Not(:sin))
track.diff_distance = diff(track_csv.distance)
track.diff_elevation = diff(track_csv.elevation)
track.slope = atand.(track.diff_elevation./track.diff_distance)


####### input
speed_kmh = 50
speed_ms = speed_kmh / 3.6

#### time manipulation
# base time, seconds driven from start of the race
time_s = track.distance ./ speed_ms
# coverting the journey time to real time
time_df = travel_time_to_real_time(time_s)

#### calculcations
# mechanical calculations are now in separate file
mechanical_power = mechanical_power_calculation(speed_ms, track.slope, track.diff_distance)

# electical losses
electrical_power = power_onboard * track.diff_distance / speed_ms
# converting mechanical work to elecctrical power and then power use
power_use = mechanical_power .+ electrical_power
power_use_accumulated = cumsum(power_use)
power_use_accumulated_wt_h = power_use_accumulated / 3600


#### plotting
plot(track.distance,power_use_accumulated_wt_h)



#### future use
# for optimization (overall list: https://www.juliaopt.org/packages/ ):
# https://github.com/robertfeldt/BlackBoxOptim.jl - looks like the case
# https://github.com/JuliaNLSolvers/Optim.jl - not sure if it handles derivative-free optimi
# https://github.com/JuliaOpt/Pajarito.jl - integer linear programming
# https://github.com/anriseth/MultiJuMP.jl - multicriterial optimiztion
# also http://julia.cookbook.tips/doku.php?id=optim , Nelder-Mead Simplex	NM
