using DataFrames
using CSV
# using Plots # default
using Plotly # not as fast, but interactive
# using PlotlyJS # a lot of dependencies, slow loading

include("mechanical.jl")

#df_track = CSV.File(file) |> DataFrame!

track_csv = CSV.read("data_australia.csv")
#plot(track_csv.distance, track_csv.elevation)

# constants
panels_efficiency = 0.228 # 910/4000
engine_efficiency = 0.87
electrics_efficiency = 0.86
battery_efficiency = 0.98
power_onboard = 40 # Wt, 0.04kWt
mass = 390 # kg

friction_1 = 0.0023;
friction_2 = 0.000041; # total friction = friction_1 + friction_2*speed


ro = 1.18 # air density
battery_capacity = 5.100 # kWt*h
g = 9.8019 # at start: [41.2646201567207,-95.9244249307473,301.540649414063];
# 9.80147 at finish: [43.9660024736000,-121.345052439700,1229.07763671875]

panels_area = 4 # m^2
panels_area_charge = 6 # m^2
frontal_area = 1 # m^2
drag = 0.18

speed_kmh = 50
speed_ms = speed_kmh / 3.6

#=
overall concept of modelling:

speeds as an input
at first compute the energy loss, since it can be done in vectorized way
energy loss and time on each sector is calculated

now, since we know the times, calculate the energy income
this also can now be done in vectorized way


there are several possible models for energy income

start with stub, developo proper models later

=#


####### preprocessing of DataFrame
# meters from km
track_csv.distance = track_csv.distance * 1000
# create new dataframe with size of n-1 (and without spoiled sin column)
track = select(track_csv[2:size(track_csv,1),:], Not(:sin))
# # delete sin, since it is calculated incorrectly
# # delete!(track,:sin) # deprecated
# select!(track, Not(:sin))



track.diff_distance = diff(track_csv.distance)
track.diff_elevation = diff(track_csv.elevation)
# calculate sin and cos from catets throug atand
# calculate slope angle through arctan in degrees
# since diff is used, there is n-1 elements, so we have to add heading zero
# track.slope = [0 ; atand.(diff(track_csv.elevation)./diff(track_csv.distance))]
track.slope = atand.(track.diff_elevation./track.diff_distance)





####### calculations

# mechanical force = drag force + friction force + gravitational force
# newtons
mechanical_force = (
    drag * frontal_area * speed_ms^2 * ro / 2 .+
    .+ mass * g * (friction_1 + friction_2 * 4 * speed_ms) * cosd.(track.slope) .+
    .+ mass * g * sind.(track.slope)
    )

# mechanical power = mechanical force * distance delta / engine efficiency
# watts * s
mechanical_power = (
    mechanical_force .* track.diff_distance / (engine_efficiency)
    )

electrical_power = power_onboard * track.diff_distance / speed_ms

power_use = mechanical_power .+ electrical_power
power_use_accumulated = cumsum(power_use)
power_use_accumulated_wt_h = power_use_accumulated / 3600 / 1000

print(length(mechanical_work(speed_ms, track.slope, track.diff_distance)))














# for optimization (overall list: https://www.juliaopt.org/packages/ ):
# https://github.com/robertfeldt/BlackBoxOptim.jl - looks like the case
# https://github.com/JuliaNLSolvers/Optim.jl - not sure if it handles derivative-free optimi
# https://github.com/JuliaOpt/Pajarito.jl - integer linear programming
# https://github.com/anriseth/MultiJuMP.jl - multicriterial optimiztion
# also http://julia.cookbook.tips/doku.php?id=optim , Nelder-Mead Simplex	NM
