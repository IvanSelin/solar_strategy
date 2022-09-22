module SolarStrategy

using DataFrames
# also consider using JuliaDB and Query
using CSV
using Plots # default
using PlotlyBase
using TimeZones
using Dates
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


#### constants
# panels_efficiency = 0.228 # 910/4000
# electrics_efficiency = 0.86
# battery_efficiency = 0.98

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

# preparing the track data
track = get_track_data("data/data_australia.csv")

# TODO: track preprocessing

# input speeds are here
# as of right now one speed for all track parts
input_speed = convert_single_kmh_speed_to_ms_vector(50, length(track.distance)) # kmh
# in the end it should be a vector
# calculating time needed to spend to travel across distance
time_df = calculate_travel_time(input_speed, track)

#### calculcations
# mechanical calculations are now in separate file
mechanical_power = mechanical_power_calculation(input_speed, track.slope, track.diff_distance)

# electical losses
electrical_power = electrical_power_calculation(track.diff_distance, input_speed)
# converting mechanical work to elecctrical power and then power use
power_use_accumulated_wt_h = calculate_power_use_accumulated(mechanical_power, electrical_power)


#### plotting
plot(track.distance,power_use_accumulated_wt_h)

# for development purposes, temporary code
data_df = generate_year_time_dataframe(100000)
data_df_with_solar = solar_radiation_pvedication_time(data_df)

# TODO: calculate proper power income

#### future use
# for optimization (overall list: https://www.juliaopt.org/packages/ ):
# https://github.com/robertfeldt/BlackBoxOptim.jl - looks like the case
# https://github.com/JuliaNLSolvers/Optim.jl - not sure if it handles derivative-free optimi
# https://github.com/JuliaOpt/Pajarito.jl - integer linear programming
# https://github.com/anriseth/MultiJuMP.jl - multicriterial optimiztion
# also http://julia.cookbook.tips/doku.php?id=optim , Nelder-Mead Simplex	NM

# https://juliahub.com/ui/Packages/LightGraphs/Xm08G/1.3.5?t=0 for graphs


end # module
