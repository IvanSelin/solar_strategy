### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fb9bd7c0-4f3c-11ed-094a-35559b7aedad
begin
	using DataFrames
	using CSV
	using Plots # default
	using PlotlyBase
	# using PlotlySave
	using TimeZones
	using Dates
	using Optim
	using LineSearches
	using PlutoUI
	using Peaks
	using ReverseDiff
	# using Alert
	# selecting a Plots backend
	# plotly(ticks=:native)
end

# ╔═╡ 7cc7b885-f841-4bcc-a82c-2f947c74de22
begin
	include("src//energy_draw.jl")
	include("src//time.jl")
	include("src//solar_radiation.jl")
	include("src//track.jl")
	include("src//utils.jl")
end

# ╔═╡ 90cd5c2a-223c-4345-8f78-498973c0fd46
PlutoUI.TableOfContents()

# ╔═╡ a7781b82-48db-4954-8e79-ab8e7864ed69
md"""
# Loading the track data
"""

# ╔═╡ 44f352a0-0c5b-41f6-a21d-56f10edae4c9
track = get_track_data("data/data_australia.csv")

# ╔═╡ 9d915c29-314b-47f9-9e4c-898ebd28f88a
# plotly(ticks=:native)

# ╔═╡ 92e47c4c-2e70-4850-af86-a997fcce5587
plot(track.distance, track.altitude, title="Altitude (m)", xlabel="Distance (m)", ylabel="Altitude (m)")

# ╔═╡ 9f2624cc-398b-42d2-bf7f-527285af6dc3
md""" # Defining the functions """

# ╔═╡ c3d41246-37cb-45f6-ac98-f64d9158087a
md""" # Playground """

# ╔═╡ 4c91701c-f90e-48cf-bb9d-d925baa85667
md""" ## Single speed graphs """

# ╔═╡ 9977c5d5-4872-41d4-8e99-1e4034493e4d
# speed in kmh
@bind speed NumberField(0.0:0.1:100.0, default=40.0)
#speed = 43.0;

# ╔═╡ 54a5bf51-1723-4ce1-8f6b-1fd199c991b5
md""" ## Optimization with different methods"""

# ╔═╡ c7434aee-3089-4098-8c9d-2d13c9df5ee9
@bind speed_chunks NumberField(0.0:0.1:100.0, default=40.0)

# ╔═╡ 8794ae20-fe98-47ab-bd80-681c09edb7d1
@bind number_of_chunks confirm(NumberField(1:50000, default=1))

# ╔═╡ 7380a326-1e9e-437c-9cb7-3aa2b54b8ec5
@bind lbfgs_m NumberField(1:50000, default=10)

# ╔═╡ 596e51be-969e-4fe6-8170-22c1aac89aca
md""" Selecting the line search and other parameters """

# ╔═╡ eeb85fd3-6721-4e0f-aed9-73731960ac35
# ls = LineSearches.HagerZhang()
ls = LineSearches.BackTracking();

# ╔═╡ 09bd67a1-a0f9-45f3-9839-2e9e592f01de
iterations_num = 10000

# ╔═╡ 264d9d64-6a6d-4e84-9af1-5795bd5bf829
begin
	lower_bound_c = fill(0.0, number_of_chunks)
	upper_bound_c = fill(100.0, number_of_chunks)
	tdc_0 = TwiceDifferentiableConstraints(lower_bound_c, upper_bound_c)
end

# ╔═╡ 655ff7ad-d7fc-47d4-bd22-0bb2c4b63cd5
@md_str " # Towards the recursive optimization! "

# ╔═╡ da6bb535-9601-4cf9-ba64-f8bfc68f3e5d
@md_str " ## The code itself"

# ╔═╡ 4f9da3c2-36e1-405a-9c83-edce6f9a30af
@md_str " ### Splitting the track"

# ╔═╡ e7445b07-20e8-404e-8d40-eccb9eb6ddc2
@md_str " ### Time calculation"

# ╔═╡ d240b89b-2d6c-4379-8132-1a47ec0ef6f3
time_size = 1000

# ╔═╡ 2b239f4b-7149-4771-8498-35b0314fb928
time_s_temp = fill(100, time_size)

# ╔═╡ 899c3a73-28ef-46e3-b08b-d2335d91f52d
@md_str "## That's where error happens!

Starting on a 8:41 i should have it still running until 16:00, and not from 16:41 to some other time"

# ╔═╡ 55157e96-b0b0-40fd-9db1-52edf69b4001
DateTime(2022,7,2,8,41,20) - DateTime(2022,7,2,0,0,0)

# ╔═╡ d2b833e2-fd0a-4ae9-8f98-eb03f658339e
Dates.second(DateTime(2022,7,2,8,41,20))

# ╔═╡ 848b59f6-d3f4-49c6-af22-7a40f5449071
(Dates.Time(DateTime(2022,7,2,8,41,20)) - Time(0,0,0)) / 10^9

# ╔═╡ cec4a60b-07f6-429a-8c15-6211ba3ded74
Dates.Date(DateTime(2022,7,2,8,41,20))

# ╔═╡ f3e06f95-0dd8-4c48-bbdc-fdac37be6548
Dates.Second.(cumsum(time_s_temp))

# ╔═╡ f0206c0b-e35f-499b-8aa0-32c69541b2e3
@md_str " ### Trip calculation"

# ╔═╡ ff9dad0d-4a00-439c-b704-486e395e5997
@md_str " ### Cost function "

# ╔═╡ 08f34cff-83b5-4e0b-a789-3c180c8d5cea
@md_str " ### Setting speeds by index"

# ╔═╡ f243f6d4-c1b8-4284-abed-d74ae47a7af5
function set_speeds_grad(speeds, track, divide_at)
	# output_speeds = fill(last(speeds), size(track.distance, 1))
	output_speeds = fill(speeds[1], divide_at[1])
	for i=2:size(divide_at,1)
		output_speeds = vcat(output_speeds, fill(speeds[i], divide_at[i] - divide_at[i-1]))
	end
	return output_speeds
end

# ╔═╡ 6baaa312-94f8-4fc3-923f-30033a4a75d4
@md_str " ### Wrapper"

# ╔═╡ cb7bf1ac-4abd-4c08-8893-e367cda83f8f
@md_str " ### Hierarchical optimization"

# ╔═╡ 9a6f42a2-7a9f-41ef-8ff2-b325a5971e42
size(track.distance, 1)

# ╔═╡ e9e72767-4fd3-4e9f-97bc-01c1f35ec916
# TODO: test hierarchical_optimization(speed, track, chunks_amount, start_energy, finish_energy, start_datetime, iteration) on some small_track

# ╔═╡ c1c47be6-b0a4-41bf-a284-26f193792748
md""" ## Experiment set-up

1. Take a short part of the track (50-100 pcs)
2. Perform a regular optimization
3. Try to run a hierarchical optimization
4. Compare results"""

# ╔═╡ eb193e5d-65ba-46b1-be15-5c399abad44b
@md_str " ### Short track "

# ╔═╡ 64ec134c-e5cf-4c97-869f-d39b91e2599e
size(keep_extremum_only_peaks(track),1)

# ╔═╡ ca583c57-d93c-4bb3-8cd3-b92184226a5f
# short_track = track[1:track_size, :]
short_track = keep_extremum_only_peaks(track)#[1:track_size, :]
# short_track = copy(track)

# ╔═╡ e75a6ae6-a09c-4c21-a421-d0124dd355c6
# track_size = 13148
track_size = size(short_track.distance, 1)
# track_size = 5000

# ╔═╡ 19d51114-43d4-4c38-9a3d-55ec982f7c56
distance_perc = last(short_track).distance / last(track.distance)

# ╔═╡ 2de331ba-bd7d-49bd-80a2-60270c769a7c
proposed_start_energy = distance_perc * 5100

# ╔═╡ 724deff0-8265-4ca9-aa3e-2dfdd6f4d293
@md_str " ### Regular optimization on short track (w/o propagation)"

# ╔═╡ 56f144c3-43ed-46d5-9e55-fd800bd24e9e
# start_energy_short = 150.0 # W*h
start_energy_short = proposed_start_energy

# ╔═╡ 7699d86d-4494-494a-b367-7c3f13fc0022
initial_speed = 40.0

# ╔═╡ f230c95a-d02a-436b-b426-364fe112cbba
@md_str "First run without optimization to check if everything is OK"

# ╔═╡ f6ed6270-5871-4dec-b2a3-5fa78caff34b
# @time result_short = optimize(td_short, fill(initial_speed, track_size), Newton(; linesearch = ls),
# 	Optim.Options(
# 		x_tol = 1e-6,
# 		f_tol = 1e-8,
# 		g_tol = 1e-6
# 	)
# )

# ╔═╡ 8d2d34de-64db-49a8-b874-cb5af3bcf6cd
minimized_inputs_short = Optim.minimizer(result_short)

# ╔═╡ c1a6f532-e95b-4fc3-b419-255038ee3589
begin
	plot(short_track.distance, short_track.altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance", right_margin = 15Plots.mm)
	plot!(twinx(), short_track.distance, minimized_inputs_short, color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance")
end

# ╔═╡ 100ab330-4b6c-4033-92ee-9d19b2f7d2e4
@md_str " ### Launching a hierarchical optimization "

# ╔═╡ cc75dbac-9f9c-4465-a08b-e543d82021fd
chunks_amount_hierarchical = 10

# ╔═╡ 0c127ac8-812e-4a99-9e0e-32b3fd316c26
start_datetime_hierarchical = DateTime(2022, 7, 1, 0, 0, 0)

# ╔═╡ a8b6ce9d-5507-4acb-954f-b43d702e1060
@md_str " ### Result speed graph"

# ╔═╡ f0d3e885-68c3-424e-8a50-d3d981bef295
@md_str " ### Result energy graph "

# ╔═╡ e1c90043-eae2-40fd-afc4-e5cd0ab4a540
@md_str "### Speed graph for 1st iteration"

# ╔═╡ 9a73fa1f-ac16-420d-8c98-91d1aeba773a
@md_str " #### And for selected bad conditioned part on 2nd iteration"

# ╔═╡ 55a9fc54-39b7-4c40-b06e-eba6ba260c50
@md_str " #### Negative energy income? WTF?!"

# ╔═╡ bf7b5b54-2f2a-4f5c-a877-1fd8f53510b0
# solar_income_2_no_speed = solar_power_income(time_2, short_track[start_index_2:end_index_2,:])

# ╔═╡ 25f0fcbb-382e-4254-b058-52af29f212f9
@md_str "## Time problems investigation"

# ╔═╡ 7d514799-a9b2-482f-aa16-c7d4a9fbfe78
@md_str " #### Why do i get negative solar income?

I need to investigate the source of the problem"

# ╔═╡ 32160090-9fef-4216-a160-76a0f0af0f0f
@md_str " #### Why does everything starts from 16:45, and not from 8:45? 

Is it the problem in the datetime calculation shift?

Check difference between UTC and LST times"

# ╔═╡ 8d6a6e0f-044b-4c6e-8432-b5da6817d019
@md_str " LST looks fine

https://www.timeanddate.com/worldclock/converter.html?iso=20220701T080000&p1=72&p2=1440

For UTC 01-07-2022T08:00 it is 17:30 in Darwin and LST is ~16:38

(caclulated by https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-time with parameters: 130.83, 9.5, 182, 17, 30)

"

# ╔═╡ 5d4ac822-380e-477d-b8a3-7d6f39b0d9ec
@md_str "#### Maybe we are riding during the wrong hours?

Since model input is UTC time, not LST time"

# ╔═╡ adc5d37e-5f66-4d11-8ba3-1426bb823230
@md_str " calculate travel time datetime produces only utc_time from speeds, track and start datetime"

# ╔═╡ 65f69edd-24a7-4a21-91b3-2e86c17df1a1
time_s_test = collect(1:100:100000)

# ╔═╡ f7b832ee-c8d5-4144-bd18-a68ed34c3e6a
@md_str "#### Check how solar income works for constant time change in fixed place"

# ╔═╡ 2e5f895b-c8b4-4bbb-8e69-88c554022be6
start_test = DateTime(2022,1,7,0,0,0)

# ╔═╡ 8f29c327-a738-4599-be47-4ff868b303e6
@md_str "where are hours 0-8 and 16-24?"

# ╔═╡ 51765acb-8459-426c-9c3a-773782f9789f
@md_str "try to simulate power income for a fixed point over a period of time"

# ╔═╡ c8a9d213-b294-4c46-9af3-3682ea171766
@md_str "data prep"

# ╔═╡ 981ee4ea-1a83-44de-81da-6e472496aa7b
track[1:5,:]

# ╔═╡ c8941591-bacb-469e-8d7f-cf15d32169d3
begin
	track_test_df = DataFrame(track[1,:])
	for i=1:999
		track_test_df = vcat(track_test_df, DataFrame(track[1,:]))
	end
end

# ╔═╡ 5bc52c1f-bf9f-447c-a7d8-795d374ca4f4
track_test_df

# ╔═╡ f4c411a4-1c12-42b6-aa9b-030d7afdccc4


# ╔═╡ 6187da23-54e1-475d-9afd-73cd90600088
speed_spi_test = fill(2.4, 1000)

# ╔═╡ e7eafaf3-35be-43ab-a8d5-30bff227664b
@md_str " #### So, maximum solar activity is at 12:00 utc_time 

Which is, on one hand, is wrong (UTC12+8:30 is 20:30 when sun is set)

But, on the other hand, it means that code is working based on utc time, like it is a local time

Which means that most likely there is an error with +8h in datetime"

# ╔═╡ 2f20704a-f2e0-4843-b458-6c3208758ab6


# ╔═╡ aa00efad-65b3-449d-9085-da6cf66c57c6


# ╔═╡ b7e2d987-463b-4d51-ba98-6cdbd3fbb01c
@md_str " #### Trying out starting speed from above iteration "

# ╔═╡ 9f247245-a14d-42ce-97e7-a60911988d9c
@md_str " #### Something is wrong with iterations/optimization design 

And again, negative energy income"

# ╔═╡ 39b4c2a0-fb33-4820-87f3-fa575ebb4e30
@md_str " ### Comparing the results "

# ╔═╡ e983aeb2-38c2-4bc4-af61-8af08f2347f5
# TODO: make a wrapper function for calling experiments with different track size and chunks amount. And which will also diplay plots

# ╔═╡ 3c7d8d1d-f975-4743-99ca-263ac269fcab
# TODO: somehow restrict huge speed changes

# ╔═╡ 97bbd950-40f6-4026-afc0-eecef2d0c784
# TODO: run hierarchical optim line-by-line

# ╔═╡ 16d51785-7cc0-4943-ba5a-0f0241315cfd
# TODO: add penalty for big speed difference

# ╔═╡ bb69fdb6-b6b8-42ab-9a27-0f926067c1a6
# TODO: maybe speed spikes are because of boundaries in hierarchical optim?

# ╔═╡ 569ab1f6-c5aa-48d8-88ae-1884ca22518a
# TODO: try to simulate bad conditioned hierarchical pieces where finish_energy not met and speeds differ a lot

# ╔═╡ 7107b365-28d2-4ec7-bea7-a93b5c89a9d8
# TODO: run hierarchical with plots depth-by-depth

# ╔═╡ 710b23ba-7648-4795-96e4-8141766479d5
# for some reason energy in control run differs from energies inside hierarchical opt#
# maybe because not every time finish energy is the same as it should?

# ╔═╡ f744cc32-907f-40b0-8586-20af362b2dc9
# build a parameter plane to see what is really happening
# make a 2-part track
# and select different speeds and see resulted cost
# make a surf plot

# ╔═╡ 873458fc-f7bd-4154-838d-0f65461f0178
track_2_pcs = short_track[25:26,:]

# ╔═╡ 4dd7879c-147e-41d6-8527-83dbac04567b


# ╔═╡ 14a03646-a5b6-4381-a3ed-39a6eb85c113
test_res = 5.

# ╔═╡ 4261573f-2207-45e2-b117-49c64079263a
# ftape_vec = ReverseDiff.GradientTape(f_test_plane_3d_grad, rand(3))
# not differiantiable - see last point in https://juliadiff.org/ReverseDiff.jl/limits/
# revert to simpler methods?

# ╔═╡ fb78a338-c2ab-499e-bddd-bf36a15ea8cc
# function cost_calc_3d_1(speed1, speed2)
# 	return solar_partial_trip_test_wrapper([30., speed1, speed2], track[25:27,:], [1, 2,3], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
# end

# ╔═╡ 56b71cfd-b90b-4d2c-9307-2b3b97a30c6f
# function cost_calc_3d_2(speed1, speed2)
# 	return solar_partial_trip_test_wrapper([speed1, 30., speed2], track[25:27,:], [1, 2,3], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
# end

# ╔═╡ c3bd2406-abc7-4e48-b8af-a231c9cda891
# function cost_calc_3d_3(speed1, speed2)
# 	return solar_partial_trip_test_wrapper([speed1, speed2, 30.], track[25:27,:], [1, 2,3], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
# end

# ╔═╡ 50f0d299-5b25-452d-9a95-adb121d6abd4
Plots.surface(5:5:80, 5:5:80, cost_calc_3d_1)

# ╔═╡ 4ffa67f3-c4bd-47e6-885a-28cab178dd61
Plots.surface(5:5:80, 5:5:80, cost_calc_3d_2)

# ╔═╡ ca232f87-1c5e-4cdf-8ed8-aef9e5a34db2
Plots.surface(5:5:80, 5:5:80, cost_calc_3d_3)

# ╔═╡ fa0a86ad-a9e7-4d53-b7dd-7866f85c000f
# make a penalty for different speeds

# ╔═╡ 167ec4cf-c863-4a85-b429-02571fad0f30
# it is better to use reverse diff

# ╔═╡ 8193ba21-d8d2-40e0-af7b-3028d26afafc
# Functions like f which map a vector to a scalar are the best case for reverse-mode automatic differentiation, but ForwardDiff may still be a good choice if x is not too large, as it is much simpler. The best case for forward-mode differentiation is a function which maps a scalar to a vector, like this g:
# from https://github.com/JuliaDiff/ForwardDiff.jl 

# ╔═╡ 98eef687-eff7-4872-8541-85baeb577b77
# so, trying out https://github.com/JuliaDiff/ReverseDiff.jl

# ╔═╡ 6e64c5d0-c9cf-4c4b-a9df-1f28225b7eff
# example here: https://github.com/JuliaDiff/ReverseDiff.jl/blob/master/examples/gradient.jl 

# ╔═╡ e2a35442-6974-4607-bbc0-b12ca32b4594
@md_str " ## Trying out solar calc without broadcasting"

# ╔═╡ b2707113-fb5b-4600-bd48-4be882759436
speed_vec_no_b = fill(40.0, size(short_track.distance, 1) )

# ╔═╡ 10b601b5-b24e-44e9-8d1a-c14aeba3dcf7
short_track

# ╔═╡ b88cc8dd-39d6-4879-98a2-cccaf6e000b7
first(short_track)

# ╔═╡ 2ab21d7e-7416-4057-a53e-8ebde10295e6
first(short_track).latitude

# ╔═╡ 84ce8c14-0f67-46a6-9386-818c8003818f
first(short_track).lstm = 12

# ╔═╡ 63775bee-b22e-4f99-a16c-d74ceaa91e65
@md_str "### Nice! It's time to understand how to make it work with arrays"

# ╔═╡ cecb4231-6cdf-42dd-81fc-7fb304e4ffc0
size(short_track)

# ╔═╡ b10c3e4b-edbc-474e-8c79-f37d9726439f
size(speed_vec_no_b)

# ╔═╡ a57b759e-9c4b-49c4-9439-4b842131f325
@md_str "Should I pass not dataframe, but separate arrays?"

# ╔═╡ e662491c-5165-4d3f-926f-fac8549e94e0
@md_str "0.007323 seconds (15 allocations: 103.297 KiB)"

# ╔═╡ 564bb277-404c-413f-9e5b-1d56df39040f
@md_str "0.006838 seconds (14 allocations: 103.234 KiB)"

# ╔═╡ 5b04879a-5029-4439-b4df-c1449fe27834
@md_str "0.013377 seconds (323 allocations: 5.050 MiB)"

# ╔═╡ 0d4cccc7-56bf-45c2-8529-16adf8eaecc7
@md_str "### Regular version makes a lot more allocations!

Time to switch to less-allocative code"

# ╔═╡ 1aff4c3c-f4fa-4cda-9810-f91514dc81a1
@md_str "trying out in-place variant"

# ╔═╡ f336bb29-2ac3-41e0-9930-623e8fea3af5
@md_str " # Trying out faster function in optimization"

# ╔═╡ 5051d08b-4a09-4e59-ae81-2caf9e555be4
@md_str "# Assembling the alloc Frankenstein "

# ╔═╡ e5b2c108-25a9-4035-bc60-5c5ca3eb54f1
@md_str " About 35-40 seconds on second run 

2-3x times faster!"

# ╔═╡ 559f597f-0d11-4212-b1e7-9ccd22661ae0
@md_str "## Time to transfer solar power income expanded alloc into the main code

Along with solar trip calculation bounds alloc, mechanical power calculation alloc, solar partial trip wrapper alloc and hierarchical optimization alloc"

# ╔═╡ 3f04f5b8-1670-4155-865e-3212b3973364
@md_str "# Trying out reverse diff yet again"

# ╔═╡ 35afbb71-6862-4c0d-9582-2e80079f39c8
+.([1,2,3],[4,5,6])

# ╔═╡ 84d06bf5-a160-4dcc-af3b-494b105d7e22
result_reverse = 1.0e8
# just some big number

# ╔═╡ 501e8567-f7c2-47e6-adaf-d3a463012912
inputs_reverse = [20., 30., 40.]

# ╔═╡ cf9aa152-1095-4d38-8641-942f0de4b860
# result_reverse2 = zeros(size(inputs_reverse, 1))
result_reverse2 = zeros(chunks_amount_hierarchical)

# ╔═╡ d6963ed3-39c3-4c60-9cfa-0d2f1fdd3a99
inputs_reverse2 = fill(40., chunks_amount_hierarchical)

# ╔═╡ 28a9dbdb-909a-43f3-a7bd-ace577293043
@md_str "### For some reason hessian is not calculating

Consider changing solar trip calculation bounds alloc for some other function"

# ╔═╡ 2606492a-89fd-4c2a-8708-131784b0a053
@inline Base.:+(x::ReverseDiff.TrackedArray{V,D}, y::ReverseDiff.TrackedArray{V,D}) where {V,D} = record_plus(x, y, D)

# ╔═╡ 09b71f55-2e31-4a7f-a6ff-f4a3ff38e4b9
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

# ╔═╡ 516e68a9-fb52-475c-b925-bba877341499
calculate_split_indexes(1,1)

# ╔═╡ 9384af58-d6ad-4abd-8df0-f4e0d5649e68
ind = calculate_split_indexes(size(track.distance, 1), 3)

# ╔═╡ 718fc4a4-8f46-4dff-a9a6-6584c134308a
split_indexes_reverse = calculate_split_indexes(size(short_track, 1), chunks_amount_hierarchical)

# ╔═╡ a5551c27-3f5c-4ea3-8c4e-9d4873c88751
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

# ╔═╡ 7ed6b68e-2cb3-4bd0-8c98-55ee42349d93
tracks = split_track_by_indexes(track, ind)

# ╔═╡ 9db98630-5796-4763-92e8-d04a4ad3845a
tracks[2].distance

# ╔═╡ b6a95a4d-6f5d-4fbc-87ed-622768e8d47a
function travel_time_to_datetime(time_s, start_datetime)
	# adjust to timezone, this is for UTC!
	utc_fix = 9 # since darwin is UTC+9:30
	daily_start_hour_time = 8
    daily_finish_hour_time = 16

    start_time_seconds = daily_start_hour_time * 60 * 60
    finish_time_seconds = daily_finish_hour_time * 60 *60
    seconds_in_a_day = 24 * 60 * 60

    # adjust travel time so it happend only between daily start hour time and daily finish hour time
    # TODO: think of in-place operations to reduce memory consumption
    day_length = finish_time_seconds - start_time_seconds
    day = time_s .÷ day_length .+ 1
    time_s .+= start_time_seconds .* day .+
    ( seconds_in_a_day .- finish_time_seconds) .* (day .- 1)

    # create a DataFrame for time information
    # adds seconds for 

	# creating new start datetime with proper timezone
	# new_start_datetime = start_datetime - Dates.Hour(utc_fix)
    return start_datetime .+ Dates.Millisecond.(round.(time_s .* 1000))
	# return new_start_datetime .+ Dates.Millisecond.(round.(time_s .* 1000))
end

# ╔═╡ 7be923ff-1be1-4496-a571-b3dd68a03ebc
travel_time_to_datetime([5, 15, 25], DateTime(2022,7,1,0,0,0))

# ╔═╡ b57ed73a-a69e-4944-b6f2-2e6bcafc2442
tt_original = travel_time_to_datetime(cumsum(time_s_temp), DateTime(2022,7,1,0,0,0))

# ╔═╡ 07cd4d1e-a293-4aa4-85fc-43cc26be9f29
plot(tt_original)

# ╔═╡ 9f924c63-9558-4915-ab1a-f8a383795906
plot(tt_original, 1:time_size)

# ╔═╡ 7779fe61-ca38-4bb6-aa6e-e7bc84870fe4
travel_time_to_datetime(cumsum(time_s_temp),DateTime(2022,7,2,8,41,20))

# ╔═╡ 249b0a50-0713-474f-a2ce-f7c012522318
plot(travel_time_to_datetime(cumsum(time_s_temp),DateTime(2022,7,2,8,41,20)), 1:time_size)

# ╔═╡ 35285182-9a52-4e37-959a-09b41927d56f
Dates.hour(DateTime(2022,7,2,8,41,20))*3600 + Dates.minute(DateTime(2022,7,2,8,41,20))*60 + Dates.second(DateTime(2022,7,2,8,41,20))

# ╔═╡ 21a36ea3-d1d8-4a64-bb19-12cb85df5da7
function travel_time_to_datetime_new(time_s, start_datetime)
	start_date = Dates.Date(start_datetime)
	start_seconds = Dates.hour(start_datetime)*3600 + Dates.minute(start_datetime)*60 + Dates.second(start_datetime)
	time_s .+= start_seconds
	
	# adjust to timezone, this is for UTC!
	utc_fix = 9 # since darwin is UTC+9:30
	daily_start_hour_time = 8
    daily_finish_hour_time = 16

    start_time_seconds = daily_start_hour_time * 60 * 60
    finish_time_seconds = daily_finish_hour_time * 60 *60
    seconds_in_a_day = 24 * 60 * 60

	

    # adjust travel time so it happend only between daily start hour time and daily finish hour time
    # TODO: think of in-place operations to reduce memory consumption
    day_length = finish_time_seconds - start_time_seconds
	days_amount = last(time_s) ÷ day_length
	if first(time_s) < start_time_seconds
		δ = start_time_seconds - first(time_s)
		time_s .+= δ
	end
	for day=1:days_amount
		# @debug day*finish_time_seconds + (day-1)*seconds_in_a_day
		time_s[time_s .> finish_time_seconds + (day-1)*seconds_in_a_day] .+= (seconds_in_a_day - finish_time_seconds) + start_time_seconds
	end
	
	# @debug days_amount
	# day = time_s .÷ day_length
	# @debug day
	# time_s .+= start_time_seconds .* day .+ ( seconds_in_a_day .- finish_time_seconds) .* (day .- 1)
	# @debug time_s

    # create a DataFrame for time information
    # adds seconds for 

	# creating new start datetime with proper timezone
	# new_start_datetime = start_datetime - Dates.Hour(utc_fix)
    return Dates.DateTime(start_date) .+ Dates.Millisecond.(round.(time_s .* 1000))
	# return Dates.DateTime(start_date) .+ Dates.Second.(round.(time_s))
	# return new_start_datetime .+ Dates.Millisecond.(round.(time_s .* 1000))
end

# ╔═╡ 805f1520-28b0-45dd-80a5-470a33bacdee
plot(travel_time_to_datetime_new(cumsum(time_s_temp),DateTime(2022,7,2,9,41,20)), 1:time_size)

# ╔═╡ 294fa952-015d-4709-bb55-c17682ffe2fb
function calculate_travel_time_datetime(speed_vector, track, start_datetime)
    time_s = calculate_travel_time_seconds(speed_vector, track)
	time_utc = travel_time_to_datetime_new(time_s, start_datetime)
	return DataFrame(utc_time=time_utc, time_s=time_s)
end

# ╔═╡ 704c9c6a-8b13-4e82-92fd-edda72397320

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


# ╔═╡ afe9cc57-f93c-4780-a592-b3d2609162f2

function solar_trip_cost(input_speed, track)
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(input_speed, track)
    cost = last(time_s) + 10 * abs(minimum(energy_in_system)) + 100 * sum(input_speed[input_speed .< 0.0])
    return cost
end

# ╔═╡ e80aee02-b99e-44c9-b503-9443c431b0e6

# consider using Union for types for multiple dispatch
# https://stackoverflow.com/questions/65094714/efficient-way-to-implement-multiple-dispatch-for-many-similar-functions 
# also consider rewriting core functions so they will work in .func mode

function solar_trip_wrapper(speed::Float64, track)
    # input speeds in km/h
    # as of right now one speed for all track parts
    input_speed = convert_single_kmh_speed_to_ms_vector(first(speed), length(track.distance))
    return solar_trip_cost(input_speed, track)
end

# ╔═╡ ae4d73ea-bd18-472b-babb-7980598a4ce9
function solar_trip_wrapper(speed::Vector{Float64}, track)
    # input in km/h
    return solar_trip_cost(convert_kmh_to_ms(speed), track)
end

# ╔═╡ b9bdb969-2f65-47d8-b1f2-a9b7e411f1c1
# function to test optimization with several big chunks to optimize
# everything under Number
function solar_trip_chunks(speeds::Vector{<:Number}, track)
# function solar_trip_chunks(speeds::Vector{Float64}, track)
    speed_ms = convert_kmh_to_ms(speeds)
    speed_vector = propagate_speeds(speed_ms, track)
    return solar_trip_cost(speed_vector, track)
end

# ╔═╡ 3ca6f786-8add-4c46-b82a-30a570828d39
function f(x)
	return solar_trip_chunks(abs.(x), track)
end

# ╔═╡ e50a7ae9-a46e-41b0-8a10-d77e9ffa7b14
d = OnceDifferentiable(f, fill(speed_chunks, number_of_chunks), autodiff = :forward)

# ╔═╡ 5c9c006c-f814-4829-8c18-108546be870b
td = TwiceDifferentiable(f, fill(speed_chunks, number_of_chunks), autodiff = :forward)

# ╔═╡ 411e63ec-b83a-4e21-9535-5d0275381039
#@time result_chunks = optimize(x -> solar_trip_chunks(x, track), fill(speed_chunks, number_of_chunks), iterations=iterations_num, LBFGS())

# @time result_chunks = optimize(d, fill(speed_chunks, number_of_chunks), LBFGS(;m=lbfgs_m, linesearch = ls))

# @time result_chunks = optimize(td, fill(speed_chunks, number_of_chunks), Newton(; linesearch = ls),
# 	Optim.Options(
# 		x_tol = 1e-6,
# 		f_tol = 1e-8,
# 		g_tol = 1e-6
# 	)
# )

@time result_chunks = optimize(td, tdc_0, fill(speed_chunks, number_of_chunks), IPNewton(),
	Optim.Options(
		x_tol = 1e-6,
		f_tol = 1e-8,
		g_tol = 1e-6
	)
)
# change convergence criteria?

# ╔═╡ 21634b70-7b3a-44b2-b629-01664ce81acf
minimized_inputs_chunks = Optim.minimizer(result_chunks)

# ╔═╡ 4e4c4794-aa95-49e9-961b-ed7c4bb81442

function solar_trip_target_cost(input_speed::Vector{Float64}, target_energy::Float64, track,
    start_energy::Float64, finish_energy::Float64)
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(input_speed, track)
    cost = last(time_s) + 10 * abs(last(energy_in_system) - target) + 100 * sum(input_speed[input_speed .< 0.0])
    return cost
end

# ╔═╡ ca2f1ec3-4ed3-4b5d-bfcd-ab43de0d2abc
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


# ╔═╡ 36e15048-b615-4ba6-a32a-093922a49243

function show_results_wrapper(input::Number, track::DataFrame)
    input_speed = convert_single_kmh_speed_to_ms_vector(first(input), length(track.distance))
    return show_result_graphs(input_speed, track)
end

# ╔═╡ 1536f51f-0776-4d87-832d-e9b2b9cc36d6

function show_results_wrapper(inputs::Vector{Float64}, track::DataFrame)
    inputs_ms = convert_kmh_to_ms(inputs);
    speed_vector = propagate_speeds(inputs_ms, track);
    return show_result_graphs(speed_vector, track)
end

# ╔═╡ 4f286021-da6e-4037-80e3-a526e880b686

function minimize_single_speed(speed::Float64, track)
    # regular optimization, Nelder-Mead, 9 seconds 
    @time result = optimize(x -> solar_trip_wrapper(x, track), [speed])
    # alert();
    minimized_inputs = Optim.minimizer(result)
    minimized_result = Optim.minimum(result)
    show_results_wrapper(first(minimized_inputs), track);
end

# ╔═╡ 1ba728ba-b6aa-40c9-b5a2-906297bd4921
function minimize_5_chunks(track)
    # get few big chunks with the same speed, Nelder-Mead, ~210 seconds
    @time result_chunks = optimize(x -> solar_trip_chunks(x, track), [41.0, 42.0, 43.0, 44.0, 45.0])
    minimized_inputs_chunks = Optim.minimizer(result_chunks)
    minimized_result_chunks = Optim.minimum(result_chunks)
    show_results_wrapper(minimized_inputs_chunks, track);
end

# ╔═╡ 5674f163-5847-4e2b-ba3f-c95348e8d1d5

function minimize(speed::Vector{Float64}, track)
    @time result = optimize(x -> solar_trip_wrapper(x, track), speed)
    minimized_inputs = Optim.minimizer(result)
    show_results_wrapper(minimized_inputs, track);
end

# ╔═╡ 422a0d48-40fe-41eb-b214-e21d009c00b2
begin
	input_speeds = convert_single_kmh_speed_to_ms_vector(speed, length(track.distance))
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(input_speeds, track);
	print()
end

# ╔═╡ 6a5416d0-39c0-431b-8add-5dbf13a1bda0
plot(track.distance, [power_use solar_power energy_in_system zeros(size(track,1))],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)

# ╔═╡ bef6c840-8dc7-4839-b2ba-623c6c46c856
plot(time.utc_time, [power_use solar_power energy_in_system zeros(size(track,1))],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph in time",
	    xlabel="Time", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
	    color=[:blue :green :cyan :red]
	    # ,ylims=[-10000, 40000]
	    )

# ╔═╡ f30abea4-a1d3-40fe-8328-3a0e5ce8a0d9
travel_time_to_datetime(time_s, DateTime(2022,7,2,8,41,20))

# ╔═╡ 8a3b49a9-4472-4b4d-944c-6b3e92a47b9f
plot(travel_time_to_datetime(time_s, DateTime(2022,7,2,8,41,20)))

# ╔═╡ 5f3a7fcf-e261-4f64-a94c-57a12093e353
begin
	inputs_ms = convert_kmh_to_ms(minimized_inputs_chunks)
	speed_vector = propagate_speeds(inputs_ms, track)
	power_use_chunks, solar_power_chunks, energy_in_system_chunks, time_chunks, time_s_chunks = solar_trip_calculation(speed_vector, track)
	last(time_s_chunks)
end

# ╔═╡ de201868-7805-4f27-81b7-f4f8204eface
begin
	plot(track.distance, track.altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance", right_margin = 15Plots.mm)
	plot!(twinx(), track.distance, speed_vector * 3.6, color=:red, ylabel="speed (km/h)", ylim=[0, 60], label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance")
end

# ╔═╡ 96a68dec-d781-4fd6-8146-649434f60919
plot(track.distance, [power_use_chunks solar_power_chunks energy_in_system_chunks zeros(size(track,1))],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance)",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)

# ╔═╡ 77d82639-dd61-46e0-b6a0-7c7400a10453
begin
	plot(time.utc_time, [power_use_chunks solar_power_chunks energy_in_system_chunks zeros(size(track,1))],
		    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"],
		    xlabel="Time", ylabel="Energy (W*h)", lw=3, size=(1200, 500),
		    color=[:blue :green :cyan :red]
		    # ,ylims=[-10000, 40000]
		, legend = :topleft 
		, right_margin = 15Plots.mm
		, title = "Energy graph (time)"
		    )
	plot!(twinx(), time.utc_time, speed_vector * 3.6, color=:red, ylabel="speed (km/h)", ylim=[0, 60], label="speed (km/h)", ymirror = true,
	title = "Energy graph (time)")
end

# ╔═╡ 69b2fd33-615a-4bf3-8936-8c9ee0651ad1
@time mechanical_power_calculation(speed_vector, track.slope, track.diff_distance)

# ╔═╡ 79afe017-7206-4c35-b5c1-84a6e4cc3517
begin
	inputs_ms_short = convert_kmh_to_ms(minimized_inputs_short)
	power_use_short, solar_power_short, energy_in_system_short, time_short, time_s_short = solar_trip_calculation(inputs_ms_short, short_track, start_energy_short)
	last(time_s_short)
end

# ╔═╡ 32bad172-97e7-485e-b64a-3c139400b471
plot(short_track.distance, [power_use_short solar_power_short energy_in_system_short zeros(track_size)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Newton",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, #size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)

# ╔═╡ 6a9d745e-e8ee-499c-ba72-fc95570ce442
time_df_no_b = calculate_travel_time_datetime(
	speed_vec_no_b,
	short_track,
	DateTime(2022,7,1,0,0,0)
)

# ╔═╡ 3ee9dae1-bccc-41b3-9980-3649683dae3c
size(time_df_no_b)

# ╔═╡ 5a64f67e-f9f0-4942-a3f9-7c7e468e9225
@time vectorized_function_result = solar_power_income(time_df_no_b, short_track, speed_vec_no_b)

# ╔═╡ c2dd681e-c318-424d-9038-f66ab9317c52
@time calculate_travel_time_datetime(speed_vector, track, DateTime(2022,7,1,0,0,0))

# ╔═╡ 7fc0e776-c9cf-4736-b113-81c31ada99c9
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

# ╔═╡ 1cd3cbae-96dd-4af2-ab6f-eb2565e6d6e1
@time power_use_test_bounds, solar_power_test_bounds, energy_in_system_test_bounds, time_test_bounds, time_s_test_bounds = solar_trip_calculation_bounds(speed_vector, track, DateTime(2022,7,1,0,0,0), 5099.0);

# ╔═╡ b872074b-094c-4d44-8d90-a1cf4908e012
@time travel_time_to_datetime(time_s_test_bounds, DateTime(2022,7,1,0,0,0))

# ╔═╡ 2b9b2782-09d8-4cfa-99fa-9c2e921efe36
function solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) + 10 * abs(finish_energy - last(energy_in_system)) #+ 0 * sum(abs.(energy_in_system[energy_in_system .< 0.0]))
	
	# + 100 * sum(abs.(speed_vector[speed_vector .< 0.0])) + 100 * sum(abs.(speed_vector[speed_vector .> 100.0 / 3.6]))
	# + 100 * abs(minimum(energy_in_system) - finish_energy) 
	
	# + 100 * sum(speed_vector[speed_vector .> 200.0]) + 100 * sum(energy_in_system[energy_in_system .< 0.0])


	# cost = last(time_s) + 10 * abs(minimum(energy_in_system)) + 100 * sum(input_speed[input_speed .< 0.0])
	return cost
end

# ╔═╡ 119c4024-6e50-4d20-a0b2-734ec9bf515d
function set_speeds(speeds, track, divide_at)
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

# ╔═╡ 360727dd-c241-4609-85fb-5f40553c1d2b
function solar_partial_trip_wrapper(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds(speeds_ms, track, indexes)
	return solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
end

# ╔═╡ 48e146a3-861b-4005-acf1-e877dbd83a50
begin
	function f_test(speed)
		return solar_partial_trip_wrapper(speed, short_track[25:28,:], [1,2,3,4], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
	end
	td_test = TwiceDifferentiable(f_test, fill(26.90116336126291, 4); autodiff = :forward)
	lower_bound_test = fill(0.0, 4)
	upper_bound_test = fill(100.0, 4)
	tdc_test = TwiceDifferentiableConstraints(lower_bound_test, upper_bound_test)
	# line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, chunks_amount),
	    #Newton(; linesearch = line_search),
	# inputs_test = fill(26.90116336126291, 4)
	inputs_test = [10., 20., 30., 40.]
	@time result_test = optimize(td_test, tdc_test, inputs_test,
		IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10,
			allow_f_increases = true,
			successive_f_tol = 50,
			show_trace = true
	    )
	)
end

# ╔═╡ 337efda5-da6c-40c2-97a7-486ba07d2175
minimized_speeds_test = Optim.minimizer(result_test)

# ╔═╡ f0fdeb9c-8791-4119-acc6-a1ec9094c590
set_speeds([1,2,3,4,5], track, ind)

# ╔═╡ a6c2c405-08e9-43ba-843a-d4bd85ded0c8
plot(set_speeds([1,2,3,4,5], track, ind))

# ╔═╡ bb76116a-5399-4d95-9e24-a158cd438619
begin
	speeds_1 = [43.564454419763095, 43.228254912151826, 43.60447581734879, 43.83531436356016, 43.987136764852366, 44.1190050188417, 44.271435203396315, 44.458800898738986, 44.61090449283989, 44.83363942288998]
	indexes_1 = [1314, 2628, 3942, 5256, 6570, 7884, 9198, 10512, 11826, 13148]
	speeds_ms_1 = convert_kmh_to_ms(speeds_1)
	minimized_speed_vector_1 = set_speeds(speeds_ms_1, short_track, indexes_1)
	power_use_1, solar_power_1, energy_in_system_1, time_1, time_s_1 = solar_trip_calculation_bounds(minimized_speed_vector_1, short_track, start_datetime_hierarchical, start_energy_short)
	
	
	plot(short_track.distance, short_track.altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance for 1st iteration", right_margin = 15Plots.mm)
	plot!(twinx(), short_track.distance, minimized_speed_vector_1 * 3.6, color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance for 1st iteration")
end

# ╔═╡ d8e62070-5595-4e5f-b066-77f9e7d816b1
plot(short_track.distance, [power_use_1 solar_power_1 energy_in_system_1 zeros(track_size)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Hierarchical iter 1",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, #size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)

# ╔═╡ e7b301c0-7c5d-4d45-a70f-91653e81ea0d
time_1

# ╔═╡ 141378f3-9b42-4b83-a87c-44b3ba3928aa
function hierarchical_optimization(speed, track, chunks_amount, start_energy, finish_energy, start_datetime, iteration, end_index)
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
		return solar_partial_trip_wrapper(speed, track, split_indexes, start_energy, finish_energy, start_datetime)
	end
	td = TwiceDifferentiable(f, fill(speed, chunks_amount); autodiff = :forward)
	lower_bound = fill(0.0, chunks_amount)
	upper_bound = fill(100.0, chunks_amount)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
	line_search = LineSearches.BackTracking();
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
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds(minimized_speed_vector, track, start_datetime, start_energy)
	# println("iteration $(iteration), speed is $(speed) planned finish energy is $(finish_energy)")
	# println("track from $(start_index) to $(end_index)")
	# println("start datetime $(start_datetime)")
	# println("stare energy $(start_energy)")
	# println("time is $(last(time_s)), cost is $(f(minimized_speeds))")
	# println("split indexes are $(split_indexes)")
	# println("distances are $(track[split_indexes, :])")
	# println("minimized speeds are: $(minimized_speeds)")
	# println("simulated finish energy is $(last(energy_in_system))")
	# println("calculated cost is $( last(time_s) + 100 * abs(last(energy_in_system) - finish_energy) + 100 * sum(abs.(energy_in_system[energy_in_system .< 0.0])) + 100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .< 0.0])) + 100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .> 100.0])) )")
	# println("finish energy difference penalty is: $(100 * abs(last(energy_in_system) - finish_energy))")
	# println("energy less than 0. penalty is: $(100 * sum(abs.(energy_in_system[energy_in_system .< 0.0])))")
	# println("speed less than 0. penalty is: $(100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .< 0.0])))")
	# println("speed more than 100. penalty is: $(100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .> 100.0 / 3.6])))")
	split_energies = energy_in_system[split_indexes]
	pushfirst!(split_energies, start_energy)
	split_times = time[split_indexes, :utc_time]
	pushfirst!(split_times, start_datetime)
	# println("split energies are $(split_energies)")
	# println("")
	
	# 5 - go though each track piece and enter function again
	# hierarchical_optimization(minimized_speeds[1], tracks[1], chunks_amount, start_energy, split_energies[1], start_datetime, iteration + 1)
	# @debug "split_energies size is $(size(split_energies, 1)), chunks_amount is $(chunks_amount)"
	result_speeds = []
	for i=1:chunks_amount
		result_speeds_chunk = hierarchical_optimization(minimized_speeds[i], tracks[i], chunks_amount, split_energies[i], split_energies[i+1], split_times[i], iteration + 1 , start_index + split_indexes[i] - 1)
		append!(result_speeds, result_speeds_chunk)
	end

	return result_speeds
end

# ╔═╡ d9e30de2-e75f-423b-8fcc-ab3847331274
@time result_hierarchical = hierarchical_optimization(initial_speed, short_track, chunks_amount_hierarchical, start_energy_short, 0., start_datetime_hierarchical, 1, track_size)

# ╔═╡ b55c819b-f312-4078-b751-cf443355be19
begin
	inputs_ms_hier = abs.(convert_kmh_to_ms(result_hierarchical))
	power_use_hier, solar_power_hier, energy_in_system_hier, time_hier, time_s_hier = solar_trip_calculation_bounds(inputs_ms_hier, short_track, start_datetime_hierarchical, start_energy_short)
	last(time_s_hier)
end

# ╔═╡ db7cc403-40ba-47d3-b27b-2a5913633ae5
plot(short_track.distance, [power_use_hier solar_power_hier energy_in_system_hier zeros(track_size)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Hierarchical",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, #size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)

# ╔═╡ 3642f56d-ac97-435a-b446-68eb7814b03c
begin
	plot(short_track.distance, short_track.altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance", right_margin = 15Plots.mm)
	plot!(twinx(), short_track.distance, result_hierarchical, color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance")
end

# ╔═╡ 6732e4c8-cd5f-454c-84fb-14aae6c02fbe
begin
	plot(short_track[1000:1050,:].distance, short_track[1000:1050,:].altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance [1000:1050]", right_margin = 15Plots.mm)
	plot!(twinx(), short_track[1000:1050,:].distance, result_hierarchical[1000:1050], color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance [1000:1050]")
end

# ╔═╡ aab805b0-2fe4-4ece-8f8c-f922226a9912
@time result_hierarchical_bad = hierarchical_optimization(26.90116336126291, short_track[25:28,:], chunks_amount_hierarchical, 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8), 1, 28)

# ╔═╡ 99a41fe8-cc03-49fe-a39f-c27070fda62e
function solar_trip_cost_energy(input_speed, track, start_energy)
	# @debug "func solar_trip_cost_energy, speeds size is $(size(input_speed, 1)), track size is $(size(track.distance, 1))"
    power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation(input_speed, track, start_energy)
    cost = last(time_s) + 10 * abs(minimum(energy_in_system)) + 100 * sum(abs.(input_speed[input_speed .< 0.0])) + 10 * abs(last(energy_in_system))
    return cost
end

# ╔═╡ 96ee2227-0386-4a86-93ad-fa41dcdbf615
function solar_trip_chunks_energy(speeds, track, start_energy)
# function solar_trip_chunks(speeds::Vector{Float64}, track)
    speed_ms = convert_kmh_to_ms(speeds)
    # speed_vector = propagate_speeds(speed_ms, track)
	# @debug "func solar_trip_chunks_energy, speeds size is $(size(speed_ms, 1)), track size is $(size(track.distance, 1))"
    return solar_trip_cost_energy(speed_ms, track, start_energy)
end

# ╔═╡ 38142a0b-1e1f-4cda-90e9-50d52e6f0227
function f_short(x)
	return solar_trip_chunks_energy(abs.(x), short_track, start_energy_short)
end

# ╔═╡ 9d53a659-612c-405f-88da-253ffe57f4a0
f_short(fill(initial_speed, track_size))

# ╔═╡ 87b00642-3d8a-4722-bb7b-b9ff12a26d8f
td_short = TwiceDifferentiable(f_short, fill(initial_speed, track_size), autodiff = :forward)

# ╔═╡ 785c4601-6d60-49b6-9470-abe90eed9d4c
begin
	# 1 - splitting the track
	# determine split indexes
	track_size_2_bad = 2628-1315+1
	start_index_2 = 1315
	end_index_2 = 2628
	chunks_amount_2 = 10
	split_indexes_2 = calculate_split_indexes(track_size_2_bad, chunks_amount_2)
	# for the case when there are less indexes than chunks
	chunks_amount_2 = size(split_indexes_2,1)

	start_energy_2 = 4990.610198979855
	start_datetime_2 = DateTime(2022,7,2,8,41,20)
	finish_energy_2 = 4990.257345152248

	# 2 - set up optimization itself
	function f_2(speed)
		return solar_partial_trip_wrapper(speed, short_track[start_index_2:end_index_2,:], split_indexes_2, start_energy_2, finish_energy_2, start_datetime_2)
	end
	td_2 = TwiceDifferentiable(f_2, fill(60., chunks_amount_2); autodiff = :forward)
	lower_bound_2 = fill(30.0, chunks_amount_2)
	upper_bound_2 = fill(100.0, chunks_amount_2)
	tdc_2 = TwiceDifferentiableConstraints(lower_bound_2, upper_bound_2)
	# result = optimize(td, fill(speed, chunks_amount),
	    #Newton(; linesearch = line_search),
	result_2 = optimize(td_2, tdc_2, fill(60., chunks_amount_2),
		IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10,
			allow_f_increases = true,
			successive_f_tol = 50,
			show_trace = false,
			extended_trace = false
	    )
	)

	# 3 - save optimized speeds
	minimized_speeds_2 = Optim.minimizer(result_2)
	println(minimized_speeds_2)
	
	# 4 - sumulate again to obtain energies and times around split indexes
	minimized_speeds_ms_2 = convert_kmh_to_ms(minimized_speeds_2/3.6)
	minimized_speed_vector_2 = set_speeds(minimized_speeds_ms_2, short_track[start_index_2:end_index_2,:], split_indexes_2)
	power_use_2, solar_power_2, energy_in_system_2, time_2, time_s_2 = solar_trip_calculation_bounds(minimized_speed_vector_2, short_track[start_index_2:end_index_2,:], start_datetime_2, start_energy_2)
	println(last(energy_in_system_2))

	plot(short_track[start_index_2:end_index_2,:].distance, short_track[start_index_2:end_index_2,:].altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance for 2nd iteration", right_margin = 15Plots.mm)
	plot!(twinx(), short_track[start_index_2:end_index_2,:].distance, minimized_speed_vector_2 * 3.6, color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance for 2nd iteration")
end

# ╔═╡ 8635aede-c512-470a-b1d0-f53d783c6179
begin
	plot(short_track[1315:2628,:].distance, [power_use_1[1315:2628] solar_power_1[1315:2628] energy_in_system_1[1315:2628] zeros(track_size_2_bad)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Hierarchical iter 1 region",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, #size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)
end

# ╔═╡ 802e662b-0922-4993-bcf6-059c34ce02b7
split_indexes_2

# ╔═╡ d5a4ad7a-2d5c-4597-947e-3f2e6b9e67e9
begin
	plot(short_track[start_index_2:end_index_2,:].distance, [power_use_2 solar_power_2 energy_in_system_2 zeros(track_size_2_bad)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Hierarchical iter 2",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, #size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)
end

# ╔═╡ ce217af4-f88e-4767-984b-557eaf007bfb
time_2

# ╔═╡ 0a54b576-9d04-42ee-a19f-5ecafed7b2d6
solar_income_2 = solar_power_income(time_2, short_track[start_index_2:end_index_2,:], minimized_speed_vector_2)

# ╔═╡ ef772852-003d-4c0c-b360-62941a8c8357
plot(short_track[start_index_2:end_index_2,:].distance, solar_income_2)

# ╔═╡ c9c142d9-2bfb-4722-8d7e-a8c6724f3351
plot(short_track[start_index_2:end_index_2,:].distance./minimized_speed_vector_2, solar_income_2)

# ╔═╡ d75bf2cc-5a06-4cd6-b4ce-ed0ce415fd79
plot(short_track[start_index_2:end_index_2,:].distance, solar_income_2_no_speed)

# ╔═╡ c7197a1e-d135-46e5-88b1-771934d58bf5
plot(time_2.time_s, solar_income_2_no_speed)

# ╔═╡ 63245989-6885-4d7c-bedd-cece121bbdef
plot(time_2.utc_time, solar_income_2_no_speed)

# ╔═╡ 9ad839d3-40cd-438b-9a93-aaa0635ec0a1
time_2

# ╔═╡ 2c6fbc64-c7f3-430d-b6a8-c4a78f49c9c4
calculate_travel_time_datetime(minimized_speed_vector_2, short_track[start_index_2:end_index_2,:], DateTime(2022,7,2,8,41,20))

# ╔═╡ 1c14c9e2-a05c-4301-bfe6-69f36d0be865
begin
	minimized_speeds_2_same = [43.228254912151826, 43.228254912151826, 43.228254912151826, 43.228254912151826, 43.228254912151826, 43.228254912151826, 43.228254912151826, 43.228254912151826, 43.228254912151826, 43.228254912151826]
	
	# 4 - sumulate again to obtain energies and times around split indexes
	minimized_speeds_ms_2_same = convert_kmh_to_ms(minimized_speeds_2_same)
	minimized_speed_vector_2_same = set_speeds(minimized_speeds_ms_2_same, short_track[start_index_2:end_index_2,:], split_indexes_2)
	power_use_2_same, solar_power_2_same, energy_in_system_2_same, time_2_same, time_s_2_same = solar_trip_calculation_bounds(minimized_speed_vector_2_same, short_track[start_index_2:end_index_2,:], start_datetime_2, start_energy_2)
	println(last(energy_in_system_2_same))

	plot(short_track[start_index_2:end_index_2,:].distance, short_track[start_index_2:end_index_2,:].altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance for 2nd iteration same speed", right_margin = 15Plots.mm)
	plot!(twinx(), short_track[start_index_2:end_index_2,:].distance, minimized_speed_vector_2_same * 3.6, color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance for 2nd iteration same speed")
end

# ╔═╡ c8d5a429-ea0e-4730-bc15-4b1e49642407
begin
	plot(short_track[start_index_2:end_index_2,:].distance, [power_use_2_same solar_power_2_same energy_in_system_2_same zeros(track_size_2_bad)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Hierarchical iter 2 same speed",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, #size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)
end

# ╔═╡ d08e3a6e-9ca6-4ab9-8228-5a24546e1d36
time_test_time = start_test .+ Dates.Millisecond.(round.(time_s_test .* 1000))

# ╔═╡ 2143a79c-efe9-47df-88c8-eda9c2e16623
plot(time_test_time)

# ╔═╡ 5629455d-a6ce-430d-862d-5467ffc50ac9
plot(Dates.hour.(time_test_time))

# ╔═╡ bb508865-1543-4621-abfd-f79a90f80db6
input_speed_solar_test = fill(10, size(time_test_time, 1)) #m/s

# ╔═╡ 93aafc03-dd14-4ca8-9226-5bbf8d42f378
size(input_speed_solar_test, 1)

# ╔═╡ a13b519f-10e2-49eb-9903-8f198676cb76
time_test_df = DataFrame(time_s=input_speed_solar_test)

# ╔═╡ 840dd334-16ad-4645-9773-393497547ae6
spi_10s = solar_power_income(time_test_df, track_test_df)

# ╔═╡ 22d50e51-2cee-474b-87f9-28f73cfd25b1
plot(time_test_df.utc_time, spi_10s)

# ╔═╡ 71291ca2-6ab5-43a8-9f54-58990e4dcbb3
spi_10s_speed = solar_power_income(time_test_df, track_test_df, speed_spi_test)

# ╔═╡ 233a5ed0-00eb-4b4e-9225-5838ceddb187
plot(time_test_df.utc_time, spi_10s_speed)

# ╔═╡ a92a0f58-0462-4a07-9053-a89caf5be9bb
time_test_df.utc_time = DateTime(2022,7,1,0,0,0) .+ Dates.Millisecond.(round.(time_s_test .* 1000))

# ╔═╡ 70e5a89c-4694-4a5a-906e-73b2adbd1336
function solar_partial_trip_test_wrapper(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds(speeds_ms, track, indexes)

	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) + 10 * abs(last(energy_in_system) - finish_energy) + sum(energy_in_system[energy_in_system .< 0.0])
	
	# + 100 * sum(abs.(speed_vector[speed_vector .< 0.0])) + 100 * sum(abs.(speed_vector[speed_vector .> 100.0 / 3.6]))
	# + 100 * abs(minimum(energy_in_system) - finish_energy) 
	
	# + 100 * sum(speed_vector[speed_vector .> 200.0]) + 100 * sum(energy_in_system[energy_in_system .< 0.0])


	# cost = last(time_s) + 10 * abs(minimum(energy_in_system)) + 100 * sum(input_speed[input_speed .< 0.0])
	return cost
end

# ╔═╡ 3a4a992b-3050-43c0-b1f2-c7965a48d202
function cost_calc(speed1, speed2)
	return solar_partial_trip_test_wrapper([speed1, speed2], track_2_pcs, [1, 2], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
end

# ╔═╡ 8ef52f08-6518-4c82-912e-b3a20f4c81a7
Plots.surface(30:50, 30:50, cost_calc)

# ╔═╡ c289432f-5827-4651-8e2c-f6bd31d8fe24
Plots.surface(5:5:80, 5:5:80, cost_calc)

# ╔═╡ 45c6e6ca-d9d5-4895-829f-3522bb4485e0
begin

	# function f_test(speed)
	# 	return solar_partial_trip_wrapper(speed, short_track[25:28,:], [1,2,3,4], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
	# end
	# td_test = TwiceDifferentiable(f_test, fill(26.90116336126291, 4); autodiff = :forward)
	# lower_bound_test = fill(0.0, 4)
	# upper_bound_test = fill(100.0, 4)
	# tdc_test = TwiceDifferentiableConstraints(lower_bound_test, upper_bound_test)

	
	function f_test_plane(speed_inp)
		return solar_partial_trip_test_wrapper(speed_inp, track_2_pcs, [1, 2], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
	end
	td_test_2d = TwiceDifferentiable(f_test_plane, [25.0, 25.0]; autodiff = :forward)
	lower_bound_test_2d = fill(0.0, 2)
	upper_bound_test_2d = fill(100.0, 2)
	tdc_test_2d = TwiceDifferentiableConstraints(lower_bound_test_2d, upper_bound_test_2d)
	# line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, chunks_amount),
	    #Newton(; linesearch = line_search),
	# inputs_test = fill(26.90116336126291, 4)
	inputs_test_2d = [25.0, 25.0]
	@time result_test_2d = optimize(td_test_2d, tdc_test_2d, inputs_test_2d,
		IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10,
			allow_f_increases = true,
			successive_f_tol = 50,
			show_trace = false
	    )
	)
end

# ╔═╡ 474bfcdf-7c28-472a-acb0-07aa7fcc719d
Optim.minimizer(result_test_2d)

# ╔═╡ efe841ab-9d93-4792-917e-6c27916381af
begin

	# function f_test(speed)
	# 	return solar_partial_trip_wrapper(speed, short_track[25:28,:], [1,2,3,4], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
	# end
	# td_test = TwiceDifferentiable(f_test, fill(26.90116336126291, 4); autodiff = :forward)
	# lower_bound_test = fill(0.0, 4)
	# upper_bound_test = fill(100.0, 4)
	# tdc_test = TwiceDifferentiableConstraints(lower_bound_test, upper_bound_test)

	
	function f_test_plane_3d(speed_inp)
		return solar_partial_trip_test_wrapper(speed_inp, track[25:27,:], [1, 2, 3], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
	end
	td_test_3d = TwiceDifferentiable(f_test_plane_3d, [25.0, 25.0, 25.]; autodiff = :forward)
	lower_bound_test_3d = fill(0.0, 3)
	upper_bound_test_3d = fill(100.0, 3)
	tdc_test_3d = TwiceDifferentiableConstraints(lower_bound_test_3d, upper_bound_test_3d)
	# line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, chunks_amount),
	    #Newton(; linesearch = line_search),
	# inputs_test = fill(26.90116336126291, 4)
	inputs_test_3d = [25.0, 25.0, 25.]
	@time result_test_3d = optimize(td_test_3d, tdc_test_3d, inputs_test_3d,
		IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10,
			allow_f_increases = true,
			successive_f_tol = 50,
			show_trace = false
	    )
	)
end

# ╔═╡ 9c58b595-b903-4bfc-bf4d-6e7cf8acf588
Optim.minimizer(result_test_3d)

# ╔═╡ c5d6352a-fd03-4045-a91a-8574c15d9989
begin
	dim = 50
	start = 20
	cost_results = zeros(dim,dim)
	for i=1:dim
		for j=1:dim
			cost_results[i, j] = solar_partial_trip_test_wrapper([start+i*1, start+j*1], track_2_pcs, [1, 2], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
		end
	end
end

# ╔═╡ 92281858-95ad-4748-8088-cc6c74ae16dc
cost_results

# ╔═╡ 5797e1bf-7bf8-4343-b740-6f09c433f6c5
minimum(cost_results)

# ╔═╡ e3ea8615-1bd2-40f3-bb79-c27e3c5375fe
argmin(cost_results)

# ╔═╡ c4bfa3e8-142b-47a9-9aa4-0457cd15f6c4
cost_calc(20+19, 20+19)

# ╔═╡ ccd02e1b-b33c-4aec-a95f-da12dbb4cda7
function solar_partial_trip_grad_wrapper(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds_grad(speeds_ms, track, indexes)

	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) + 10 * abs(last(energy_in_system) - finish_energy) + sum(energy_in_system[energy_in_system .< 0.0])
	
	# + 100 * sum(abs.(speed_vector[speed_vector .< 0.0])) + 100 * sum(abs.(speed_vector[speed_vector .> 100.0 / 3.6]))
	# + 100 * abs(minimum(energy_in_system) - finish_energy) 
	
	# + 100 * sum(speed_vector[speed_vector .> 200.0]) + 100 * sum(energy_in_system[energy_in_system .< 0.0])


	# cost = last(time_s) + 10 * abs(minimum(energy_in_system)) + 100 * sum(input_speed[input_speed .< 0.0])
	return cost
end

# ╔═╡ 68b614f2-dc13-4c46-9649-273179fbe27c
function f_test_plane_3d_grad(speed_inp)
		return solar_partial_trip_grad_wrapper(speed_inp, track[25:27,:], [1, 2, 3], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
	end

# ╔═╡ 2a19189a-8da9-4e1a-a76b-6d4a09e0ad0f
function solar_radiation_no_broadcasting(data_df, track_dataframe)
    # starts here: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-outside-the-earths-atmosphere

    # copy generate_year_time_dataframe(100000)
    # data_df = copy(time_dataframe)
    data_df.latitude = track_dataframe.latitude
    data_df.longitude = track_dataframe.longitude
    data_df.altitude = track_dataframe.altitude

    # for debug purposes
    # data_df = copy(time_df)
    # data_df.latitude = track.latitude
    # data_df.longitude = track.longitude
    # data_df.altitude = track.altitude

    # setting up needed values for calculcations
    # data_df.day = Dates.day.(data_df.datetime)

    # solar constant
    H_constant = 1353 # W/m^2
    # AM0 = 1366 # air mass zero, W/m^2

    # radiant power density otside the Erath's athmosphere in W/m^2
    # H = H_constant * (1 .+ 0.033*cosd.( (360 * (Dates.dayofyear.(data_df.utc_time) .- 2)) / 365 ) )
    # plot(Dates.dayofyear.(data_df.utc_time),H, title = "Radial power density W/m^2")

    # TODO: air mass calculation with phi? and theta?

    # time shift for current latitude, hours
    # local_standard_time_meridian = 15 * Dates.Hour(local_time - UTC_time).value
    # very rough estimation
    # proper solution like https://stackoverflow.com/questions/5584602/determine-timezone-from-latitude-longitude-without-using-web-services-like-geona
    # data_df.lstm = 15 * ceil(data_df.longitude * 24 / 360)

    # equation of time
    B = 360 / 365 * (Dates.dayofyear(data_df.utc_time) - 81) # in degrees
    equation_of_time = 9.87 * sind(2B) - 7.53 * cosd(B) - 1.5 * sind(B) # minutes
    # plot(
    #     equation_of_time,
    #     title = "Equation of time",
    #     label = "Equation of time",
    #     xlabel = "Day of the Year",
    #     ylabel = "Minutes"
    # )

    # TODO: resume from here
    # because local standard time meridian is not real sun time , minutes
    # time_correction_factor = 4 * (data_df.longitude - data_df.lstm) + equation_of_time
    # plot(time_correction_factor, title="Time correction factor")

    # local_solar_time = local_time + Dates.Second(round(time_correction_factor * 60.0))
    # # LST = 4 longitude + 60 UTC_minutes +EoT
    # # so, we only need to know UTC time and longitude
    # local_solar_time = UTC_time +
    #     Dates.Second( round( (4*longitude + equation_of_time[day]) * 60) )

    local_solar_time = data_df.utc_time +
        Dates.Second( round(4 * data_df.longitude + equation_of_time) * 60)
    
    # TODO: should hour angle be only in integer values? 
    # minutes_from_start_of_the_day = Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time);
    hour_angle = 15 * ((Dates.hour(data_df.utc_time) * 60 + Dates.minute(data_df.utc_time)) /60 - 12)
    # TODO: hour angle calculation for UTC time
    # plot(data_df.utc_time, hour_angle, title="Hour angle vs UTC time")
    # TODO: hour angle calculation for local solar time ??? check if needed

    # TODO: continue to sun declination angle
    # TODO: continous sun declination angle calculation? (as for now quantified by days)
    sun_declination_angle = -23.45 * cosd(360 / 365 * (Dates.dayofyear(data_df.utc_time) + 10))
    # plot(sun_declination_angle, title = "Sun declination angle")

    # elevation (altitude) angle
    # elevation_angle = 90 .+ data_df.latitude .- sun_declination_angle
    # plot(data_df.utc_time, elevation_angle, title="Elevation angle (deg)")

    # elevation?
    elevation = asin(
        sind(sun_declination_angle) * sind(data_df.latitude) +
        cosd(sun_declination_angle) * cosd(data_df.latitude) * cosd.(hour_angle)
    )
    # plot(elevation, title="Elevation of the sun in rad")
    # when elevation angle is negative it means that sun is below horizon 
    # so, we should nullify these elements
    elevation_filtered = copy(elevation)
	if elevation_filtered <= 0
		elevation_filtered = 0
	end
    # elevation_filtered[elevation_filtered .<=0] .= 0
	
    # plot(data_df.utc_time, elevation_filtered, title="elevation angle when sun above horizon")

    # TODO: sunrise and sunset hours

    # sunrise_hour = 12 .- 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60
    # # plot(data_df.utc_time, sunrise_hour, title = "sunrise hour")

    # sunset_hour = 12 .+ 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60

    # plot(data_df.utc_time, sunset_hour, title = "sunset hour")

    # calculating zenith angle, needed for optical air mass
    zenith_angle = pi/2 - elevation_filtered
    # plot(data_df.utc_time, zenith_angle, title="Zenith angle of the sun in rad")
    cos_zenith_angle = cos(zenith_angle)

    # solar radiation consists of direct and diffuse radiation
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass 

    # calculate air mass - how much air light must pass through
    air_mass_full = 1 / (cos_zenith_angle + 0.50572 * (96.07995 - deg2rad(zenith_angle))^-1.6364 )

    # direct intensity, not actually used
    # intensity_direct = H_constant .* 0.7 .^ (air_mass_full .^ 0.678) # W/m^2
    # H_constant * 70% of radiation ^ air_mass_full ^ coefficient to fit the experimental data
    # plot(data_df.utc_time, intensity_direct, title="Direct intensity at Earth ground, W/m^2")

    # with altitude (from sea level)
    intensity_direct_altitude = H_constant * ( (1 - 0.14 * data_df.altitude / 1000) * 0.7 ^ (air_mass_full ^ 0.678) + 0.14 * data_df.altitude / 1000) # W/m^2
    # plot(data_df.utc_time, intensity_direct_altitude, title="Direct intensity at altitude, W/m^2")
    # diffuse radiation
    intensity_diffuse = intensity_direct_altitude * 0.1
    # global irradiance
    intensity_global = intensity_direct_altitude + intensity_diffuse
    # plot(data_df.utc_time, intensity_global, title = "Global irradiance at altitude")
    # irradiance * cos(zenith_angle) is incident radiation ?

    
    s_incident = intensity_global
    # TODO: calculate radiation on a tilted surface
    # can be done through perpendecular (incident) or horizontal
    # s_horizontal = s_incident .* sin.(elevation)
    # module_angle_rad - at what angle to surface panels are located
    module_angle_rad = 0.0 # just on the ground
    # azimuth_angle = 0.0 # just something, since it is on the ground, will not be used
    # TODO: calculate sun azimuth angle according to https://www.pveducation.org/pvcdrom/properties-of-sunlight/azimuth-angle 
    azimuth_cos = ( sind(sun_declination_angle) * cos(data_df.latitude) - 
    cosd(sun_declination_angle) * sin(data_df.latitude) * cos(hour_angle) ) / cos(elevation)
	if azimuth_cos > 1
		azimuth_cos = 1
	elseif azimuth_cos < -1
		azimuth_cos = -1
	end
    # azimuth_cos[azimuth_cos .> 1] .= 1
    # azimuth_cos[azimuth_cos .< -1] .= -1
    azimuth_angle = acos( azimuth_cos )
    lst = Dates.hour(local_solar_time)
	if Dates.hour(local_solar_time) > 12
		azimuth_angle = 2*pi - azimuth_angle
	end
    # data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ] .= 2*pi .- data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ]
		
    # module azimuth angle
    module_azimuth_angle = azimuth_angle
    # s_module stand for solar module, tilted surface
    # s_module = s_incident .* sin.(elevation .+ data_df.module_angle_rad)
    # plot(data_df.utc_time, s_module, title="Solar intensity without module azimuth")
    # OR https://www.pveducation.org/pvcdrom/properties-of-sunlight/arbitrary-orientation-and-tilt 
    s_module_azimuth = s_incident * (
        cos(elevation) * sin(module_angle_rad) * cos(azimuth_angle - module_azimuth_angle) +
        sin(elevation) * cos(module_angle_rad)
        )
    # plot(data_df.utc_time, s_module_azimuth, title="Solar intensity with module azimuth")

    # next - simulate the race
    return s_module_azimuth
end

# ╔═╡ 9020198b-e68e-4253-a6ed-31a8d76944a6
function solar_power_income_no_broadcasting(time_df, track_df, speed_vector)
    electrics_efficiency = 0.86
    solar_panels_efficiency = 0.228
    panels_area = 4 # m^2
    solar_intensity = solar_radiation_no_broadcasting(time_df, track_df)
    power_income = electrics_efficiency * solar_panels_efficiency * solar_intensity * panels_area * track_df.diff_distance / speed_vector # W*s
    power_income_wt_h = power_income / 3600
    # power_income(i) = electrics_efficiency * panels_efficiency * solar_rad_h * panels_area * ...
    # ( dist_diff(i) * 3.6 / ( speed * 3600 ) ); % kWt*h
    return power_income_wt_h 
end

# ╔═╡ c22bfed6-38cf-454f-93b1-6516e10e3b96
solar_power_income_no_broadcasting(time_df_no_b[1,:], short_track[1,:], speed_vec_no_b[1])

# ╔═╡ a37676f5-e278-4df8-aa88-4f3ec1e3c414
solar_power_income_no_broadcasting(first(time_df_no_b), first(short_track), first(speed_vec_no_b))

# ╔═╡ 45a74d26-5ad9-404d-af14-33ae4f863a3d
solar_power_income_no_broadcasting(time_df_no_b, short_track, speed_vec_no_b)

# ╔═╡ 7f7ab4bd-69d4-498e-b6f1-a01bb5857ce5
broadcast(solar_power_income_no_broadcasting, time_df_no_b, short_track, speed_vec_no_b)

# ╔═╡ 15186cf0-a932-47b6-9f4c-3378b123de04
function solar_radiation_no_broadcasting_arrays(latitude, longitude, altitude, utc_time)
    # starts here: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-outside-the-earths-atmosphere

    # copy generate_year_time_dataframe(100000)

    # for debug purposes
    # data_df = copy(time_df)
    # data_df.latitude = track.latitude
    # data_df.longitude = track.longitude
    # data_df.altitude = track.altitude

    # setting up needed values for calculcations
    # data_df.day = Dates.day.(data_df.datetime)

    # solar constant
    H_constant = 1353 # W/m^2
    # AM0 = 1366 # air mass zero, W/m^2

    # radiant power density otside the Erath's athmosphere in W/m^2
    # H = H_constant * (1 .+ 0.033*cosd.( (360 * (Dates.dayofyear.(data_df.utc_time) .- 2)) / 365 ) )
    # plot(Dates.dayofyear.(data_df.utc_time),H, title = "Radial power density W/m^2")

    # TODO: air mass calculation with phi? and theta?

    # time shift for current latitude, hours
    # local_standard_time_meridian = 15 * Dates.Hour(local_time - UTC_time).value
    # very rough estimation
    # proper solution like https://stackoverflow.com/questions/5584602/determine-timezone-from-latitude-longitude-without-using-web-services-like-geona
    # data_df.lstm = 15 * ceil(data_df.longitude * 24 / 360)

    # equation of time
    B = 360 / 365 * (Dates.dayofyear(utc_time) - 81) # in degrees
    equation_of_time = 9.87 * sind(2B) - 7.53 * cosd(B) - 1.5 * sind(B) # minutes
    # plot(
    #     equation_of_time,
    #     title = "Equation of time",
    #     label = "Equation of time",
    #     xlabel = "Day of the Year",
    #     ylabel = "Minutes"
    # )

    # TODO: resume from here
    # because local standard time meridian is not real sun time , minutes
    # time_correction_factor = 4 * (data_df.longitude - data_df.lstm) + equation_of_time
    # plot(time_correction_factor, title="Time correction factor")

    # local_solar_time = local_time + Dates.Second(round(time_correction_factor * 60.0))
    # # LST = 4 longitude + 60 UTC_minutes +EoT
    # # so, we only need to know UTC time and longitude
    # local_solar_time = UTC_time +
    #     Dates.Second( round( (4*longitude + equation_of_time[day]) * 60) )

    local_solar_time = utc_time +
        Dates.Second( round(4 * longitude + equation_of_time) * 60)
    
    # TODO: should hour angle be only in integer values? 
    # minutes_from_start_of_the_day = Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time);
    hour_angle = 15 * ((Dates.hour(utc_time) * 60 + Dates.minute(utc_time)) /60 - 12)
    # TODO: hour angle calculation for UTC time
    # plot(data_df.utc_time, hour_angle, title="Hour angle vs UTC time")
    # TODO: hour angle calculation for local solar time ??? check if needed

    # TODO: continue to sun declination angle
    # TODO: continous sun declination angle calculation? (as for now quantified by days)
    sun_declination_angle = -23.45 * cosd(360 / 365 * (Dates.dayofyear(utc_time) + 10))
    # plot(sun_declination_angle, title = "Sun declination angle")

    # elevation (altitude) angle
    # elevation_angle = 90 .+ data_df.latitude .- sun_declination_angle
    # plot(data_df.utc_time, elevation_angle, title="Elevation angle (deg)")

    # elevation?
    elevation = asin(
        sind(sun_declination_angle) * sind(latitude) +
        cosd(sun_declination_angle) * cosd(latitude) * cosd.(hour_angle)
    )
    # plot(elevation, title="Elevation of the sun in rad")
    # when elevation angle is negative it means that sun is below horizon 
    # so, we should nullify these elements
    elevation_filtered = copy(elevation)
	if elevation_filtered <= 0
		elevation_filtered = 0
	end
    # elevation_filtered[elevation_filtered .<=0] .= 0
	
    # plot(data_df.utc_time, elevation_filtered, title="elevation angle when sun above horizon")

    # TODO: sunrise and sunset hours

    # sunrise_hour = 12 .- 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60
    # # plot(data_df.utc_time, sunrise_hour, title = "sunrise hour")

    # sunset_hour = 12 .+ 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60

    # plot(data_df.utc_time, sunset_hour, title = "sunset hour")

    # calculating zenith angle, needed for optical air mass
    zenith_angle = pi/2 - elevation_filtered
    # plot(data_df.utc_time, zenith_angle, title="Zenith angle of the sun in rad")
    cos_zenith_angle = cos(zenith_angle)

    # solar radiation consists of direct and diffuse radiation
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass 

    # calculate air mass - how much air light must pass through
    air_mass_full = 1 / (cos_zenith_angle + 0.50572 * (96.07995 - deg2rad(zenith_angle))^-1.6364 )

    # direct intensity, not actually used
    # intensity_direct = H_constant .* 0.7 .^ (air_mass_full .^ 0.678) # W/m^2
    # H_constant * 70% of radiation ^ air_mass_full ^ coefficient to fit the experimental data
    # plot(data_df.utc_time, intensity_direct, title="Direct intensity at Earth ground, W/m^2")

    # with altitude (from sea level)
    intensity_direct_altitude = H_constant * ( (1 - 0.14 * altitude / 1000) * 0.7 ^ (air_mass_full ^ 0.678) + 0.14 * altitude / 1000) # W/m^2
    # plot(data_df.utc_time, intensity_direct_altitude, title="Direct intensity at altitude, W/m^2")
    # diffuse radiation
    intensity_diffuse = intensity_direct_altitude * 0.1
    # global irradiance
    intensity_global = intensity_direct_altitude + intensity_diffuse
    # plot(data_df.utc_time, intensity_global, title = "Global irradiance at altitude")
    # irradiance * cos(zenith_angle) is incident radiation ?

    
    s_incident = intensity_global
    # TODO: calculate radiation on a tilted surface
    # can be done through perpendecular (incident) or horizontal
    # s_horizontal = s_incident .* sin.(elevation)
    # module_angle_rad - at what angle to surface panels are located
    module_angle_rad = 0.0 # just on the ground
    # azimuth_angle = 0.0 # just something, since it is on the ground, will not be used
    # TODO: calculate sun azimuth angle according to https://www.pveducation.org/pvcdrom/properties-of-sunlight/azimuth-angle 
    azimuth_cos = ( sind(sun_declination_angle) * cos(latitude) - 
    cosd(sun_declination_angle) * sin(latitude) * cos(hour_angle) ) / cos(elevation)
	if azimuth_cos > 1
		azimuth_cos = 1
	elseif azimuth_cos < -1
		azimuth_cos = -1
	end
    # azimuth_cos[azimuth_cos .> 1] .= 1
    # azimuth_cos[azimuth_cos .< -1] .= -1
    azimuth_angle = acos( azimuth_cos )
    lst = Dates.hour(local_solar_time)
	if Dates.hour(local_solar_time) > 12
		azimuth_angle = 2*pi - azimuth_angle
	end
    # data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ] .= 2*pi .- data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ]
		
    # module azimuth angle
    module_azimuth_angle = azimuth_angle
    # s_module stand for solar module, tilted surface
    # s_module = s_incident .* sin.(elevation .+ data_df.module_angle_rad)
    # plot(data_df.utc_time, s_module, title="Solar intensity without module azimuth")
    # OR https://www.pveducation.org/pvcdrom/properties-of-sunlight/arbitrary-orientation-and-tilt 
    s_module_azimuth = s_incident * (
        cos(elevation) * sin(module_angle_rad) * cos(azimuth_angle - module_azimuth_angle) +
        sin(elevation) * cos(module_angle_rad)
        )
    # plot(data_df.utc_time, s_module_azimuth, title="Solar intensity with module azimuth")

    # next - simulate the race
    return s_module_azimuth
end

# ╔═╡ 01382010-99aa-4393-b6e0-7eaefd2b2ace
function solar_power_income_expanded(latitude, longitude, altitude, utc_time, diff_distance, speed_vector)
    electrics_efficiency = 0.86
    solar_panels_efficiency = 0.228
    panels_area = 4 # m^2
    solar_intensity = solar_radiation_no_broadcasting_arrays(latitude, longitude, altitude, utc_time)
    power_income = electrics_efficiency * solar_panels_efficiency * solar_intensity * panels_area * diff_distance / speed_vector # W*s
    power_income_wt_h = power_income / 3600
    # power_income(i) = electrics_efficiency * panels_efficiency * solar_rad_h * panels_area * ...
    # ( dist_diff(i) * 3.6 / ( speed * 3600 ) ); % kWt*h
    return power_income_wt_h 
end

# ╔═╡ 89b6d6ee-1b73-4e2b-aa06-eda5a40979f7
solar_power_income_expanded(
	first(time_df_no_b.latitude),
	first(time_df_no_b.longitude), 
	first(time_df_no_b.altitude), 
	first(time_df_no_b.utc_time),
	first(short_track.diff_distance),
	first(speed_vec_no_b)
)

# ╔═╡ ff5c50ef-19df-4123-9da2-2efd43229b6f
solar_power_income_expanded.(
	time_df_no_b.latitude,
	time_df_no_b.longitude, 
	time_df_no_b.altitude, 
	time_df_no_b.utc_time,
	short_track.diff_distance,
	speed_vec_no_b
)

# ╔═╡ f46e89ce-dfa7-4eb4-91d4-6f8e97ab50ba
@time broadcasted_function_result = broadcast(solar_power_income_expanded, time_df_no_b.latitude,
	time_df_no_b.longitude, 
	time_df_no_b.altitude, 
	time_df_no_b.utc_time,
	short_track.diff_distance,
	speed_vec_no_b
)

# ╔═╡ 178b9ad6-7784-4984-9684-3241ba7217dd
@time broadcasted_function_result_dot = solar_power_income_expanded.(
	time_df_no_b.latitude,
	time_df_no_b.longitude, 
	time_df_no_b.altitude, 
	time_df_no_b.utc_time,
	short_track.diff_distance,
	speed_vec_no_b
)

# ╔═╡ cbe4aeeb-ad2a-4607-b010-20b34569773d
function solar_radiation_pvedication_alloc(data_df, track_dataframe)
    # starts here: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-outside-the-earths-atmosphere

    # copy generate_year_time_dataframe(100000)
    # data_df = copy(time_dataframe)
    data_df.latitude = track_dataframe.latitude
    data_df.longitude = track_dataframe.longitude
    data_df.altitude = track_dataframe.altitude

    # for debug purposes
    # data_df = copy(time_df)
    # data_df.latitude = track.latitude
    # data_df.longitude = track.longitude
    # data_df.altitude = track.altitude

    # setting up needed values for calculcations
    # data_df.day = Dates.day.(data_df.datetime)

    # solar constant
    H_constant = 1353 # W/m^2
    # AM0 = 1366 # air mass zero, W/m^2

    # radiant power density otside the Erath's athmosphere in W/m^2
    # H = H_constant * (1 .+ 0.033*cosd.( (360 * (Dates.dayofyear.(data_df.utc_time) .- 2)) / 365 ) )
    # plot(Dates.dayofyear.(data_df.utc_time),H, title = "Radial power density W/m^2")

    # TODO: air mass calculation with phi? and theta?

    # time shift for current latitude, hours
    # local_standard_time_meridian = 15 * Dates.Hour(local_time - UTC_time).value
    # very rough estimation
    # proper solution like https://stackoverflow.com/questions/5584602/determine-timezone-from-latitude-longitude-without-using-web-services-like-geona


    # equation of time
    B = 360 / 365 * (Dates.dayofyear.(data_df.utc_time) .- 81) # in degrees
    equation_of_time = 9.87 * sind.(2B) - 7.53 * cosd.(B) - 1.5 * sind.(B) # minutes
    # plot(
    #     equation_of_time,
    #     title = "Equation of time",
    #     label = "Equation of time",
    #     xlabel = "Day of the Year",
    #     ylabel = "Minutes"
    # )

    # TODO: resume from here
    # because local standard time meridian is not real sun time , minutes
    # time_correction_factor = 4 * (data_df.longitude - data_df.lstm) + equation_of_time
    # plot(time_correction_factor, title="Time correction factor")

    # local_solar_time = local_time + Dates.Second(round(time_correction_factor * 60.0))
    # # LST = 4 longitude + 60 UTC_minutes +EoT
    # # so, we only need to know UTC time and longitude
    # local_solar_time = UTC_time +
    #     Dates.Second( round( (4*longitude + equation_of_time[day]) * 60) )

    data_df.local_solar_time = data_df.utc_time +
        Dates.Second.( round.(4 * data_df.longitude + equation_of_time) * 60)
    
    # TODO: should hour angle be only in integer values? 
    # minutes_from_start_of_the_day = Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time);
    hour_angle = 15 * ((Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time)) ./60 .- 12)
    # TODO: hour angle calculation for UTC time
    # plot(data_df.utc_time, hour_angle, title="Hour angle vs UTC time")
    # TODO: hour angle calculation for local solar time ??? check if needed

    # TODO: continue to sun declination angle
    # TODO: continous sun declination angle calculation? (as for now quantified by days)
    sun_declination_angle = -23.45 * cosd.(360 / 365 * (Dates.dayofyear.(data_df.utc_time) .+ 10))
    # plot(sun_declination_angle, title = "Sun declination angle")

    # elevation (altitude) angle
    # elevation_angle = 90 .+ data_df.latitude .- sun_declination_angle
    # plot(data_df.utc_time, elevation_angle, title="Elevation angle (deg)")

    # elevation?
    elevation = asin.(
        sind.(sun_declination_angle) .* sind.(data_df.latitude) +
        cosd.(sun_declination_angle) .* cosd.(data_df.latitude) .* cosd.(hour_angle)
    )
    # plot(elevation, title="Elevation of the sun in rad")
    # when elevation angle is negative it means that sun is below horizon 
    # so, we should nullify these elements
    # elevation_filtered = copy(elevation)
	# elevation_filtered = elevation
    elevation[elevation .<=0] .= 0
    # plot(data_df.utc_time, elevation_filtered, title="elevation angle when sun above horizon")

    # TODO: sunrise and sunset hours

    # sunrise_hour = 12 .- 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60
    # # plot(data_df.utc_time, sunrise_hour, title = "sunrise hour")

    # sunset_hour = 12 .+ 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60

    # plot(data_df.utc_time, sunset_hour, title = "sunset hour")

    # calculating zenith angle, needed for optical air mass
    zenith_angle = pi/2 .- elevation
    # plot(data_df.utc_time, zenith_angle, title="Zenith angle of the sun in rad")
    cos_zenith_angle = cos.(zenith_angle)

    # solar radiation consists of direct and diffuse radiation
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass 

    # calculate air mass - how much air light must pass through
    air_mass_full = 1 ./ (cos_zenith_angle .+ 0.50572 .* (96.07995 .- deg2rad.(zenith_angle)).^-1.6364 )

    # direct intensity, not actually used
    # intensity_direct = H_constant .* 0.7 .^ (air_mass_full .^ 0.678) # W/m^2
    # H_constant * 70% of radiation ^ air_mass_full ^ coefficient to fit the experimental data
    # plot(data_df.utc_time, intensity_direct, title="Direct intensity at Earth ground, W/m^2")

    # with altitude (from sea level)
    intensity_direct_altitude = H_constant .* ( (1 .- 0.14 .* data_df.altitude ./ 1000) .* 0.7 .^ (air_mass_full .^ 0.678) .+ 0.14 .* data_df.altitude ./ 1000) # W/m^2
    # plot(data_df.utc_time, intensity_direct_altitude, title="Direct intensity at altitude, W/m^2")
    # diffuse radiation
    # intensity_diffuse = intensity_direct_altitude * 0.1
    # global irradiance
    intensity_global = intensity_direct_altitude * 1.1
    # plot(data_df.utc_time, intensity_global, title = "Global irradiance at altitude")
    # irradiance * cos(zenith_angle) is incident radiation ?

    
    # s_incident = intensity_global
    # TODO: calculate radiation on a tilted surface
    # can be done through perpendecular (incident) or horizontal
    # s_horizontal = s_incident .* sin.(elevation)
    # module_angle_rad - at what angle to surface panels are located
    # data_df.module_angle_rad .= 0.0 # just on the ground
    # data_df.azimuth_angle .= 0.0 # just something, since it is on the ground, will not be used
    # TODO: calculate sun azimuth angle according to https://www.pveducation.org/pvcdrom/properties-of-sunlight/azimuth-angle 
    # azimuth_cos = ( sind.(sun_declination_angle) .* cos.(data_df.latitude) .- 
    # cosd.(sun_declination_angle) .* sin.(data_df.latitude) .* cos.(hour_angle) ) ./ cos.(elevation)
    # azimuth_cos[azimuth_cos .> 1] .= 1
    # azimuth_cos[azimuth_cos .< -1] .= -1
    # data_df.azimuth_angle .= acos.( azimuth_cos )
    # data_df.lst .= Dates.hour.(data_df.local_solar_time)
    # data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ] .= 2*pi .- data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ]
    # module azimuth angle
    # module_azimuth_angle = data_df.azimuth_angle
    # s_module stand for solar module, tilted surface
    # s_module = s_incident .* sin.(elevation .+ data_df.module_angle_rad)
    # plot(data_df.utc_time, s_module, title="Solar intensity without module azimuth")
    # OR https://www.pveducation.org/pvcdrom/properties-of-sunlight/arbitrary-orientation-and-tilt 
    # s_module_azimuth = s_incident .* (
    #     cos.(elevation) .* sin.(data_df.module_angle_rad) .* cos.(data_df.azimuth_angle .- module_azimuth_angle) .+
    #     sin.(elevation) .* cos.(data_df.module_angle_rad)
    #     )
	s_module_azimuth = intensity_global .* sin.(elevation)
    # plot(data_df.utc_time, s_module_azimuth, title="Solar intensity with module azimuth")

    # next - simulate the race
    return s_module_azimuth
end


# ╔═╡ 45c19245-4dbb-4047-b0ed-13b27231c667
function solar_power_income_alloc(time_df, track_df, speed_vector)
    electrics_efficiency = 0.86
    solar_panels_efficiency = 0.228
    panels_area = 4 # m^2
    # solar_intensity = solar_radiation_pvedication_alloc(time_df, track_df)
    # power_income = electrics_efficiency .* solar_panels_efficiency .* solar_intensity .* panels_area .* track_df.diff_distance ./ speed_vector # W*s
    # power_income_wt_h = power_income ./ 3600
    # power_income(i) = electrics_efficiency * panels_efficiency * solar_rad_h * panels_area * ...
    # ( dist_diff(i) * 3.6 / ( speed * 3600 ) ); % kWt*h
    # return power_income_wt_h 
	return solar_radiation_pvedication_alloc(time_df, track_df) .* electrics_efficiency .* solar_panels_efficiency .* panels_area .* track_df.diff_distance ./ speed_vector ./ 3600
end

# ╔═╡ b1c4718b-3829-4b39-a84b-2445d605f08d
@time reduced_alloc_function_result = solar_power_income_alloc(time_df_no_b, short_track, speed_vec_no_b)

# ╔═╡ daee61c5-6a65-4688-8f1d-7ee6c5a8d884
function solar_radiation_no_broadcasting_arrays_alloc(latitude, longitude, altitude, utc_time)
    # starts here: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-outside-the-earths-atmosphere

    # copy generate_year_time_dataframe(100000)

    # for debug purposes
    # data_df = copy(time_df)
    # data_df.latitude = track.latitude
    # data_df.longitude = track.longitude
    # data_df.altitude = track.altitude

    # setting up needed values for calculcations
    # data_df.day = Dates.day.(data_df.datetime)

    # solar constant
    H_constant = 1353 # W/m^2
    # AM0 = 1366 # air mass zero, W/m^2

    # radiant power density otside the Erath's athmosphere in W/m^2
    # H = H_constant * (1 .+ 0.033*cosd.( (360 * (Dates.dayofyear.(data_df.utc_time) .- 2)) / 365 ) )
    # plot(Dates.dayofyear.(data_df.utc_time),H, title = "Radial power density W/m^2")

    # TODO: air mass calculation with phi? and theta?

    # time shift for current latitude, hours
    # local_standard_time_meridian = 15 * Dates.Hour(local_time - UTC_time).value
    # very rough estimation
    # proper solution like https://stackoverflow.com/questions/5584602/determine-timezone-from-latitude-longitude-without-using-web-services-like-geona
    # data_df.lstm = 15 * ceil(data_df.longitude * 24 / 360)

    # equation of time
    B = 360 / 365 * (Dates.dayofyear(utc_time) - 81) # in degrees
    equation_of_time = 9.87 * sind(2B) - 7.53 * cosd(B) - 1.5 * sind(B) # minutes
    # plot(
    #     equation_of_time,
    #     title = "Equation of time",
    #     label = "Equation of time",
    #     xlabel = "Day of the Year",
    #     ylabel = "Minutes"
    # )

    # TODO: resume from here
    # because local standard time meridian is not real sun time , minutes
    # time_correction_factor = 4 * (data_df.longitude - data_df.lstm) + equation_of_time
    # plot(time_correction_factor, title="Time correction factor")

    # local_solar_time = local_time + Dates.Second(round(time_correction_factor * 60.0))
    # # LST = 4 longitude + 60 UTC_minutes +EoT
    # # so, we only need to know UTC time and longitude
    # local_solar_time = UTC_time +
    #     Dates.Second( round( (4*longitude + equation_of_time[day]) * 60) )

    local_solar_time = utc_time +
        Dates.Second( round(4 * longitude + equation_of_time) * 60)
    
    # TODO: should hour angle be only in integer values? 
    # minutes_from_start_of_the_day = Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time);
    hour_angle = 15 * ((Dates.hour(utc_time) * 60 + Dates.minute(utc_time)) / 60 - 12)
    # TODO: hour angle calculation for UTC time
    # plot(data_df.utc_time, hour_angle, title="Hour angle vs UTC time")
    # TODO: hour angle calculation for local solar time ??? check if needed

    # TODO: continue to sun declination angle
    # TODO: continous sun declination angle calculation? (as for now quantified by days)
    sun_declination_angle = -23.45 * cosd(360 / 365 * (Dates.dayofyear(utc_time) + 10))
    # plot(sun_declination_angle, title = "Sun declination angle")

    # elevation (altitude) angle
    # elevation_angle = 90 .+ data_df.latitude .- sun_declination_angle
    # plot(data_df.utc_time, elevation_angle, title="Elevation angle (deg)")

    # elevation?
    elevation = asin(
        sind(sun_declination_angle) * sind(latitude) +
        cosd(sun_declination_angle) * cosd(latitude) * cosd.(hour_angle)
    )
    # plot(elevation, title="Elevation of the sun in rad")
    # when elevation angle is negative it means that sun is below horizon 
    # so, we should nullify these elements
    # elevation_filtered = copy(elevation)
	if elevation <= 0
		elevation = 0
	end
    # elevation_filtered[elevation_filtered .<=0] .= 0
	
    # plot(data_df.utc_time, elevation_filtered, title="elevation angle when sun above horizon")

    # TODO: sunrise and sunset hours

    # sunrise_hour = 12 .- 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60
    # # plot(data_df.utc_time, sunrise_hour, title = "sunrise hour")

    # sunset_hour = 12 .+ 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60

    # plot(data_df.utc_time, sunset_hour, title = "sunset hour")

    # calculating zenith angle, needed for optical air mass
    zenith_angle = pi/2 - elevation
    # plot(data_df.utc_time, zenith_angle, title="Zenith angle of the sun in rad")
    cos_zenith_angle = cos(zenith_angle)

    # solar radiation consists of direct and diffuse radiation
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass 

    # calculate air mass - how much air light must pass through
    air_mass_full = 1 / (cos_zenith_angle + 0.50572 * (96.07995 - deg2rad(zenith_angle))^-1.6364 )

    # direct intensity, not actually used
    # intensity_direct = H_constant .* 0.7 .^ (air_mass_full .^ 0.678) # W/m^2
    # H_constant * 70% of radiation ^ air_mass_full ^ coefficient to fit the experimental data
    # plot(data_df.utc_time, intensity_direct, title="Direct intensity at Earth ground, W/m^2")

    # with altitude (from sea level)
    intensity_direct_altitude = H_constant * ( (1 - 0.14 * altitude / 1000) * 0.7 ^ (air_mass_full ^ 0.678) + 0.14 * altitude / 1000) # W/m^2
    # plot(data_df.utc_time, intensity_direct_altitude, title="Direct intensity at altitude, W/m^2")
    # diffuse radiation
    # intensity_diffuse = intensity_direct_altitude * 0.1
    # global irradiance
    intensity_global = intensity_direct_altitude * 1.1
    # plot(data_df.utc_time, intensity_global, title = "Global irradiance at altitude")
    # irradiance * cos(zenith_angle) is incident radiation ?
    s_module_azimuth = intensity_global * sin(elevation)
    # plot(data_df.utc_time, s_module_azimuth, title="Solar intensity with module azimuth")

    # next - simulate the race
    return s_module_azimuth
end

# ╔═╡ cf80712e-bf9a-4fde-a058-c0404011dca2
function solar_power_income_expanded_alloc(latitude, longitude, altitude, utc_time, diff_distance, speed_vector)
    electrics_efficiency = 0.86
    solar_panels_efficiency = 0.228
    panels_area = 4 # m^2
    # solar_intensity = solar_radiation_no_broadcasting_arrays_alloc(latitude, longitude, altitude, utc_time)
    # power_income = electrics_efficiency * solar_panels_efficiency * solar_intensity * panels_area * diff_distance / speed_vector # W*s
    # power_income_wt_h = power_income / 3600
    # # power_income(i) = electrics_efficiency * panels_efficiency * solar_rad_h * panels_area * ...
    # # ( dist_diff(i) * 3.6 / ( speed * 3600 ) ); % kWt*h
    # return power_income_wt_h 
	return solar_radiation_no_broadcasting_arrays_alloc(latitude, longitude, altitude, utc_time) * electrics_efficiency * solar_panels_efficiency * panels_area * diff_distance / speed_vector / 3600.0
end

# ╔═╡ ad211fed-2077-4869-bdcb-23c9804db113
@time broadcasted_function_result_dot_alloc = solar_power_income_expanded_alloc.(
	time_df_no_b.latitude,
	time_df_no_b.longitude, 
	time_df_no_b.altitude, 
	time_df_no_b.utc_time,
	short_track.diff_distance,
	speed_vec_no_b
)

# ╔═╡ 6038ed31-c1c5-4a3a-8809-9e6a39fd0e5b
function mechanical_power_calculation_alloc(speed_ms, slope, diff_distance)
    drag = 0.18
    frontal_area = 1 # m^2
    ro = 1.18 # air density

    mass = 390 # kg
    g = 9.8019 # at start: [41.2646201567207,-95.9244249307473,301.540649414063];
    # 9.80147 at finish: [43.9660024736000,-121.345052439700,1229.07763671875]
    friction_1 = 0.0023;
    friction_2 = 0.000041; # total friction = friction_1 + friction_2*speed

    engine_efficiency = 0.87

    # mechanical force = drag force + friction force + gravitational force
    # newtons
    mechanical_force = (
        drag * frontal_area * speed_ms ^ 2 * ro / 2. +
        mass * g * (friction_1 + friction_2 * 4 * speed_ms) * cosd(slope) +
        mass * g * sind(slope)
        )

    # mechanical power = mechanical force * distance delta / engine efficiency
    # watts * s
    mechanical_power = (
        mechanical_force * diff_distance / engine_efficiency
        )

    # TODO: get rid of return, or at least make it type-stable
    # see https://docs.julialang.org/en/v1/manual/faq/#Types,-type-declarations,-and-constructors-1
    return mechanical_power
end

# ╔═╡ c4d4fb36-fcd4-4d1c-b5f3-f208ad48fd3f
@time mechanical_power_calculation_alloc.(speed_vector, track.slope, track.diff_distance)

# ╔═╡ 3028c7b7-c633-47b9-bb5d-f15101375b4f
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
    # power_use_accumulated_wt_h = calculate_power_use_accumulated(mechanical_power, electrical_power)
	power_use_accumulated_wt_h = mechanical_power + electrical_power
	# cumsum!(power_use_accumulated_wt_h, power_use_accumulated_wt_h)
	power_use_accumulated_wt_h_cumsum = cumsum(power_use_accumulated_wt_h)
	power_use_accumulated_wt_h_cumsum = power_use_accumulated_wt_h_cumsum / 3600.

    # get solar energy income
	# @debug "track size is $(size(track.latitude, 1))"
    solar_power = solar_power_income_expanded_alloc.(
		track.latitude,
		track.longitude, 
		track.altitude, 
		time_df.utc_time,
		track.diff_distance,
		input_speed
	)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power)
    # TODO: night charging with additional solar panels

    # #### plotting
    # plot(track.distance, power_use, title="Power spent on toute")
    # plot(track.distance, solar_power, title="Power gained on the route")

    # plot(track.distance, power_use_accumulated_wt_h, title="Power spent on the route, accumulated")
    # plot(track.distance, solar_power_accumulated, title="Power gained on the route, accumulated")

    # plot(track.distance, solar_power_accumulated - power_use_accumulated_wt_h, title="Power balance w/o battery")
    # battery_capacity = 5100 # wt, to be used later for physical constraints
    energy_in_system = start_energy .+ solar_power_accumulated .- power_use_accumulated_wt_h_cumsum
    # plot(track.distance, energy_in_system, title="Power balance with battery")

    # TODO: calculate night charging - do it later since it is not critical as of right now
    # TODO: block overcharging - cost function?
    # at first do the black-box optimization, then gradient one
    # will start with Optim
    # TODO: find an optimal single speed - make a loss function and start optimization process
    # time_seconds = calculate_travel_time_seconds(input_speed, track)
	time_seconds = track.diff_distance ./ input_speed
	# cumsum!(time_seconds, time_seconds)
	time_seconds_cumsum = cumsum(time_seconds)
    # TODO: find an optimal speed vector
    return power_use_accumulated_wt_h_cumsum, solar_power_accumulated, energy_in_system, time_df, time_seconds_cumsum
end

# ╔═╡ 49158b6b-be82-4689-9a7b-4de7e707e08c
@time power_use_test_alloc, solar_power_test_alloc, energy_in_system_test_alloc, time_test_alloc, time_s_test_alloc = solar_trip_calculation_bounds_alloc(speed_vector, track, DateTime(2022,7,1,0,0,0), 5099.0);

# ╔═╡ 54daafd9-2842-43f4-bde3-eabdd6403767
function solar_partial_trip_wrapper_alloc(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds(speeds_ms, track, indexes)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds_alloc(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) + 10 * abs(finish_energy - last(energy_in_system))
	return cost
	# return solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
end

# ╔═╡ 2a00bd6b-a0ce-4b54-bf88-d6cc2b2204c8
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
	line_search = LineSearches.BackTracking();
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
	# println("time is $(last(time_s)), cost is $(f(minimized_speeds))")
	# println("split indexes are $(split_indexes)")
	# # println("distances are $(track[split_indexes, :])")
	# # println("minimized speeds are: $(minimized_speeds)")
	# println("simulated finish energy is $(last(energy_in_system))")
	# # println("calculated cost is $( last(time_s) + 100 * abs(last(energy_in_system) - finish_energy) + 100 * sum(abs.(energy_in_system[energy_in_system .< 0.0])) + 100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .< 0.0])) + 100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .> 100.0])) )")
	# # println("finish energy difference penalty is: $(100 * abs(last(energy_in_system) - finish_energy))")
	# # println("energy less than 0. penalty is: $(100 * sum(abs.(energy_in_system[energy_in_system .< 0.0])))")
	# # println("speed less than 0. penalty is: $(100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .< 0.0])))")
	# # println("speed more than 100. penalty is: $(100 * sum(abs.(minimized_speed_vector[minimized_speed_vector .> 100.0 / 3.6])))")
	split_energies = energy_in_system[split_indexes]
	pushfirst!(split_energies, start_energy)
	split_times = time[split_indexes, :utc_time]
	pushfirst!(split_times, start_datetime)
	# println("split energies are $(split_energies)")
	# println("")
	
	# 5 - go though each track piece and enter function again
	# hierarchical_optimization(minimized_speeds[1], tracks[1], chunks_amount, start_energy, split_energies[1], start_datetime, iteration + 1)
	# @debug "split_energies size is $(size(split_energies, 1)), chunks_amount is $(chunks_amount)"
	result_speeds = []
	for i=1:chunks_amount
		result_speeds_chunk = hierarchical_optimization_alloc(minimized_speeds[i], tracks[i], chunks_amount, split_energies[i], split_energies[i+1], split_times[i], iteration + 1 , start_index + split_indexes[i] - 1)
		append!(result_speeds, result_speeds_chunk)
	end

	return result_speeds
end

# ╔═╡ 28f216c0-f5b3-4f22-a5eb-192423b980e9
@time result_hierarchical_alloc = hierarchical_optimization_alloc(initial_speed, short_track, chunks_amount_hierarchical, start_energy_short, 0., start_datetime_hierarchical, 1, track_size)

# ╔═╡ 3198a48f-36e2-452a-9acd-e478a15057a3
begin
	inputs_ms_hier_alloc = abs.(convert_kmh_to_ms(result_hierarchical_alloc))
	power_use_hier_alloc, solar_power_hier_alloc, energy_in_system_hier_alloc, time_hier_alloc, time_s_hier_alloc = solar_trip_calculation_bounds_alloc(inputs_ms_hier_alloc, short_track, start_datetime_hierarchical, start_energy_short)
	last(time_s_hier_alloc)
end

# ╔═╡ a80b460a-efac-4ecc-ad30-3a00e200f31e
last(time_s_hier_alloc) - last(time_s_hier)

# ╔═╡ 9e921706-7eda-4574-b8c5-ad368490b4e1
plot(short_track.distance, [power_use_hier_alloc solar_power_hier_alloc energy_in_system_hier_alloc zeros(track_size)],
	    label=["Energy use" "Energy income" "Energy in system" "Failure threshold"], title="Energy graph (distance) for short track Hierarchical",
	    xlabel="Distance (m)", ylabel="Energy (W*h)", lw=3, #size=(1200, 500),
	    color=[:blue :green :cyan :red] # ,ylims=[-10000, 40000]
)

# ╔═╡ e13c04ce-bb90-480f-978c-2f02f01d31de
plot(short_track.distance, power_use_hier_alloc - power_use_hier)

# ╔═╡ bd7d392b-f98f-4e98-9a5d-155d1d6fd199
begin
	plot(short_track.distance, short_track.altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance (alloc)", right_margin = 15Plots.mm)
	plot!(twinx(), short_track.distance, result_hierarchical_alloc, color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance (alloc)")
end

# ╔═╡ 67791160-75c3-476f-b7e0-6d264a9dc1a4
begin
	plot(short_track[1000:1050,:].distance, short_track[1000:1050,:].altitude, label="altitude", ylabel="altitude", title="Speed (km/h) vs distance [1000:1050]", right_margin = 15Plots.mm)
	plot!(twinx(), short_track[1000:1050,:].distance, result_hierarchical_alloc[1000:1050], color=:red, ylabel="speed (km/h)", label="speed (km/h)", ymirror = true, title="Speed (km/h) vs distance [1000:1050]")
end

# ╔═╡ fdec6c16-d095-4914-86bb-b94983b04465
plot(short_track.distance, result_hierarchical_alloc - result_hierarchical)

# ╔═╡ a81bfab4-dc7a-4fff-a72f-99a8dfd73a5d
function solar_trip_calculation_bounds_alloc_reverse(input_speed, track, start_datetime,
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
    # power_use_accumulated_wt_h = calculate_power_use_accumulated(mechanical_power, electrical_power)
	power_use_accumulated_wt_h = +(mechanical_power, electrical_power)
	# cumsum!(power_use_accumulated_wt_h, power_use_accumulated_wt_h)
	power_use_accumulated_wt_h_cumsum = cumsum(power_use_accumulated_wt_h)
	power_use_accumulated_wt_h_cumsum = power_use_accumulated_wt_h_cumsum / 3600.

    # get solar energy income
	# @debug "track size is $(size(track.latitude, 1))"
    solar_power = solar_power_income_expanded_alloc.(
		track.latitude,
		track.longitude, 
		track.altitude, 
		time_df.utc_time,
		track.diff_distance,
		input_speed
	)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power)
    # TODO: night charging with additional solar panels

    # #### plotting
    # plot(track.distance, power_use, title="Power spent on toute")
    # plot(track.distance, solar_power, title="Power gained on the route")

    # plot(track.distance, power_use_accumulated_wt_h, title="Power spent on the route, accumulated")
    # plot(track.distance, solar_power_accumulated, title="Power gained on the route, accumulated")

    # plot(track.distance, solar_power_accumulated - power_use_accumulated_wt_h, title="Power balance w/o battery")
    # battery_capacity = 5100 # wt, to be used later for physical constraints
    energy_in_system = start_energy .+ solar_power_accumulated .- power_use_accumulated_wt_h_cumsum
    # plot(track.distance, energy_in_system, title="Power balance with battery")

    # TODO: calculate night charging - do it later since it is not critical as of right now
    # TODO: block overcharging - cost function?
    # at first do the black-box optimization, then gradient one
    # will start with Optim
    # TODO: find an optimal single speed - make a loss function and start optimization process
    # time_seconds = calculate_travel_time_seconds(input_speed, track)
	time_seconds = track.diff_distance ./ input_speed
	# cumsum!(time_seconds, time_seconds)
	time_seconds_cumsum = cumsum(time_seconds)
    # TODO: find an optimal speed vector
    return power_use_accumulated_wt_h_cumsum, solar_power_accumulated, energy_in_system, time_df, time_seconds_cumsum
end

# ╔═╡ 72160b43-5668-49f8-81cc-076f3ab5b300
function solar_partial_trip_grad_wrapper_alloc_reverse(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds_grad(speeds_ms, track, indexes)

	# power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds(speed_vector, track, start_datetime, start_energy)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds_alloc_reverse(speed_vector, track, start_datetime, start_energy)

	cost = last(time_s) + 10 * abs(finish_energy - last(energy_in_system))
	# return cost
	
	# cost = last(time_s) + 10 * abs(last(energy_in_system) - finish_energy) + sum(energy_in_system[energy_in_system .< 0.0])
	
	# + 100 * sum(abs.(speed_vector[speed_vector .< 0.0])) + 100 * sum(abs.(speed_vector[speed_vector .> 100.0 / 3.6]))
	# + 100 * abs(minimum(energy_in_system) - finish_energy) 
	
	# + 100 * sum(speed_vector[speed_vector .> 200.0]) + 100 * sum(energy_in_system[energy_in_system .< 0.0])


	# cost = last(time_s) + 10 * abs(minimum(energy_in_system)) + 100 * sum(input_speed[input_speed .< 0.0])
	return cost
end

# ╔═╡ 0cc7f382-5424-4205-8b4f-5e1fc5e3bae3
function f_reverse(speed)
	# return solar_partial_trip_wrapper_alloc(speed, track, split_indexes, start_energy, finish_energy, start_datetime)
	return solar_partial_trip_grad_wrapper_alloc_reverse(speed, track[25:27,:], [1, 2, 3], 1.5990611991330004, 1.1239726167707227, DateTime(2022,7,1,16,0,8))
end

# ╔═╡ 01d54fa2-cd95-4150-af88-f1fbf8f1ae08
ftape_vec = ReverseDiff.GradientTape(f_reverse, rand(3))
# not differiantiable - see last point in https://juliadiff.org/ReverseDiff.jl/limits/
# revert to simpler methods?

# ╔═╡ 990755c3-7cd5-49f9-9ddd-b86d6b267cb5
compiled_tape = ReverseDiff.compile(ftape_vec)

# ╔═╡ ed32a5b0-a5ee-4eec-9093-657ce7f907a0
f_reverse(inputs_reverse)

# ╔═╡ 7ac8a793-352b-4b73-aa3d-d93b0c7fce9e
@time ReverseDiff.gradient(f_reverse, inputs_reverse)

# ╔═╡ 4894c24f-9729-4ebc-80a9-72f6e0763a4b
function f_reverse2(speed)
	# return solar_partial_trip_wrapper_alloc(speed, track, split_indexes, start_energy, finish_energy, start_datetime)
	return solar_partial_trip_grad_wrapper_alloc_reverse(speed, short_track, split_indexes_reverse, start_energy_short, 0., DateTime(2022,7,1,16,0,8))
end

# ╔═╡ d33d71f3-e16c-4001-b561-239e00aaa102
ftape_vec2 = ReverseDiff.GradientTape(f_reverse2, inputs_reverse2)

# ╔═╡ 7c982555-a277-49bf-add5-1b27ec066d10
compiled_tape2 = ReverseDiff.compile(ftape_vec2)

# ╔═╡ e4441a0e-55ad-47bd-8b30-004417d9dc71
ftape_jac = ReverseDiff.JacobianTape(f_reverse2, inputs_reverse2)

# ╔═╡ db3eedbf-ad89-487f-a8f9-f018f52cd2f0
ftape_hess = ReverseDiff.HessianTape(f_reverse2, inputs_reverse2)

# ╔═╡ 61b50198-e582-4171-8752-7108f5021d30
@time ReverseDiff.gradient!(result_reverse2, compiled_tape2, inputs_reverse2)

# ╔═╡ 1f56f505-ba6f-4f59-9985-ca844874e336


# ╔═╡ 04587beb-c645-49cb-b74e-bfa51e434fdc
function grad!(storage, inputs)
	ReverseDiff.gradient!(storage, compiled_tape2, inputs)
end

# ╔═╡ 496eba70-c366-4ecf-99f1-850a0bb64ecf
function hess!(storage, inputs)
	ReverseDiff.hessian!(storage, compiled_tape2, inputs)
end

# ╔═╡ b829bbdd-84ef-4c66-af88-fd928518a73f
@time result_chunks_reverse = optimize(f_reverse2, grad!, hess!, fill(40., chunks_amount_hierarchical), Newton(),
	Optim.Options(
		x_tol = 1e-6,
		f_tol = 1e-8,
		g_tol = 1e-6
	)
)

# ╔═╡ 08fd7243-13d8-47af-b4d1-2aa410927e93
@md_str "## Not really working because of cumsum"

# ╔═╡ 6d29b8d8-bb2e-4fc6-ae72-9f97265f1451
@md_str "## WORKING!!!

TODO: trying out hierarchical optim with ReverseDiff"

# ╔═╡ 6e06a810-186a-4d91-b63a-259c4ef32300


# ╔═╡ dffe90b9-3b9a-4f9c-a7c7-365d894777ea


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlotlyBase = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
TimeZones = "f269a46b-ccf7-5d73-abea-4c690281aa53"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.4.1"
LineSearches = "~7.2.0"
Optim = "~1.7.3"
Peaks = "~0.4.1"
PlotlyBase = "~0.8.19"
Plots = "~1.35.4"
PlutoUI = "~0.7.48"
ReverseDiff = "~1.14.4"
TimeZones = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "4a00228acf5c4ed037e08a4e454ed29f6326b0e6"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e9f7992287edfc27b3cbe0046c544bace004ca5b"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.22"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "84259bb6172806304b9101094a7cc4bc6f56dbc6"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "46d2680e618f8abd007bce0c3026cb0c4a8f2032"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.12.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "558078b0b78278683a7445c626ee78c86b9bb000"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "8b7a4d23e22f5d44883671da70865ca98f2ebf9d"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "5a2cff9b6b77b33b89f3d97a4d367747adce647e"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.15.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "187198a4ed8ccd7b5d99c41b69c679269ea2b2d4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.32"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "00a9d4abadc05b9476e937a5557fcce476b9e547"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.69.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "3cdd8948c55d8b53b5323f23c9581555dc2e30e1"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.5.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "db619c421554e1e7e07491b85a8f4b96b3f04ca0"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "6872f9594ff273da6d13c7c1a1545d5a8c7d0c1c"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.6"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "3c3c4a401d267b04942545b1e964a20279587fd7"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "b9fe76d1a39807fdcf790b991981a922de0c3050"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.3"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

[[deps.Peaks]]
deps = ["Compat"]
git-tree-sha1 = "5f1390b0a0ef6d6411f9a9a37c4444d6a7e44780"
uuid = "18e31ff7-3703-566c-8e60-38913d67486b"
version = "0.4.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "041704a5182f25cdcbb1369f13d9d9f94a86b5fd"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.35.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "460d9e154365e058c4d886f6f7d6df5ffa1ea80e"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.1.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "d12e612bba40d189cead6ff857ddb67bd2e6a387"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "9b1c0c8e9188950e66fc28f40bfe0f8aac311fe0"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.7"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ReverseDiff]]
deps = ["ChainRulesCore", "DiffResults", "DiffRules", "ForwardDiff", "FunctionWrappers", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NaNMath", "Random", "SpecialFunctions", "StaticArrays", "Statistics"]
git-tree-sha1 = "afc870db2b2c2df1ba3f7b199278bb071e4f6f90"
uuid = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
version = "1.14.4"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "f86b3a049e5d05227b10e15dbb315c5b90f14988"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.9"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Scratch", "Unicode"]
git-tree-sha1 = "d634a3641062c040fc8a7e2a3ea17661cc159688"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.9.0"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═fb9bd7c0-4f3c-11ed-094a-35559b7aedad
# ╠═7cc7b885-f841-4bcc-a82c-2f947c74de22
# ╠═90cd5c2a-223c-4345-8f78-498973c0fd46
# ╠═a7781b82-48db-4954-8e79-ab8e7864ed69
# ╠═44f352a0-0c5b-41f6-a21d-56f10edae4c9
# ╠═9d915c29-314b-47f9-9e4c-898ebd28f88a
# ╠═92e47c4c-2e70-4850-af86-a997fcce5587
# ╟─9f2624cc-398b-42d2-bf7f-527285af6dc3
# ╠═704c9c6a-8b13-4e82-92fd-edda72397320
# ╟─afe9cc57-f93c-4780-a592-b3d2609162f2
# ╟─4e4c4794-aa95-49e9-961b-ed7c4bb81442
# ╟─e80aee02-b99e-44c9-b503-9443c431b0e6
# ╟─ae4d73ea-bd18-472b-babb-7980598a4ce9
# ╟─b9bdb969-2f65-47d8-b1f2-a9b7e411f1c1
# ╟─36e15048-b615-4ba6-a32a-093922a49243
# ╟─1536f51f-0776-4d87-832d-e9b2b9cc36d6
# ╟─ca2f1ec3-4ed3-4b5d-bfcd-ab43de0d2abc
# ╟─4f286021-da6e-4037-80e3-a526e880b686
# ╟─1ba728ba-b6aa-40c9-b5a2-906297bd4921
# ╟─5674f163-5847-4e2b-ba3f-c95348e8d1d5
# ╟─c3d41246-37cb-45f6-ac98-f64d9158087a
# ╟─4c91701c-f90e-48cf-bb9d-d925baa85667
# ╟─9977c5d5-4872-41d4-8e99-1e4034493e4d
# ╠═422a0d48-40fe-41eb-b214-e21d009c00b2
# ╠═6a5416d0-39c0-431b-8add-5dbf13a1bda0
# ╠═bef6c840-8dc7-4839-b2ba-623c6c46c856
# ╟─54a5bf51-1723-4ce1-8f6b-1fd199c991b5
# ╠═c7434aee-3089-4098-8c9d-2d13c9df5ee9
# ╠═8794ae20-fe98-47ab-bd80-681c09edb7d1
# ╠═7380a326-1e9e-437c-9cb7-3aa2b54b8ec5
# ╟─596e51be-969e-4fe6-8170-22c1aac89aca
# ╠═eeb85fd3-6721-4e0f-aed9-73731960ac35
# ╠═09bd67a1-a0f9-45f3-9839-2e9e592f01de
# ╠═3ca6f786-8add-4c46-b82a-30a570828d39
# ╠═e50a7ae9-a46e-41b0-8a10-d77e9ffa7b14
# ╠═5c9c006c-f814-4829-8c18-108546be870b
# ╠═264d9d64-6a6d-4e84-9af1-5795bd5bf829
# ╠═411e63ec-b83a-4e21-9535-5d0275381039
# ╠═21634b70-7b3a-44b2-b629-01664ce81acf
# ╠═5f3a7fcf-e261-4f64-a94c-57a12093e353
# ╠═de201868-7805-4f27-81b7-f4f8204eface
# ╠═96a68dec-d781-4fd6-8146-649434f60919
# ╠═77d82639-dd61-46e0-b6a0-7c7400a10453
# ╠═655ff7ad-d7fc-47d4-bd22-0bb2c4b63cd5
# ╠═da6bb535-9601-4cf9-ba64-f8bfc68f3e5d
# ╠═4f9da3c2-36e1-405a-9c83-edce6f9a30af
# ╠═09b71f55-2e31-4a7f-a6ff-f4a3ff38e4b9
# ╠═a5551c27-3f5c-4ea3-8c4e-9d4873c88751
# ╠═e7445b07-20e8-404e-8d40-eccb9eb6ddc2
# ╠═b6a95a4d-6f5d-4fbc-87ed-622768e8d47a
# ╠═7be923ff-1be1-4496-a571-b3dd68a03ebc
# ╠═d240b89b-2d6c-4379-8132-1a47ec0ef6f3
# ╠═2b239f4b-7149-4771-8498-35b0314fb928
# ╠═b57ed73a-a69e-4944-b6f2-2e6bcafc2442
# ╠═07cd4d1e-a293-4aa4-85fc-43cc26be9f29
# ╠═9f924c63-9558-4915-ab1a-f8a383795906
# ╠═7779fe61-ca38-4bb6-aa6e-e7bc84870fe4
# ╠═249b0a50-0713-474f-a2ce-f7c012522318
# ╠═899c3a73-28ef-46e3-b08b-d2335d91f52d
# ╠═55157e96-b0b0-40fd-9db1-52edf69b4001
# ╠═d2b833e2-fd0a-4ae9-8f98-eb03f658339e
# ╠═848b59f6-d3f4-49c6-af22-7a40f5449071
# ╠═35285182-9a52-4e37-959a-09b41927d56f
# ╠═cec4a60b-07f6-429a-8c15-6211ba3ded74
# ╠═21a36ea3-d1d8-4a64-bb19-12cb85df5da7
# ╠═f3e06f95-0dd8-4c48-bbdc-fdac37be6548
# ╠═805f1520-28b0-45dd-80a5-470a33bacdee
# ╠═294fa952-015d-4709-bb55-c17682ffe2fb
# ╠═f0206c0b-e35f-499b-8aa0-32c69541b2e3
# ╠═7fc0e776-c9cf-4736-b113-81c31ada99c9
# ╠═ff9dad0d-4a00-439c-b704-486e395e5997
# ╠═2b9b2782-09d8-4cfa-99fa-9c2e921efe36
# ╠═08f34cff-83b5-4e0b-a789-3c180c8d5cea
# ╠═119c4024-6e50-4d20-a0b2-734ec9bf515d
# ╠═f243f6d4-c1b8-4284-abed-d74ae47a7af5
# ╠═6baaa312-94f8-4fc3-923f-30033a4a75d4
# ╠═360727dd-c241-4609-85fb-5f40553c1d2b
# ╠═cb7bf1ac-4abd-4c08-8893-e367cda83f8f
# ╠═141378f3-9b42-4b83-a87c-44b3ba3928aa
# ╠═9a6f42a2-7a9f-41ef-8ff2-b325a5971e42
# ╠═516e68a9-fb52-475c-b925-bba877341499
# ╠═9384af58-d6ad-4abd-8df0-f4e0d5649e68
# ╠═7ed6b68e-2cb3-4bd0-8c98-55ee42349d93
# ╠═f0fdeb9c-8791-4119-acc6-a1ec9094c590
# ╠═a6c2c405-08e9-43ba-843a-d4bd85ded0c8
# ╠═9db98630-5796-4763-92e8-d04a4ad3845a
# ╠═e9e72767-4fd3-4e9f-97bc-01c1f35ec916
# ╠═c1c47be6-b0a4-41bf-a284-26f193792748
# ╠═eb193e5d-65ba-46b1-be15-5c399abad44b
# ╠═64ec134c-e5cf-4c97-869f-d39b91e2599e
# ╠═ca583c57-d93c-4bb3-8cd3-b92184226a5f
# ╠═e75a6ae6-a09c-4c21-a421-d0124dd355c6
# ╠═19d51114-43d4-4c38-9a3d-55ec982f7c56
# ╠═2de331ba-bd7d-49bd-80a2-60270c769a7c
# ╠═724deff0-8265-4ca9-aa3e-2dfdd6f4d293
# ╠═99a41fe8-cc03-49fe-a39f-c27070fda62e
# ╠═96ee2227-0386-4a86-93ad-fa41dcdbf615
# ╠═56f144c3-43ed-46d5-9e55-fd800bd24e9e
# ╠═7699d86d-4494-494a-b367-7c3f13fc0022
# ╠═38142a0b-1e1f-4cda-90e9-50d52e6f0227
# ╠═f230c95a-d02a-436b-b426-364fe112cbba
# ╠═9d53a659-612c-405f-88da-253ffe57f4a0
# ╠═87b00642-3d8a-4722-bb7b-b9ff12a26d8f
# ╠═f6ed6270-5871-4dec-b2a3-5fa78caff34b
# ╠═8d2d34de-64db-49a8-b874-cb5af3bcf6cd
# ╟─79afe017-7206-4c35-b5c1-84a6e4cc3517
# ╟─c1a6f532-e95b-4fc3-b419-255038ee3589
# ╟─32bad172-97e7-485e-b64a-3c139400b471
# ╠═100ab330-4b6c-4033-92ee-9d19b2f7d2e4
# ╠═cc75dbac-9f9c-4465-a08b-e543d82021fd
# ╠═0c127ac8-812e-4a99-9e0e-32b3fd316c26
# ╠═d9e30de2-e75f-423b-8fcc-ab3847331274
# ╠═b55c819b-f312-4078-b751-cf443355be19
# ╠═a8b6ce9d-5507-4acb-954f-b43d702e1060
# ╠═3642f56d-ac97-435a-b446-68eb7814b03c
# ╠═6732e4c8-cd5f-454c-84fb-14aae6c02fbe
# ╠═f0d3e885-68c3-424e-8a50-d3d981bef295
# ╠═db7cc403-40ba-47d3-b27b-2a5913633ae5
# ╠═e1c90043-eae2-40fd-afc4-e5cd0ab4a540
# ╠═bb76116a-5399-4d95-9e24-a158cd438619
# ╠═d8e62070-5595-4e5f-b066-77f9e7d816b1
# ╠═e7b301c0-7c5d-4d45-a70f-91653e81ea0d
# ╠═8635aede-c512-470a-b1d0-f53d783c6179
# ╠═9a73fa1f-ac16-420d-8c98-91d1aeba773a
# ╠═785c4601-6d60-49b6-9470-abe90eed9d4c
# ╠═802e662b-0922-4993-bcf6-059c34ce02b7
# ╠═d5a4ad7a-2d5c-4597-947e-3f2e6b9e67e9
# ╠═55a9fc54-39b7-4c40-b06e-eba6ba260c50
# ╠═ce217af4-f88e-4767-984b-557eaf007bfb
# ╠═0a54b576-9d04-42ee-a19f-5ecafed7b2d6
# ╠═ef772852-003d-4c0c-b360-62941a8c8357
# ╠═c9c142d9-2bfb-4722-8d7e-a8c6724f3351
# ╠═bf7b5b54-2f2a-4f5c-a877-1fd8f53510b0
# ╠═d75bf2cc-5a06-4cd6-b4ce-ed0ce415fd79
# ╠═c7197a1e-d135-46e5-88b1-771934d58bf5
# ╠═63245989-6885-4d7c-bedd-cece121bbdef
# ╠═25f0fcbb-382e-4254-b058-52af29f212f9
# ╠═7d514799-a9b2-482f-aa16-c7d4a9fbfe78
# ╠═9ad839d3-40cd-438b-9a93-aaa0635ec0a1
# ╠═32160090-9fef-4216-a160-76a0f0af0f0f
# ╠═8d6a6e0f-044b-4c6e-8432-b5da6817d019
# ╠═5d4ac822-380e-477d-b8a3-7d6f39b0d9ec
# ╠═2c6fbc64-c7f3-430d-b6a8-c4a78f49c9c4
# ╠═adc5d37e-5f66-4d11-8ba3-1426bb823230
# ╠═f30abea4-a1d3-40fe-8328-3a0e5ce8a0d9
# ╠═8a3b49a9-4472-4b4d-944c-6b3e92a47b9f
# ╠═65f69edd-24a7-4a21-91b3-2e86c17df1a1
# ╠═f7b832ee-c8d5-4144-bd18-a68ed34c3e6a
# ╠═2e5f895b-c8b4-4bbb-8e69-88c554022be6
# ╠═d08e3a6e-9ca6-4ab9-8228-5a24546e1d36
# ╠═2143a79c-efe9-47df-88c8-eda9c2e16623
# ╠═5629455d-a6ce-430d-862d-5467ffc50ac9
# ╠═8f29c327-a738-4599-be47-4ff868b303e6
# ╠═51765acb-8459-426c-9c3a-773782f9789f
# ╠═c8a9d213-b294-4c46-9af3-3682ea171766
# ╠═bb508865-1543-4621-abfd-f79a90f80db6
# ╠═981ee4ea-1a83-44de-81da-6e472496aa7b
# ╠═93aafc03-dd14-4ca8-9226-5bbf8d42f378
# ╠═c8941591-bacb-469e-8d7f-cf15d32169d3
# ╠═5bc52c1f-bf9f-447c-a7d8-795d374ca4f4
# ╠═a13b519f-10e2-49eb-9903-8f198676cb76
# ╠═a92a0f58-0462-4a07-9053-a89caf5be9bb
# ╠═840dd334-16ad-4645-9773-393497547ae6
# ╠═22d50e51-2cee-474b-87f9-28f73cfd25b1
# ╠═f4c411a4-1c12-42b6-aa9b-030d7afdccc4
# ╠═6187da23-54e1-475d-9afd-73cd90600088
# ╠═71291ca2-6ab5-43a8-9f54-58990e4dcbb3
# ╠═233a5ed0-00eb-4b4e-9225-5838ceddb187
# ╠═e7eafaf3-35be-43ab-a8d5-30bff227664b
# ╠═2f20704a-f2e0-4843-b458-6c3208758ab6
# ╠═aa00efad-65b3-449d-9085-da6cf66c57c6
# ╠═b7e2d987-463b-4d51-ba98-6cdbd3fbb01c
# ╠═1c14c9e2-a05c-4301-bfe6-69f36d0be865
# ╠═c8d5a429-ea0e-4730-bc15-4b1e49642407
# ╠═9f247245-a14d-42ce-97e7-a60911988d9c
# ╠═39b4c2a0-fb33-4820-87f3-fa575ebb4e30
# ╠═e983aeb2-38c2-4bc4-af61-8af08f2347f5
# ╠═3c7d8d1d-f975-4743-99ca-263ac269fcab
# ╠═97bbd950-40f6-4026-afc0-eecef2d0c784
# ╠═16d51785-7cc0-4943-ba5a-0f0241315cfd
# ╠═bb69fdb6-b6b8-42ab-9a27-0f926067c1a6
# ╠═569ab1f6-c5aa-48d8-88ae-1884ca22518a
# ╠═aab805b0-2fe4-4ece-8f8c-f922226a9912
# ╠═48e146a3-861b-4005-acf1-e877dbd83a50
# ╠═337efda5-da6c-40c2-97a7-486ba07d2175
# ╠═7107b365-28d2-4ec7-bea7-a93b5c89a9d8
# ╠═710b23ba-7648-4795-96e4-8141766479d5
# ╠═f744cc32-907f-40b0-8586-20af362b2dc9
# ╠═873458fc-f7bd-4154-838d-0f65461f0178
# ╠═70e5a89c-4694-4a5a-906e-73b2adbd1336
# ╠═c5d6352a-fd03-4045-a91a-8574c15d9989
# ╠═3a4a992b-3050-43c0-b1f2-c7965a48d202
# ╠═92281858-95ad-4748-8088-cc6c74ae16dc
# ╠═8ef52f08-6518-4c82-912e-b3a20f4c81a7
# ╠═c289432f-5827-4651-8e2c-f6bd31d8fe24
# ╠═4dd7879c-147e-41d6-8527-83dbac04567b
# ╠═5797e1bf-7bf8-4343-b740-6f09c433f6c5
# ╠═e3ea8615-1bd2-40f3-bb79-c27e3c5375fe
# ╠═c4bfa3e8-142b-47a9-9aa4-0457cd15f6c4
# ╠═45c6e6ca-d9d5-4895-829f-3522bb4485e0
# ╠═474bfcdf-7c28-472a-acb0-07aa7fcc719d
# ╠═efe841ab-9d93-4792-917e-6c27916381af
# ╠═14a03646-a5b6-4381-a3ed-39a6eb85c113
# ╠═ccd02e1b-b33c-4aec-a95f-da12dbb4cda7
# ╠═68b614f2-dc13-4c46-9649-273179fbe27c
# ╠═4261573f-2207-45e2-b117-49c64079263a
# ╠═9c58b595-b903-4bfc-bf4d-6e7cf8acf588
# ╠═fb78a338-c2ab-499e-bddd-bf36a15ea8cc
# ╠═56b71cfd-b90b-4d2c-9307-2b3b97a30c6f
# ╠═c3bd2406-abc7-4e48-b8af-a231c9cda891
# ╠═50f0d299-5b25-452d-9a95-adb121d6abd4
# ╠═4ffa67f3-c4bd-47e6-885a-28cab178dd61
# ╠═ca232f87-1c5e-4cdf-8ed8-aef9e5a34db2
# ╠═fa0a86ad-a9e7-4d53-b7dd-7866f85c000f
# ╠═167ec4cf-c863-4a85-b429-02571fad0f30
# ╠═8193ba21-d8d2-40e0-af7b-3028d26afafc
# ╠═98eef687-eff7-4872-8541-85baeb577b77
# ╠═6e64c5d0-c9cf-4c4b-a9df-1f28225b7eff
# ╠═e2a35442-6974-4607-bbc0-b12ca32b4594
# ╠═2a19189a-8da9-4e1a-a76b-6d4a09e0ad0f
# ╠═9020198b-e68e-4253-a6ed-31a8d76944a6
# ╠═b2707113-fb5b-4600-bd48-4be882759436
# ╠═6a9d745e-e8ee-499c-ba72-fc95570ce442
# ╠═10b601b5-b24e-44e9-8d1a-c14aeba3dcf7
# ╠═b88cc8dd-39d6-4879-98a2-cccaf6e000b7
# ╠═2ab21d7e-7416-4057-a53e-8ebde10295e6
# ╠═84ce8c14-0f67-46a6-9386-818c8003818f
# ╠═c22bfed6-38cf-454f-93b1-6516e10e3b96
# ╠═a37676f5-e278-4df8-aa88-4f3ec1e3c414
# ╠═63775bee-b22e-4f99-a16c-d74ceaa91e65
# ╠═45a74d26-5ad9-404d-af14-33ae4f863a3d
# ╠═7f7ab4bd-69d4-498e-b6f1-a01bb5857ce5
# ╠═3ee9dae1-bccc-41b3-9980-3649683dae3c
# ╠═cecb4231-6cdf-42dd-81fc-7fb304e4ffc0
# ╠═b10c3e4b-edbc-474e-8c79-f37d9726439f
# ╠═a57b759e-9c4b-49c4-9439-4b842131f325
# ╠═15186cf0-a932-47b6-9f4c-3378b123de04
# ╠═01382010-99aa-4393-b6e0-7eaefd2b2ace
# ╠═89b6d6ee-1b73-4e2b-aa06-eda5a40979f7
# ╠═ff5c50ef-19df-4123-9da2-2efd43229b6f
# ╠═f46e89ce-dfa7-4eb4-91d4-6f8e97ab50ba
# ╠═e662491c-5165-4d3f-926f-fac8549e94e0
# ╠═178b9ad6-7784-4984-9684-3241ba7217dd
# ╠═564bb277-404c-413f-9e5b-1d56df39040f
# ╠═5a64f67e-f9f0-4942-a3f9-7c7e468e9225
# ╠═5b04879a-5029-4439-b4df-c1449fe27834
# ╠═0d4cccc7-56bf-45c2-8529-16adf8eaecc7
# ╠═1aff4c3c-f4fa-4cda-9810-f91514dc81a1
# ╠═cbe4aeeb-ad2a-4607-b010-20b34569773d
# ╠═daee61c5-6a65-4688-8f1d-7ee6c5a8d884
# ╠═45c19245-4dbb-4047-b0ed-13b27231c667
# ╠═cf80712e-bf9a-4fde-a058-c0404011dca2
# ╠═b1c4718b-3829-4b39-a84b-2445d605f08d
# ╠═ad211fed-2077-4869-bdcb-23c9804db113
# ╠═f336bb29-2ac3-41e0-9930-623e8fea3af5
# ╠═6038ed31-c1c5-4a3a-8809-9e6a39fd0e5b
# ╠═69b2fd33-615a-4bf3-8936-8c9ee0651ad1
# ╠═c4d4fb36-fcd4-4d1c-b5f3-f208ad48fd3f
# ╠═c2dd681e-c318-424d-9038-f66ab9317c52
# ╠═3028c7b7-c633-47b9-bb5d-f15101375b4f
# ╠═1cd3cbae-96dd-4af2-ab6f-eb2565e6d6e1
# ╠═49158b6b-be82-4689-9a7b-4de7e707e08c
# ╠═b872074b-094c-4d44-8d90-a1cf4908e012
# ╠═5051d08b-4a09-4e59-ae81-2caf9e555be4
# ╠═54daafd9-2842-43f4-bde3-eabdd6403767
# ╠═2a00bd6b-a0ce-4b54-bf88-d6cc2b2204c8
# ╠═28f216c0-f5b3-4f22-a5eb-192423b980e9
# ╠═e5b2c108-25a9-4035-bc60-5c5ca3eb54f1
# ╠═3198a48f-36e2-452a-9acd-e478a15057a3
# ╠═a80b460a-efac-4ecc-ad30-3a00e200f31e
# ╠═bd7d392b-f98f-4e98-9a5d-155d1d6fd199
# ╠═67791160-75c3-476f-b7e0-6d264a9dc1a4
# ╠═9e921706-7eda-4574-b8c5-ad368490b4e1
# ╠═fdec6c16-d095-4914-86bb-b94983b04465
# ╠═e13c04ce-bb90-480f-978c-2f02f01d31de
# ╠═559f597f-0d11-4212-b1e7-9ccd22661ae0
# ╠═3f04f5b8-1670-4155-865e-3212b3973364
# ╠═a81bfab4-dc7a-4fff-a72f-99a8dfd73a5d
# ╠═35afbb71-6862-4c0d-9582-2e80079f39c8
# ╠═72160b43-5668-49f8-81cc-076f3ab5b300
# ╠═0cc7f382-5424-4205-8b4f-5e1fc5e3bae3
# ╠═718fc4a4-8f46-4dff-a9a6-6584c134308a
# ╠═4894c24f-9729-4ebc-80a9-72f6e0763a4b
# ╠═01d54fa2-cd95-4150-af88-f1fbf8f1ae08
# ╠═990755c3-7cd5-49f9-9ddd-b86d6b267cb5
# ╠═84d06bf5-a160-4dcc-af3b-494b105d7e22
# ╠═501e8567-f7c2-47e6-adaf-d3a463012912
# ╠═ed32a5b0-a5ee-4eec-9093-657ce7f907a0
# ╠═7ac8a793-352b-4b73-aa3d-d93b0c7fce9e
# ╠═cf9aa152-1095-4d38-8641-942f0de4b860
# ╠═d6963ed3-39c3-4c60-9cfa-0d2f1fdd3a99
# ╠═d33d71f3-e16c-4001-b561-239e00aaa102
# ╠═7c982555-a277-49bf-add5-1b27ec066d10
# ╠═e4441a0e-55ad-47bd-8b30-004417d9dc71
# ╠═28a9dbdb-909a-43f3-a7bd-ace577293043
# ╠═2606492a-89fd-4c2a-8708-131784b0a053
# ╠═db3eedbf-ad89-487f-a8f9-f018f52cd2f0
# ╠═61b50198-e582-4171-8752-7108f5021d30
# ╠═1f56f505-ba6f-4f59-9985-ca844874e336
# ╠═04587beb-c645-49cb-b74e-bfa51e434fdc
# ╠═496eba70-c366-4ecf-99f1-850a0bb64ecf
# ╠═b829bbdd-84ef-4c66-af88-fd928518a73f
# ╠═08fd7243-13d8-47af-b4d1-2aa410927e93
# ╠═6d29b8d8-bb2e-4fc6-ae72-9f97265f1451
# ╠═6e06a810-186a-4d91-b63a-259c4ef32300
# ╠═dffe90b9-3b9a-4f9c-a7c7-365d894777ea
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
