### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ cff82de0-3d15-11ee-37ff-035acc6e8e4a
begin
	using DataFrames
	using CSV
	using Plots # default
	using PlotlyBase
	using PlutoUI
	using Dates
	using Optim
	using Peaks
	using Metrics
	using ProgressMeter
	using BenchmarkTools
	using CurveFit
	using ProgressLogging
	using Plots.PlotMeasures
end

# ╔═╡ 9ce6e2c1-69e4-47c6-830c-3c0eb13d784e
begin
	include("..//src//energy_draw.jl")
	include("..//src//time.jl")
	include("..//src//solar_radiation.jl")
	include("..//src//track.jl")
	include("..//src//utils.jl")
	include("..//src//strategy_calculation.jl")
end

# ╔═╡ 7c89a14e-f6c7-4a09-8825-a81d1ba9b824
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 1500px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 39d80a54-3675-47ff-82bc-deb1daf03a59
TableOfContents()

# ╔═╡ 1b605aae-c293-49b0-a06d-fcff2ef3c12a
md"# Loading data"

# ╔═╡ d6fb7d82-6bd9-4da9-b183-67fb16ab9f66
track,segments = get_track_and_segments("..//data//data_australia.csv");

# ╔═╡ aef61af9-31df-4757-95bd-7fdfc0edbc1a
track

# ╔═╡ dbe59867-4be2-4a6a-9236-35dc96d2d387
plot(track.distance/1e3, track.altitude, title="Original track data", legend=false, size=(300,250), xlabel="Distance (km)", ylabel="Altitude (m)")

# ╔═╡ 3327e2e2-4398-4aed-83c1-eaeb7565bd48
plot(
	track.distance[1:100]/1e3,
	track.altitude[1:100],
	title="Original track data",
	legend=false, size=(300,250),
	xlabel="Distance (km)",
	ylabel="Altitude (m)",
	markershape=:diamond
)

# ╔═╡ e53d252c-31b5-430c-a241-b41159243bc7
track_peaks, segments_peaks, points_peaks = keep_extremum_only_peaks_segments_with_points(track);

# ╔═╡ 164055ba-18bb-482d-816e-e278fc76a0d0
plot(track_peaks.distance, track_peaks.altitude, title="Peaks track data")

# ╔═╡ 1b6f993d-e4cc-4973-b5cf-269f01c7f409
md"# Parametrized merging"

# ╔═╡ 16405fa9-f44a-4b60-9b87-84bfd8aa8ed2
thresholds = 0:0.05:3.

# ╔═╡ 7513698b-e1cf-4dc1-a461-ede112a76180
collect(thresholds)

# ╔═╡ 57d53d47-733a-4baf-b2c4-f71d97429966
tracks = Dict()

# ╔═╡ c2da7fb6-e6f0-494c-9f6b-6c20e7127b52
segments_dict = Dict()

# ╔═╡ 5f90dc5b-1146-4085-8ff3-4dbe33854aad
points = Dict()

# ╔═╡ b5c8451c-03ba-42d2-98aa-be25daa2729b
begin
	for thr in thresholds
		tracks[thr], points[thr] = parametrized_track_simplification(track, thr);
		segments_dict[thr] = get_segments_for_track(tracks[thr]);
	end
end

# ╔═╡ 24f87654-2894-42a8-b13a-3d2a8eacc53b
tracks[0.05]

# ╔═╡ 51ae1ae6-98cb-4d76-9ce2-a0a95d88478d
plot(tracks[0.05].distance, tracks[0.05].altitude, title="0.05 threshold")

# ╔═╡ e9130295-ff5e-4cd5-92f0-844a4dfdcd8e
tracks[0.5]

# ╔═╡ 418a5c92-3a0e-4640-93db-24fc0dd98c51
plot(tracks[0.5].distance, tracks[0.5].altitude, title="0.5 threshold")

# ╔═╡ 7f4daafa-8a3f-46a2-81f4-fd3d6d3ade37
md"# Optimizing"

# ╔═╡ 1f039515-b3d7-49ae-aef5-f5951934778d
start_energy = 5100.

# ╔═╡ 5008c2a6-b61b-4150-a8e2-225fcbad8680
start_datetime =  DateTime(2023,1,1,10,0,0)

# ╔═╡ 0c7cbc87-96b7-4855-91e9-860996b2a5ba
speeds = minimize_single_speed(track, segments, start_energy,start_datetime, 50.)

# ╔═╡ 36d52b29-260c-427d-86f2-5d45d792e10a
opt_speed = speeds[1]

# ╔═╡ eb2f6b3b-6bdb-48b6-8c49-eb7f3a579b57
simulate_run(
	fill(opt_speed, size(segments,1)),
	track,
	segments,
	start_energy,
	start_datetime
)

# ╔═╡ b98aed7d-53f3-4ee4-91d6-9e5b3a9fc67b
md"Это был прогон на одной скорости на исходной трассе (полные данные)

Что надо сделать теперь:
Для каждого из упрощений трассы с параметром, прогнать симуляцию и собрать энергии/времена. Запомнить их
Взять исходную трассу и её энергии/времена. Взять только те точки, которые есть в упрощённой трассе. Показать разницу энергий/времён
В конце добавить ещё и peaks трассу"

# ╔═╡ 459d3b2c-0954-4bd8-9369-ac294e1adcb7
md"# Расчёт энергии для трассы"

# ╔═╡ 04a05ccd-b5d5-4ae9-934f-b9d2290a2d58
power_use_track, solar_power_track, time_s_track = solar_trip_boundaries(
		fill(opt_speed, size(segments,1)) / 3.6, segments, start_datetime
);

# ╔═╡ 316d420f-b698-4d3c-988b-732e55e5089a
energy_in_system_track = start_energy .+ solar_power_track .- power_use_track

# ╔═╡ 77aec1eb-7e4c-4c0e-a31f-f8f7036792df
pushfirst!(energy_in_system_track, start_energy);

# ╔═╡ 277fcf00-4e82-4ddc-a603-a53bd6695cf7
function simulate_run_energies(speed, track, segments, start_energy, start_datetime)
	speeds = fill(speed / 3.6, size(segments, 1));
	power_use, solar_power, time_s = solar_trip_boundaries(
		speeds, segments, start_datetime
	)
	energy_in_system = start_energy .+ solar_power .- power_use
	pushfirst!(energy_in_system, start_energy)
	return energy_in_system, time_s
end

# ╔═╡ 2f02f496-5548-46d5-9ba3-bca71b0b22b6
function simulate_run_income_use_time(speed, segments, start_datetime)
	speeds = fill(speed / 3.6, size(segments, 1));
	power_use, solar_power, time_s = solar_trip_boundaries(
		speeds, segments, start_datetime
	)
	return solar_power, power_use, time_s
end

# ╔═╡ 3d672dd6-82dd-402f-847c-67795f87ac8a
function calculate_energy_in_system(start_energy, solar_power, power_use)
	energy_in_system = start_energy .+ solar_power .- power_use
	pushfirst!(energy_in_system, start_energy)
	return energy_in_system
end

# ╔═╡ 6081d5ae-583d-4572-8930-facd41137bb4
function segment_data_to_track_data(values, first_value = 0.)
	new_value = [];
	push!(new_value, first_value);
	append!(new_value, values);
	return new_value;
end

# ╔═╡ b9b6d1ec-208c-4684-a6d4-3e0c89a0277d
md"sanity check. energy difference should be 0 for track and track"

# ╔═╡ a23c2215-531a-47f0-8241-6e0ac06f6d46
energy_track, time_track = simulate_run_energies(opt_speed, track, segments, start_energy, start_datetime)

# ╔═╡ 4d48b247-499b-42f6-8e24-085aee654de4
plot(energy_track-energy_in_system_track)

# ╔═╡ 611dec1f-6d0e-4b6d-9f1e-a547aec312ac
md"# Эксперименты с сокращенными трассами"

# ╔═╡ 14e1a205-8cd3-4444-957c-6902e7c41d6a
tracks[0.05]

# ╔═╡ 466171d6-a749-4d18-af72-099e099dd5bd
size(tracks[0.05], 1)

# ╔═╡ 0131849a-9f90-4097-ab3d-9bd47449bdc2
size(segments_dict[0.05], 1)

# ╔═╡ f7717337-d999-4150-b888-fb0644fd00d1
points[0.05]

# ╔═╡ 5ce6ed99-e0c1-4100-a5a8-b92e2bf10fc0
track[points[0.05],:]

# ╔═╡ 700120b8-78dd-4f11-8b91-e0bbd267be60
size(track[points[0.05],:], 1)

# ╔═╡ 64f4a5bc-7d13-45eb-a2c8-3325166cc564
function compare_track_energies(full_track, full_segments, reduced_track, reduced_segments, reduced_points, speed, start_energy, start_datetime)
	full_energy, full_time = simulate_run_energies(speed, full_track, full_segments, start_energy, start_datetime)
	reduced_energy, reduced_time = simulate_run_energies(speed, reduced_track, reduced_segments, start_energy, start_datetime)
	energy_diff = full_energy[reduced_points] - reduced_energy;
	# time_diff = full_time[reduced_points] - reduced_time;
	# return energy_diff, time_diff
	return energy_diff
end

# ╔═╡ 8f4fc0ed-77e4-4d30-a5f4-d2d6af16bc56
function compare_track_energies_income_use(full_track, full_segments, reduced_track, reduced_segments, reduced_points, speed, start_energy, start_datetime)
	full_income, full_use, full_time = simulate_run_income_use_time(speed, full_segments, start_datetime);
	reduced_income, reduced_use, reduced_time = simulate_run_income_use_time(speed, reduced_segments, start_datetime);

	full_energy = calculate_energy_in_system(start_energy, full_income, full_use);
	reduced_energy = calculate_energy_in_system(start_energy, reduced_income, reduced_use);

	full_income_track = segment_data_to_track_data(full_income, 0.);
	reduced_income_track = segment_data_to_track_data(reduced_income, 0.);
	full_use_track = segment_data_to_track_data(full_use, 0.);
	reduced_use_track = segment_data_to_track_data(reduced_use, 0.);
	full_time_track = segment_data_to_track_data(full_time, 0.);
	reduced_time_track = segment_data_to_track_data(reduced_time, 0.);

	# everything here is in segments. check segments legths
	income_diff = full_income_track[reduced_points] - reduced_income_track;
	use_diff = full_use_track[reduced_points] - reduced_use_track;
	time_diff = full_time_track[reduced_points] - reduced_time_track;
	energy_diff = full_energy[reduced_points] - reduced_energy;
	return income_diff, use_diff, time_diff, energy_diff
end

# ╔═╡ 9e2eb0bc-0ffa-4e10-be7a-df9aebc6b609
test_en = compare_track_energies(
	track, segments,
	tracks[0.05], segments_dict[0.05], points[0.05],
	opt_speed, start_energy, start_datetime
)

# ╔═╡ 4e2f9485-9131-40e6-afa5-e6134a638750
track[points[0.05],:]

# ╔═╡ aa894177-a0a5-4488-a784-3752edaeddcb
tracks[0.05]

# ╔═╡ 97690c57-f630-4762-bdb8-6630cc14c3a7
test_income, test_use, test_time, test_energy = compare_track_energies_income_use(
	track, segments,
	tracks[0.05], segments_dict[0.05], points[0.05],
	opt_speed, start_energy, start_datetime
)

# ╔═╡ f61f4da9-0689-4752-ae25-27415a95c10d
function plot_differences(income,use,time,energy, track, points)
	track_points = track[points, :];
	income_plot = plot(track_points.distance, income, title="Разница в получении энергии (Вт*ч)");
	use_plot = plot(track_points.distance, use, title="Разница в тратах энергии (Вт*ч)");
	time_plot = plot(track_points.distance, time, title="Разница во времени (с)", xlabel="Дистанция (м)");
	energy_plot = plot(track_points.distance, energy, title="Итоговая разница в энергии (Вт*ч)");

	plot(income_plot, use_plot, energy_plot, time_plot, layout=(4,1), size=(1000,700), legend=false)
end

# ╔═╡ dadc2e8c-3aae-4944-b5e6-cc99b9c4e1f5
function plot_differences_en(income,use,time,energy, track, points)
	track_points = track[points, :];
	income_plot = plot(track_points.distance, income, title="Energy income difference");
	use_plot = plot(track_points.distance, use, title="Energy drain difference");
	time_plot = plot(track_points.distance, time, title="Time difference (seconds)", xlabel="Distance (meters)");
	energy_plot = plot(track_points.distance, energy, title="Energy difference (Wt*h)");

	plot(
		# income_plot, use_plot,
		energy_plot, time_plot, layout=(2,1), size=(600,400), legend=false)
end

# ╔═╡ 0aeeb368-896e-4bc8-ad1c-d92b0c324392
plot_differences(test_income, test_use, test_time, test_energy, track, points[0.05])

# ╔═╡ cdf8794f-67ed-41df-bd9e-b1c57c82133b
plot(track[points[0.05],:].distance, [test_income, test_use, test_time, test_energy], title="Difference", labels=["income" "use" "time" "energy"])

# ╔═╡ eb2ddfbc-a8ab-4b88-ac46-37af7c563d42
function compare_track_energies_thr(full_track, full_segments, thr, speed, start_energy, start_datetime)
	return compare_track_energies(
	full_track, full_segments,
	tracks[thr], segments_dict[thr], points[thr],
	speed, start_energy, start_datetime
)
end

# ╔═╡ 501fc5bb-8f50-4a19-b2a0-3fafe9a7738e
function compare_track_energies_thr(thr, speed)
	return compare_track_energies(
	track, segments,
	tracks[thr], segments_dict[thr], points[thr],
	speed, start_energy, start_datetime
)
end

# ╔═╡ 322e9eb5-7da5-492f-a44d-2f7ee8c2ffcb
energy_diff_0_05 = compare_track_energies_thr(track, segments, 0.05, opt_speed, start_energy, start_datetime)

# ╔═╡ fcbeb07f-3c3f-4973-ac1a-f133c7bdd4f3
plot(energy_diff_0_05)

# ╔═╡ 3cd8689b-5f6c-481b-9e5b-ec0356597598
plot(
	compare_track_energies_thr(
		0.05,
		opt_speed
	)
)

# ╔═╡ 35102b6b-1a3f-4c99-b03a-86fb5ad00e14
plot(
	compare_track_energies_thr(
		0.2,
		60.
	)
)

# ╔═╡ 6950f324-b23c-4b49-b391-586bda1b872e
function compare_track_energies_plot_thr(track, segments, speed, thr)

	test_income, test_use, test_time, test_energy = compare_track_energies_income_use(
		track, segments,
		tracks[thr], segments_dict[thr], points[thr],
		opt_speed, start_energy, start_datetime
	)

	return plot_differences(test_income, test_use, test_time, test_energy, track, points[thr])
end

# ╔═╡ a8f3f8a9-673c-497b-a311-89798970f358
function compare_track_energies_plot_thr_en(track, segments, speed, thr)

	test_income, test_use, test_time, test_energy = compare_track_energies_income_use(
		track, segments,
		tracks[thr], segments_dict[thr], points[thr],
		opt_speed, start_energy, start_datetime
	)

	return plot_differences_en(test_income, test_use, test_time, test_energy, track, points[thr])
end

# ╔═╡ ce77f956-bd0a-48c1-a2b2-1f26a49877d8
compare_track_energies_plot_thr(track, segments, opt_speed, 0.)

# ╔═╡ 425a892c-b550-486e-b376-8566ae086f6a
compare_track_energies_plot_thr(track, segments, opt_speed, 0.05)

# ╔═╡ d9a12d12-8c19-4cd3-88ba-847e1c273850
compare_track_energies_plot_thr(track, segments, opt_speed, 0.25)

# ╔═╡ f84953c2-bed0-4375-a9d1-0bb646c832e9
compare_track_energies_plot_thr(track, segments, opt_speed, 0.5)

# ╔═╡ c26613ea-827d-46b8-9bd2-15ddf6cd885b
compare_track_energies_plot_thr(track, segments, opt_speed, 0.75)

# ╔═╡ b845f4f8-22c0-49cb-9580-d51c19e115bd
compare_track_energies_plot_thr(track, segments, opt_speed, 1.)

# ╔═╡ d95e9d51-867e-484a-8006-e1fb2d1e0fbf
compare_track_energies_plot_thr_en(track, segments, opt_speed, 1.75)

# ╔═╡ baa774fc-d955-4b69-8a32-7bf7b6a92997
md"# Сравнение с peaks"

# ╔═╡ da13556d-62b0-49ce-867f-f64b8f227599
p_income, p_use, p_time, p_energy = compare_track_energies_income_use(
	track, segments,
	track_peaks, segments_peaks, points_peaks,
	opt_speed, start_energy, start_datetime
)

# ╔═╡ ecc4c5e7-7b9f-47b2-9f35-85f6c241f4ea
plot_differences(
	p_income, p_use, p_time, p_energy, track, points_peaks
)

# ╔═╡ 691e1f9a-5d10-475f-9ede-8d687f9d07b1
md"# Новое получение peaks"

# ╔═╡ f17bfc9f-01b5-4535-b05e-a4251e6a3594
peak_points = get_peak_points_plateau(track.altitude)

# ╔═╡ 1c54eea3-b89a-4c7c-a7ab-7759f1bd4242
track_peaks_pl, segments_peaks_pl = get_track_and_segments_for_selected_points(track, peak_points)

# ╔═╡ 89fdba6b-0c3f-4dd2-ae5b-a26dc1a90368
pl_income, pl_use, pl_time, pl_energy = compare_track_energies_income_use(
	track, segments,
	track_peaks_pl, segments_peaks_pl, peak_points,
	opt_speed, start_energy, start_datetime
)

# ╔═╡ 0f248c56-5ab5-4183-9260-634be0fc9736
plot_differences(
	pl_income, pl_use, pl_time, pl_energy, track, peak_points
)

# ╔═╡ 3e84943e-1f80-44d6-b3fc-9f7c0553c4c1
plot(track.distance[peak_points,:], pl_use)

# ╔═╡ 1bb03b35-2d86-4fc4-a7ed-d2bcf7120d7d
md"# Более точная высота в peaks"

# ╔═╡ 67ac9252-66d6-465b-9857-abd0611fb81c
md"В чём суть: сейчас высота сегмента считается как среднее между ограничивающими его точками.

То есть банально (y1+y2)/2. Но это работает только для одного участка. Для двух и более это работать не будет, потому как не ограничивающие точки не учитываются, что неправильно.


Из-за этого получение энергии будет считаться немного неправильно."

# ╔═╡ 950d2611-9730-4552-af0e-72e9cfc801d3
md"Например, посмотрим на следующий график:

Крайние точки одни и те же, но средние на нескольких участках будут разные. Поэтому и результат должен отличаться"

# ╔═╡ 5f5d1144-b83b-4c67-a20e-669387b93493
begin
	x_dots = [1,2,4]
	y_low = [7,4,1]
	y_high = [7,6,1]
	y_mid = [7,5,1]
end

# ╔═╡ a472bc02-13d7-4ce1-b74a-2a783f0a7965
plot(
	x_dots,
	[y_low y_high y_mid],
	linestyle=[:solid :solid :dot],
	label=["Ниже прямой, среднее=3.5" "Выше прямой, среднее=4.5" "Прямая между ТИ, среднее=4.0"],
	xlabel="Дистанция",
	ylabel="Высота",
	# markershapes=:diamond
	markershapes=[:diamond :square :circle],
	size=(500, 300)
)

# ╔═╡ ae90c511-3fcd-405f-b6fb-7814da8893aa
plot(
	x_dots,
	[y_low y_high y_mid],
	linestyle=[:solid :solid :dot],
	label=["Lower than straight, avg=3.5" "Higher than straight, avg=4.5" "Straight line between POI, avg=4.0"],
	markershapes=[:diamond :square :circle],
	xlabel="Distance",
	ylabel="Altitude",
	size=(500, 250),
	legend=:bottomleft
)

# ╔═╡ 171cc89f-05a6-49df-bd9a-deab45341778
md"Что надо делать - считать среднее значение функции на нескольких участках.

Это можно сделать, посчитав площадь под графиком (проинтегрировав) и поделив на длину участка.

То есть вместо get mean data сделать последовательное численное интегрирование. Методом средних прямоугольников, так по идее будет нулевая ошибка"

# ╔═╡ 1a4155e3-a791-4f5c-bc46-687dd17c3742
function get_average_on_segment(data_x, data_y, from, to)
	integrated_value = 0
	# for i=1:length(data_y)-1
	for i=from:to-1
		integrated_value += (data_y[i] + data_y[i+1]) / 2 * (data_x[i+1] - data_x[i])
	end
	result = integrated_value / (data_x[to] - data_x[from])
	return result
end

# ╔═╡ aa9edc8e-6357-4c3e-a8c2-f38241de7f38
md"Проверим на данных с графика"

# ╔═╡ fbdbf4dc-b1b4-4eef-836f-558c276b83d1
get_average_on_segment(x_dots, y_mid, 1, 3)

# ╔═╡ b031cc68-9306-4265-973b-e427b801f439
get_average_on_segment(x_dots, y_low, 1, 3)

# ╔═╡ e6cbc2b6-7c6a-4405-9899-4af0f7ba979c
get_average_on_segment(x_dots, y_high, 1, 3)

# ╔═╡ 170467ab-6e26-4e24-9e2a-c244bdf08b9d
md"Похоже на правду"

# ╔═╡ 17a4b5b0-3aa7-4f60-b843-df81823fa710
get_average_on_segment(x_dots, [4,4,1], 1, 3)

# ╔═╡ 4110b150-d44b-4e2d-9769-2edac7c1a81c
get_average_on_segment(x_dots, [4,2,1], 1, 3)

# ╔═╡ 4036eb31-d5a1-404c-aa5d-ecb26badaa70
get_average_on_segment(x_dots, [4,3,1], 1, 3)

# ╔═╡ d7d2142b-b2e6-44e8-aaa9-c293b2a9b70c
get_average_on_segment(x_dots, [4,2,1], 1, 1)

# ╔═╡ 58436584-cddc-44ce-b19e-10a1fc7873cf
get_average_on_segment(x_dots, [4,3,1], 1, 2)

# ╔═╡ ef6dfaa8-987f-4f79-9957-16db862aa1ff
get_average_on_segment(x_dots, [4,3,1], 2, 3)

# ╔═╡ 0aaabc2f-aac7-4f03-90e3-a42e95f08e62
md"Да, всё сходится. Теперь пора сделать это для формирования сегментов

Опционально, сделать код векторизируемым. (Но зачем? и так будет быстро работать)"

# ╔═╡ 84f39121-69bb-4a04-8ca4-b89af0f6a0f8
md"Но остаётся вопрос: а что с широтой, долготой и углом наклона (slope)?

По идее коррдинаты надо таким же образом считать.

С углом наклона непонятно, сейчас это разница высоты / разница дистанции"

# ╔═╡ 6c33a6a3-808a-4086-914d-6a23d749ee38
md"По идее так же надо считать с учётом профиля. Только здесь не среднее значение функции на интервале, а что-то другое.

Или не надо. Стоит подумать"

# ╔═╡ 6ef045bd-2426-41f6-9f61-4a17f136ac7a
function get_track_and_segments_for_selected_points_modified(track, points)
	# constructing proper lat, lon and alt values with numerical integration
	new_altitude = [];
	new_longitude = [];
	new_latitude = [];
	for i=1:length(points)-1
		push!(new_altitude, get_average_on_segment(track.distance, track.altitude, points[i], points[i+1]))
		push!(new_longitude,  get_average_on_segment(track.distance, track.longitude, points[i], points[i+1]))
		push!(new_latitude,  get_average_on_segment(track.distance, track.latitude, points[i], points[i+1]))
	end

	new_track = copy(track)
	new_track.index = 1:size(new_track, 1)
	new_track = new_track[points,:]
	
	new_segments = DataFrame(
        from = new_track.index[1:size(new_track, 1) - 1],
        to = new_track.index[2:size(new_track, 1)],
        diff_distance = diff(new_track.distance),
        diff_altitude = diff(new_track.altitude)
    )

	new_segments.slope = atand.(new_segments.diff_altitude ./ new_segments.diff_distance)

	# new_segments.latitude = get_mean_data(new_track.latitude)
    # new_segments.longitude = get_mean_data(new_track.longitude)
	new_segments.latitude = new_latitude
	new_segments.longitude = new_longitude
	# new_segments.altitude = get_mean_data(new_track.altitude)
	new_segments.altitude = new_altitude
	new_segments.weather_coeff .= 1.0
	
	return new_track, new_segments
end

# ╔═╡ 4642901d-e71b-41bc-8d34-5e1f01f83306
track_new_peaks, segments_new_peaks = get_track_and_segments_for_selected_points_modified(track, peak_points)

# ╔═╡ 6e697d1a-5823-4f9f-a3e0-c649adc0b449
new_peaks_income, new_peaks_use, new_peaks_time, new_peaks_energy = compare_track_energies_income_use(
	track, segments,
	track_new_peaks, segments_new_peaks, peak_points,
	opt_speed, start_energy, start_datetime
)

# ╔═╡ 4d6b8b1c-081b-4db2-bede-6d67cb41fbaf
plot_differences(
	new_peaks_income, new_peaks_use, new_peaks_time, new_peaks_energy, track, peak_points
)

# ╔═╡ 28275f01-96c4-4988-8d2c-19648ea388b6
md"Надо сравнить с предыдущей версией peaks"

# ╔═╡ 9b5b321a-ed7f-4b88-8519-498f2d19fb1c
segments_new_peaks.altitude - segments_peaks_pl.altitude

# ╔═╡ 09dddfde-100f-43d9-ad48-65eecb39923b
md"Как ни странно, получается хуже

Не хуже, а по-другому. Выходят другие высоты"

# ╔═╡ 91cae49c-bd36-473f-a23b-d15ecd1ac073
sum(segments_new_peaks.altitude - segments_peaks_pl.altitude)

# ╔═╡ dea20b86-88d4-4337-8f58-25606c7adff1
md"Разница в минус по высоте, поэтому и меньше энергии получается"

# ╔═╡ 9d390be8-b231-4fc8-827a-f561293e82b9
md"Ещё надо посмотреть насколько резкие перепады в power use и почему"

# ╔═╡ 1550e304-31d0-4d4c-b5a2-8e4d72561f05
size(track_new_peaks, 1)

# ╔═╡ 55f67256-7b17-40a3-a739-d502ac5864fb
size(track, 1) / size(track_new_peaks, 1)

# ╔═╡ 1e7a691c-d3cd-4d4c-a202-311e222a6320
md"## А теперь для k сравним"

# ╔═╡ 5da7a24e-b075-4d6f-ab20-e7bfd4a2aed9
track_new_1_75, segments_new_1_75 = get_track_and_segments_for_selected_points_modified(track, points[1.75])

# ╔═╡ b62e5863-7bf9-4e96-9891-267dc4baa737
reg_1_75_income, reg_1_75_use, reg_1_75_time, reg_1_75_energy = compare_track_energies_income_use(
	track, segments,
	tracks[1.75], segments_dict[1.75], points[1.75],
	opt_speed, start_energy, start_datetime
)

# ╔═╡ 9ab0a9a7-c071-49ba-bc0d-077b44c5c041
new_1_75_income, new_1_75_use, new_1_75_time, new_1_75_energy = compare_track_energies_income_use(
	track, segments,
	track_new_1_75, segments_new_1_75, points[1.75],
	opt_speed, start_energy, start_datetime
)

# ╔═╡ 26c9ff7f-2ff7-4fe8-9abe-d72f4d869c5c
md"Сперва без численного интегрирования"

# ╔═╡ d6645aef-70a6-4af3-8cb4-4254d80a7281
plot_differences(
	reg_1_75_income, reg_1_75_use, reg_1_75_time, reg_1_75_energy, track, points[1.75]
)

# ╔═╡ 95fd37d2-f7b0-479b-9779-093efcc11066
md"Потом с численным интегрированием"

# ╔═╡ e42e57dd-5513-4ca8-9ddf-5816f94f5707
plot_differences(
	new_1_75_income, new_1_75_use, new_1_75_time, new_1_75_energy, track, points[1.75]
)

# ╔═╡ 205cc77c-2f01-4552-a660-08bd750d9b43
size(track,1)/size(points[1.75],1)

# ╔═╡ aa4596d6-a281-4ed3-8115-f5f3e30d1f35
md"Разницы нет?" 

# ╔═╡ 9af6b676-d183-4f52-b92a-8017a0df4c85
last(new_1_75_energy)

# ╔═╡ aa4c6f7a-d6ca-430e-8297-5e07e9e9e2b3
last(reg_1_75_energy)

# ╔═╡ 81eac0f6-9c2f-4123-886f-c0317d70444b
md"Похоже что немного есть таки в пользу правильного усреднения"

# ╔═╡ 2790174f-7a2e-46c9-967e-065c4f9a12fd
last(new_1_75_use)

# ╔═╡ d36a4805-e59b-45ba-9433-54fd26733670
last(reg_1_75_use)

# ╔═╡ 8695b7a4-8127-4766-a7f5-fcc7d7f5b331
md"Одинаково по использованию энергии (логично, среднее не используется)"

# ╔═╡ a9bad1e6-7a9c-49c0-a0aa-33aaa6776691
last(new_1_75_income)

# ╔═╡ 0059a826-9daa-4d9c-9213-b1acefb8b92f
last(reg_1_75_income)

# ╔═╡ cbc854d5-d287-495c-9921-408ce23d1f30
md"А вот без усреднения хуже приход энергии" 

# ╔═╡ af6550be-7453-41ef-bcbe-dd41312eec76
md"Итого, integrated k=1.75 all the way!" 

# ╔═╡ 5b770bfc-64d1-4b8e-b30d-e40cf4eb5397
md"# Анализ power use" 

# ╔═╡ f9dc21cb-d389-454d-99da-383208e268f9
md"Или вообще не рассматривать эту проблему, т.к. она небольшая?" 

# ╔═╡ 76438295-379c-4d08-a57c-ac3b341f3a90
md"# Сравнение параметрического построения циклом"

# ╔═╡ 51bd9a66-5e73-4fe4-b7cf-47ed699d8011
md"## Сперва без интегрированной altitude"

# ╔═╡ cd865e2f-c24f-4b6d-a8ed-45de9867e7b9
begin
	energies = Dict()
	p=plot()
	for thr in 0:0.1:1#thresholds
		energies[thr] = compare_track_energies(
			track, segments,
			tracks[thr], segments_dict[thr], points[thr],
			opt_speed, start_energy, start_datetime
		)
		plot!(track.distance[points[thr]], energies[thr], label="Energy $thr")
	end
	plot(p)
end

# ╔═╡ ffbb006a-47be-404a-a564-d171e70c0590
function draw_comparison_for_parametrized_peaks(thresholds, track, segments)
	energies = Dict()
	p=plot()
	for thr in thresholds
		track_thr, points_thr = parametrized_track_simplification(track, thr);
		segments_thr = get_segments_for_track(track_thr);
		energies_thr = compare_track_energies(
			track, segments,
			track_thr, segments_thr, points_thr,
			opt_speed, start_energy, start_datetime
		)
		plot!(track.distance[points_thr], energies_thr, label="Energy $thr")
	end
	plot(p)
end

# ╔═╡ 40e187fb-8f63-4b09-bb0e-110633a1ba27
draw_comparison_for_parametrized_peaks([0., .01, .02, .03, .04, .05, .1, .25, .5, .75, 1.], track, segments)

# ╔═╡ cfb3b295-61a3-4851-b511-6b9d6b162ded
draw_comparison_for_parametrized_peaks([0., .01, .02, .03, .04, .05], track, segments)

# ╔═╡ 256c6718-883a-479b-838a-2e939f438c3b
draw_comparison_for_parametrized_peaks([0., .01, .02, .041], track, segments)

# ╔═╡ 58cc06b0-14d1-476f-a5a4-ca56ef022695
md"## Peaks"

# ╔═╡ 697ff289-724d-4f3b-a89d-4ec3a96567bc
function draw_comparison_for_regular_peaks(track, segments, track_peaks, segments_peaks, points_peaks)
	energies_thr = compare_track_energies(
		track, segments,
		track_peaks, segments_peaks, points_peaks,
		opt_speed, start_energy, start_datetime
	)
	plot(track.distance[points_peaks], energies_thr, label="Peaks")
end

# ╔═╡ cde20eed-7581-43f1-9abc-39e6a0555478
draw_comparison_for_regular_peaks(
	track, segments,
	track_peaks, segments_peaks, points_peaks
)

# ╔═╡ 1133a183-3631-4ef8-b2b8-1b697a29b75f
draw_comparison_for_parametrized_peaks([0.1, 0.5, 0.7, 0.75, 1], track, segments)

# ╔═╡ 6f342cb3-c8e2-4358-86b9-b45ad05d2551
md"Похоже что peaks это что-то близкое к 0.7"

# ╔═╡ 151dfb5e-f19a-4f7e-8764-12f99185094b
md"## Сравнение качество количество"

# ╔═╡ eede70d8-37e2-4082-a8ae-cd1491e768fa
md"Сделать таблицу сравнения, где будет количество участков и насколько снижена точность. Для thr и для peaks

Качество по r_squared, mse?"

# ╔═╡ 6b27c66e-873b-4ad9-8b9b-8987ee09c03b
function make_comparison_df(track, segments, thresholds, speed, start_energy, start_datetime)
	# 1. reference result
	full_energy, full_time = simulate_run_energies(speed, track, segments, start_energy, start_datetime)
	res_df = DataFrame(Threshold=Float64[], Finish_diff=Float64[], MAE=Float64[], MSE=Float64[], RMSE=Float64[], R2=Float64[], Length=Int64[], ExecTime=Float64[])
	
	for thr in thresholds
		track_thr, points_thr = parametrized_track_simplification(track, thr);
		segments_thr = get_segments_for_track(track_thr);
		t_start = time()
		reduced_energy, reduced_time = simulate_run_energies(speed, track_thr, segments_thr, start_energy, start_datetime)
		t_finish = time()
		exec_time = t_finish - t_start
		energy_diff = full_energy[points_thr] - reduced_energy;
		source = full_energy[points_thr];
		new_energy = reduced_energy;
		# это вся разница
		# считаем от неё метрики: r^2, MSE, RMSE разница на финише
		# sqrt(sum((arr1 .- arr2).^2) / length(arr1)) 
		mae_val = mae(source, new_energy)
		mse_val = mse(source, new_energy)
		rmse_val = sqrt(mse_val)
		r2_val = r2_score(source, new_energy)
		last_diff = last(source) - last(new_energy)
		number_of_segments = length(points_thr)

		push!(res_df, (thr, last_diff, mae_val, mse_val, rmse_val, r2_val, number_of_segments, exec_time))
		
		# energies_thr = compare_track_energies(
		# 	track, segments,
		# 	track_thr, segments_thr, points_thr,
		# 	speed, start_energy, start_datetime
		# )
	end

	# посчитать в отдельной функции
	# # добавить результаты обычного peaks (не параметрического)
	# track_peaks, segments_peaks, points_peaks2 = keep_extremum_only_peaks_segments_with_points(track);
	# energies_peaks = compare_track_energies(
	# 	track, segments,
	# 	track_peaks, segments_peaks, points_peaks2,
	# 	speed, start_energy, start_datetime
	# )

	# # добавить результаты нового peaks (параметрического)
	# track_mod_peaks, segments_mod_peaks = get_track_and_segments_for_selected_points_modified(track, points_peaks2)
	# energies_mod_peaks = compare_track_energies(
	# 	track, segments,
	# 	track_mod_peaks, segments_mod_peaks, points_peaks2,
	# 	speed, start_energy, start_datetime
	# )
	
	return res_df
end

# ╔═╡ e593df22-ebbe-4096-894e-4414b4013d0b
function make_comparison_peaks(track, segments, speed, start_energy, start_datetime)
	# 1. reference result
	full_energy, full_time = simulate_run_energies(speed, track, segments, start_energy, start_datetime)

	res_df = DataFrame(Series=String[], Finish_diff=Float64[], MAE=Float64[], MSE=Float64[], RMSE=Float64[], R2=Float64[], Length=Int64[], ExecTime=Float64[])
	
	# добавить результаты обычного peaks (не параметрического)
	track_peaks, segments_peaks, points_peaks = keep_extremum_only_peaks_segments_with_points(track);
	t_start = time()
	reduced_energy, reduced_time = simulate_run_energies(speed, track_peaks, segments_peaks, start_energy, start_datetime)
	t_finish = time()
	exec_time = t_finish - t_start
	# energy_diff_peaks = full_energy[points_peaks] - reduced_energy;

	source = full_energy[points_peaks];
	new_energy = reduced_energy;
	mae_val = mae(source, new_energy)
	mse_val = mse(source, new_energy)
	rmse_val = sqrt(mse_val)
	r2_val = r2_score(source, new_energy)
	last_diff = last(source) - last(new_energy)
	number_of_segments = length(points_peaks)

	# push!(res_df, ("Regular peaks", last_diff, mae_val, mse_val, rmse_val, r2_val, number_of_segments))

	# добавить результаты нового peaks (параметрического)
	track_mod_peaks, segments_mod_peaks = get_track_and_segments_for_selected_points_modified(track, points_peaks)
	reduced_mod_energy, reduced_mod_time = simulate_run_energies(speed, track_mod_peaks, segments_mod_peaks, start_energy, start_datetime)
	# energy_diff_peaks = full_energy[reduced_points] - reduced_mod_energy;

	source = full_energy[points_peaks];
	new_energy = reduced_mod_energy;
	mae_val = mae(source, new_energy)
	mse_val = mse(source, new_energy)
	rmse_val = sqrt(mse_val)
	r2_val = r2_score(source, new_energy)
	last_diff = last(source) - last(new_energy)
	number_of_segments = length(points_peaks)

	push!(res_df, ("Экстремумы", last_diff, mae_val, mse_val, rmse_val, r2_val, number_of_segments, exec_time))
end

# ╔═╡ e46f1433-f940-4b3c-a6c5-49a70c07f671
md"Завершить создание сравнительной таблички)"

# ╔═╡ cc352146-3a5e-4431-919e-649995ea7fed
thr_df = make_comparison_df(track, segments, [.0:.01:.1; .15:.025:1.5], opt_speed, start_energy, start_datetime)

# ╔═╡ e8944e65-f0e0-4169-bcb1-3d6aa2b30bf5
plot(
	thr_df.Threshold,
	Matrix(thr_df)[:,2:end-2],
	labels=["Finish diff" "MAE" "MSE" "RMSE" "R^2"]
)

# ╔═╡ 6ee8254c-be59-46a2-8b9c-5e47b50237f7
md"Ура! ошибка очень мала

Посмотрим что будет, если чрезмерно упрощать"

# ╔═╡ 39456f2f-700d-4476-9c96-ad818c058363
thr_df_large = make_comparison_df(track, segments, .0:.25:5, opt_speed, start_energy, start_datetime)

# ╔═╡ 9a70e256-3fbe-4660-ae01-cd02d5136a89
plot(
	thr_df_large.Threshold,
	[thr_df_large.Finish_diff thr_df_large.MAE thr_df_large.RMSE thr_df_large.R2 thr_df_large.Length],
	labels=["Finish diff" "MAE" "RMSE" "R^2" "Number of segments"]
)

# ╔═╡ dd8cf5c2-ccb9-4e55-8b1a-7c5cde97c36a
plot(
	thr_df_large.Length,
	[thr_df_large.Finish_diff thr_df_large.MAE thr_df_large.RMSE thr_df_large.R2 ],
	labels=["Finish diff" "MAE" "RMSE" "R^2" ],
	xscale=:log10,
	xlabel="Length, segments"
)

# ╔═╡ b2c8b3cb-7f6f-4f1e-9ccd-605bd110756a
md"Больше 2 threshold брать точно не следует

По-хорошему вообще больше 1, исходя из разницы на финише"

# ╔═╡ fc619c3a-556b-4bf7-bff9-28f9a9bc0fb2
peaks_df = make_comparison_peaks(track, segments, opt_speed, start_energy, start_datetime)

# ╔═╡ 15bc9488-d29d-468e-8995-1ff5df9487e6
md"Свести оба фрейма вместе с thr в виде стринги"

# ╔═╡ 0be83fa6-81c6-4d32-a2ff-db833750a3f1
names_series = vcat(string.(thr_df_large.Threshold[1:3]), peaks_df.Series, string.(thr_df_large.Threshold[4:end]))

# ╔═╡ b7c1cde6-7d63-4497-8a3e-2df959ca6918
total_df = vcat(thr_df_large[1:3,2:end], peaks_df[:,2:end], thr_df_large[4:end,2:end])

# ╔═╡ 9de59850-3bd2-46a9-a3f9-3b73055bdf1d
total_df.name = names_series

# ╔═╡ 7b07a2b8-eb69-4da7-ada5-4d44503e4274
Plots.scatter(
	total_df.name,
	[
		total_df.Finish_diff total_df.MAE total_df.RMSE total_df.R2 total_df.Length
	],
	labels = [
		"Finish diff" "MAE" "RMSE" "R^2" "Number of segments"
	],
	xlabel="K"
)

# ╔═╡ 5d239b82-f341-446b-9aaf-64076385723f
begin
	Plots.scatter(
		total_df.name,
		total_df.Length,
		label="Length",
		xlabel="K",
		color=:red,
		ylabel="Количество участков",
		# legend = :outertopright
		legend=false
	)
	Plots.scatter!(
		twinx(),
		total_df.name,
		total_df.R2,
		label="R^2",
		# xlabel="K",
		ylabel="R^2",
		# legend = :outerbottomright
		legend=false
	)
end

# ╔═╡ 6a4bf9a3-fbb9-494a-b48e-000553cf5460
begin
	Plots.scatter(
		total_df.name,
		total_df.Length,
		xlabel="K",
		color=:red,
		ylabel="Количество участков",
		# legend = :outertopright
		legend=false
	)
	Plots.scatter!(
		twinx(),
		total_df.name,
		total_df.RMSE,
		# xlabel="K",
		ylabel="RMSE",
		# legend = :outerbottomright
		legend=false
	)
end

# ╔═╡ 4b17b7fd-3a5d-4560-bcaf-bffb1d096ab5
Plots.scatter(
	total_df.name,
	[
		abs.(total_df.Finish_diff) total_df.Length
	],
	labels = [
		"Ошибка на финише (Вт*ч)" "Количество участков"
	],
	xlabel="K"
)

# ╔═╡ aac0c842-14b6-48e3-bde1-30e32713046f
Plots.scatter(
	total_df.name,
	[
		total_df.Finish_diff./start_energy total_df.Length
	],
	labels = [
		"Ошибка на финише (нормировано)" "Количество участков"
	],
	xlabel="K"
)

# ╔═╡ 27ffc2b1-9239-4fa4-b575-bab651196be6
begin
	length_k_plot = Plots.scatter(
		total_df.name,
		total_df.Length,
		xlabel="K",
		# color=:red,
		ylabel=" Количество участков",
		# legend = :outertopright
		legend=false
	)
	r_squared_k_plot = Plots.scatter(
		total_df.name,
		total_df.R2,
		xlabel="K",
		color=:red,
		ylabel="R^2",
		# legend = :outertopright
		legend=false
	)
	finish_diff_k_plot = Plots.scatter(
		total_df.name,
		total_df.Finish_diff./start_energy*100,
		xlabel="K",
		color=:yellow,
		ylabel=" Относительная ошибка %",
		# yscale=:log10,
		# legend = :outertopright
		legend=false
		# title=""
	)
	
	plot(length_k_plot, r_squared_k_plot, finish_diff_k_plot, layout=(3,1), size=(900,700), legend=false, left_margin = 20px)
end

# ╔═╡ 5fccfd98-c270-4cb5-b8e4-cce7e0e5f04b
begin
	length_k_plot2 = Plots.scatter(
		total_df.name,
		total_df.Length,
		# xlabel="K",
		# color=:red,
		ylabel="# of segments  ",
		# legend = :outertopright
		legend=false
	)
	r_squared_k_plot2 = Plots.scatter(
		total_df.name,
		total_df.R2,
		# xlabel="K",
		color=:red,
		ylabel="R^2",
		# legend = :outertopright
		legend=false
	)
	finish_diff_k_plot2 = Plots.scatter(
		total_df.name,
		abs.(total_df.Finish_diff)./start_energy*100,
		xlabel="k",
		color=:yellow,
		ylabel="Finish energy error %",
		# yscale=:log10,
		# legend = :outertopright
		legend=false
		# title=""
	)
	
	plot(length_k_plot2, r_squared_k_plot2, finish_diff_k_plot2, layout=(3,1), size=(500,520), legend=false, left_margin = 20px)
end

# ╔═╡ 63c9378a-df62-41d1-b3a0-db5a1be3b3c5
sort!(total_df, [:Length, :name], rev=[true, false]);

# ╔═╡ 0902e42e-8435-4a45-ab7f-9d7969dce188
total_df_short = copy(total_df[total_df.Length .> 500,:])

# ╔═╡ 9e6bd014-b57b-41a4-a170-9e466325d51e
md"А теперь попробуем сделать так, чтобы каждая точка была отдельной серией"

# ╔═╡ 2f655bd9-ad59-4114-95b5-d4fdfb6447e6
Plots.scatter(
	total_df_short.Length',
	total_df_short.RMSE',
	labels=permutedims(total_df_short.name),
	markershapes=:auto,
	xlabel="Количество участков",
	ylabel="Энергия в системе, RMSE",
	size=(400,350)
)

# ╔═╡ e8edab4f-c2cd-41dd-9aa8-2140f1641e56
Plots.scatter(
	total_df_short.Length',
	total_df_short.Finish_diff',
	labels=permutedims(total_df_short.name),
	markershapes=:auto,
	xlabel="Количество участков",
	ylabel="Разница в энрегии на финише",
	size=(400,350),
	legend=:topright,
	ylims=(-0.2, 2.25)
)

# ╔═╡ 4ba9ce08-ef82-415a-bdcc-5c70b8d4a16e
md"Теперь такая же картинка для полноты картины с большим количество thr"

# ╔═╡ b50bbc59-dc26-4cec-b674-6d733f573421
total_df_shorter = copy(total_df[total_df.Length .> 100,:])

# ╔═╡ a7e2adf8-154b-4ad1-bc53-aed37c9846ad
Plots.scatter(
	total_df_shorter.Length',
	total_df_shorter.RMSE',
	labels=permutedims(total_df_shorter.name),
	markershapes=:auto,
	xlabel="Track Length, pieces",
	ylabel="RMSE",
)

# ╔═╡ d679e948-8d30-4620-a4fb-f202c04117b8
Plots.scatter(
	total_df_shorter.Length',
	total_df_shorter.Finish_diff',
	labels=permutedims(total_df_shorter.name),
	markershapes=:auto,
	xlabel="Track Length, pieces",
	ylabel="Energy diff on finish",
)

# ╔═╡ b018504c-bcfc-4c01-9b05-43936b43b615
md"Видно, что брать threshold больше 2.5 смысла не имеет, т.к. начинается значительная потеря точности"

# ╔═╡ 4f2abf90-3e70-422a-9a1f-a85083269a9a
md"Можно взять вариант peaks, как наиболее логичный"

# ╔═╡ 3304d896-71d8-4e81-9fb9-eb84a63d52c3
md"### А теперь для трассы высотой в 10 раз больше"

# ╔═╡ ca2d64b3-4398-4f5d-a5c0-70d73999ae02
begin
	track_high = copy(track)
	track_high.altitude = track_high.altitude .* 10;
	segments_high = get_segments_for_track(track_high);
end

# ╔═╡ 2d082f91-5208-45ac-93e4-f45ce7d5fe66
md"Сперва надо узнать какая скорость оптимальная"

# ╔═╡ d7b60c09-d97f-4a78-9c1e-98dbf9da1372
opt_speed_high = 29.5

# ╔═╡ 407f122d-01ae-4e6b-9ae0-b25b15cf23f8
simulate_run_finish_time(
	fill(opt_speed_high, size(segments_high, 1)),
	track_high,
	segments_high,
	start_energy,
	start_datetime
)

# ╔═╡ e7134856-8988-4ae9-bec2-53ab57124541
thr_df_high_large = make_comparison_df(track_high, segments_high, .0:.5:20, opt_speed_high, start_energy, start_datetime)

# ╔═╡ b47c0f30-933c-47a4-8691-46c4d3d48ed7
md"Прям сильно отличается. Здесь уже до k=4.5 можно упрощать"

# ╔═╡ b6887cf2-dfd0-43f4-860f-31b18935caef
md"## Сравнение самих трасс"

# ╔═╡ a8b3097c-d965-435e-8712-2df777def9b6
# исходная трасса
plot(track.distance, track.altitude, title="Track data")

# ╔═╡ 342e32ee-5343-4b4b-b008-9f2cdb44e9b6
md"Сперва найдём какие-нибудь точки, от и до которых будем смотреть трассу (по длине)" 

# ╔═╡ 4efd0ae0-1322-4488-b80b-d7c5aa067d9d
from_distance = 500

# ╔═╡ 8cb42e3c-b5b9-4d0e-ae84-f7b9304920b0
to_distance = 2000

# ╔═╡ cdd5190e-1cb4-4cc6-96af-09d149a06f82
track_dist = track[from_distance .< track.distance .< to_distance,:]

# ╔═╡ c41e20bf-70b3-433d-add3-97ed29c98c84
plot(track_dist.distance, track_dist.altitude, title="Track data from to")

# ╔═╡ 2ff0a042-f411-4994-b777-fa9bb2ce8c6b
size(track, 1)

# ╔═╡ a0a42d56-1dc3-4701-80ee-ff80078e7734
size(track_peaks, 1)

# ╔═╡ ccef9339-201e-4ea3-9eda-06972be2919a
size(track, 1) / size(track_peaks, 1)

# ╔═╡ a8f4478d-37b8-4483-b8f1-d541f65f7b5d
track_peaks_dist = track_peaks[from_distance .< track_peaks.distance .< to_distance,:]

# ╔═╡ 5f4819b9-3b5b-4700-8b1f-230e38a1cf28
plot(track_peaks_dist.distance, track_peaks_dist.altitude, title="Track peaks data from to")

# ╔═╡ ed85fb2c-b35b-42fa-8933-b740ceb8f392
track_1_75, points_1_75 = parametrized_track_simplification(track, 1.75)

# ╔═╡ 5d6f15de-3c56-446b-8afa-9a688804ce4e
md"Всего 1032 точки, как-то маловато

Надо проверить, нормальные ли получаются результаты" 

# ╔═╡ 092a1984-5515-4de6-b936-06344281d932
draw_comparison_for_parametrized_peaks([0., 1.75], track, segments)

# ╔═╡ f6bea570-3934-4c3b-8a73-d196e105479b
draw_comparison_for_parametrized_peaks([0.,2.], track, segments)

# ╔═╡ d3937387-d7df-4228-8f3d-ff325a1b3d53
segments_1_75 = get_segments_for_track(track_1_75);

# ╔═╡ 76428346-cbff-4745-afd0-31a479ca2494
compare_track_energies_plot_thr(track, segments, opt_speed, 1.75)

# ╔═╡ ffe7d1a9-6c55-45a1-abd6-c4610d922dd0
compare_track_energies_plot_thr(track, segments, opt_speed, 2.25)

# ╔═╡ 8c1bc363-0825-496a-bc35-d03d34a54220
md"Нормально выходит, погрешность меньше процента для 1.75"

# ╔═╡ fad5c87a-2611-40ea-a6d7-a98975ddd7ed
track_1_75_dist = track_1_75[from_distance .< track_1_75.distance .< to_distance,:]

# ╔═╡ 12ec6787-9a6b-4797-9296-fd8cbeb3c1b0
plot(track_1_75_dist.distance, track_1_75_dist.altitude, title="Track 1.75 data from to")

# ╔═╡ e2342368-eaa0-41e1-be5f-0cfce2403c34
begin
	xlim_left=3000
	xlim_right=4000
	ylim_bottom=minimum(track[xlim_left .< track.distance .< xlim_right, :].altitude)
	ylim_top=maximum(track[xlim_left .< track.distance .< xlim_right, :].altitude)
	original_track_plot = plot(
		track.distance, track.altitude,
		title="Исходное представление маршрута",
		# xlabel="Дистанция (м)",
		ylabel="Высота (м)",
		legend=false,
		xlim=[xlim_left, xlim_right], 
		ylim=[ylim_bottom, ylim_top]
	)

	peaks_track_plot = plot(
		track_peaks.distance, track_peaks.altitude,
		title="Сокращённое представление, по подъёмам",
		# xlabel="Дистанция (м)",
		ylabel="Высота (м)",
		legend=false,
		xlim=[xlim_left, xlim_right], 
		ylim=[ylim_bottom, ylim_top]
	)

	k_1_75_track_plot = plot(
		track_1_75.distance, track_1_75.altitude,
		title="Сокращённое представление, параметрическое, k=1.75",
		xlabel="Дистанция (м)",
		ylabel="Высота (м)",
		legend=false,
		xlim=[xlim_left, xlim_right], 
		ylim=[ylim_bottom, ylim_top]
	)
	
	plot(original_track_plot, peaks_track_plot, k_1_75_track_plot, layout=(3,1),
		size=(1000,700),
		legend=false,
		left_margin=20px,
		right_margin=20px
	)
end

# ╔═╡ 91fc63c9-cab2-4151-be14-f856ffdeb201
begin
	# xlim_left=3000
	# xlim_right=4000
	# ylim_bottom=minimum(track[xlim_left .< track.distance .< xlim_right, :].altitude)
	# ylim_top=maximum(track[xlim_left .< track.distance .< xlim_right, :].altitude)
	original_track_plot2 = plot(
		track.distance, track.altitude,
		title="Initial track representation",
		# xlabel="Дистанция (м)",
		# ylabel="Altitude (meters)",
		legend=false,
		xlim=[xlim_left, xlim_right], 
		ylim=[ylim_bottom, ylim_top]
	)

	peaks_track_plot2 = plot(
		track_peaks.distance, track_peaks.altitude,
		title="Extemum points only",
		# xlabel="Дистанция (м)",
		ylabel="Altitude (meters)",
		legend=false,
		xlim=[xlim_left, xlim_right], 
		ylim=[ylim_bottom, ylim_top]
	)

	k_1_75_track_plot2 = plot(
		track_1_75.distance, track_1_75.altitude,
		title="Parametric merging, k=1.75",
		xlabel="Distance (meters)",
		# ylabel="Altitude (meters)",
		legend=false,
		xlim=[xlim_left, xlim_right], 
		ylim=[ylim_bottom, ylim_top]
	)
	
	plot(original_track_plot2, peaks_track_plot2, k_1_75_track_plot2, layout=(3,1),
		size=(600,400),
		legend=false,
		left_margin=20px,
		right_margin=20px
	)
end

# ╔═╡ cb6f30e9-bd84-494a-856e-7dc284124a22
md"# Влияние длины трассы на производительность"

# ╔═╡ 66c7f41c-ba84-4522-9b9a-27b6cc3a6893
md"## Для симуляции"

# ╔═╡ d08e325e-5d2a-486c-87e5-cb0833e3882b
@time simulate_run_finish_time(
	fill(opt_speed, size(segments_dict[0.05], 1)),
	tracks[0.05],
	segments_dict[0.05],
	start_energy,
	start_datetime
)

# ╔═╡ 53ce657b-1dd9-49c6-a152-fd3358d17b8f
function make_comparison_times_df(track, segments, thresholds, speed, start_energy, start_datetime)
	# 1. reference result
	full_energy, full_time = simulate_run_energies(speed, track, segments, start_energy, start_datetime)
	res_df = DataFrame(Threshold=Float64[], Finish_diff=Float64[], MAE=Float64[], MSE=Float64[], RMSE=Float64[], R2=Float64[], Length=Int64[], MedTime=Float64[], MeanTime=Float64[])

	N = length(thresholds);
	# p = Progress(N);
	# update!(p,0)
	jj = Threads.Atomic{Int}(0)
	l = Threads.SpinLock()
	# @withprogress name="iterating" begin
	# @Threads.threads for thr in thresholds
	for thr in thresholds
		track_thr, points_thr = parametrized_track_simplification(track, thr);
		segments_thr = get_segments_for_track(track_thr);
		# t_start = time()
		reduced_energy, reduced_time = simulate_run_energies(speed, track_thr, segments_thr, start_energy, start_datetime)
		bench_info = @benchmark tmp1, tmp2 = simulate_run_energies($speed, $track_thr, $segments_thr, $start_energy, $start_datetime) samples=50 evals=10
		# t_finish = time()
		# exec_time = t_finish - t_start
		med_time = median(bench_info.times)
		mean_time = mean(bench_info.times)
		energy_diff = full_energy[points_thr] - reduced_energy;
		source = full_energy[points_thr];
		new_energy = reduced_energy;
		# это вся разница
		# считаем от неё метрики: r^2, MSE, RMSE разница на финише
		# sqrt(sum((arr1 .- arr2).^2) / length(arr1)) 
		mae_val = mae(source, new_energy)
		mse_val = mse(source, new_energy)
		rmse_val = sqrt(mse_val)
		r2_val = r2_score(source, new_energy)
		last_diff = last(source) - last(new_energy)
		number_of_segments = length(points_thr)

		# Threads.atomic_add!(jj, 1)
		# Threads.lock(l)
		push!(res_df, (thr, last_diff, mae_val, mse_val, rmse_val, r2_val, number_of_segments, med_time, mean_time))
		# index = findfirst(x -> x==thr, thresholds)
		# @logprogress index/length(thresholds)
		
		# update!(p, jj[])
		# @logprogress jj[]/N
		# Threads.unlock(l) 
		
		# energies_thr = compare_track_energies(
		# 	track, segments,
		# 	track_thr, segments_thr, points_thr,
		# 	speed, start_energy, start_datetime
		# )
	end
	# end

	

	# посчитать в отдельной функции
	# # добавить результаты обычного peaks (не параметрического)
	# track_peaks, segments_peaks, points_peaks2 = keep_extremum_only_peaks_segments_with_points(track);
	# energies_peaks = compare_track_energies(
	# 	track, segments,
	# 	track_peaks, segments_peaks, points_peaks2,
	# 	speed, start_energy, start_datetime
	# )

	# # добавить результаты нового peaks (параметрического)
	# track_mod_peaks, segments_mod_peaks = get_track_and_segments_for_selected_points_modified(track, points_peaks2)
	# energies_mod_peaks = compare_track_energies(
	# 	track, segments,
	# 	track_mod_peaks, segments_mod_peaks, points_peaks2,
	# 	speed, start_energy, start_datetime
	# )
	
	return res_df
end

# ╔═╡ 95f02f1c-914b-4363-a100-eb34512ee911
thr_df_large_times = make_comparison_times_df(track, segments, .0:.25:5, opt_speed, start_energy, start_datetime)

# ╔═╡ df9ec094-a2ae-428d-aacb-f53a91c4eff8
Plots.scatter(
	thr_df_large_times.Length,
	[
		thr_df_large_times.MedTime thr_df_large_times.MeanTime
	],
	labels=["Median time" "Mean time"],
	xlabel="Length (segments)",
	ylabel="Time (ns)"
)

# ╔═╡ f1dd564a-1a34-4945-bba4-f142760ea533
mean_linear_fit = curve_fit(Polynomial, thr_df_large_times.Length, thr_df_large_times.MedTime, 1)

# ╔═╡ 283b19e1-2b1b-4d7f-9e9a-3aaeb0eca290
mean_fit_sim = curve_fit(Polynomial, thr_df_large_times.Length, thr_df_large_times.MeanTime, 1)

# ╔═╡ 87d04f23-87cc-4c8d-a364-57e221ace9ce
plot(
	thr_df_large_times.Length,
	[ thr_df_large_times.MeanTime/1e9  mean_fit_sim.(thr_df_large_times.Length)/1e9],
	# labels=["Mean" "Median" "Mean (fit):"*text(mean_fit).str "Median (fit):"*text(median_fit).str],
	labels=["Среднее время по замеру" "Аппроксимация:"*text(mean_fit_sim).str],
	xlabel="Количество участков",
	ylabel="Время (с)",
	seriestypes=[:scatter :path ],
	title="Время моделирования"
)

# ╔═╡ 23dffd9a-4fd3-41ce-9645-943b6166d7f2
Plots.scatter(
	thr_df_large_times.Length,
	[
		thr_df_large_times.MeanTime/1e9
	],
	label="Время моделирования",
	xlabel="Количество участков",
	ylabel="Время (с)"
)

# ╔═╡ 71463a7f-35e0-4196-aa3c-0a0df5734696
Plots.scatter(
	thr_df_large_times.Length,
	[
		thr_df_large_times.MeanTime/1e9
	],
	label="Modeling time",
	xlabel="Amount of segments",
	ylabel="Time (seconds)"
)

# ╔═╡ b59f861e-8ced-4900-9fdf-d59eef82e9e3
md"Линейная зависимость скорости симуляции"

# ╔═╡ 5ef24d5a-e279-4ba0-83d0-4ac348df6290
md"Теперь надо посчитать фактический выигрыш"

# ╔═╡ e426497b-55c1-4792-803f-012636c7559f
thr_df_large_times[thr_df_large_times.Threshold .== 0.0, :].MedTime / thr_df_large_times[thr_df_large_times.Threshold .== 1.75, :].MedTime

# ╔═╡ 20f36ba5-02a8-4883-b9ad-e9339b139b73
md"40.25 раз, абсолютно линейный рост производительности, ровно во сколько сократилось число участков" 

# ╔═╡ 84b3df87-33ae-4c0f-836d-a8b225e52a4e
md"Теперь добавим замер с peaks" 

# ╔═╡ 24ee0cfe-af10-4287-8265-4638ba20d042
bench_peaks_info = @benchmark tmp1, tmp2 = simulate_run_energies(opt_speed, track_peaks, segments_peaks, start_energy, start_datetime) samples=50 evals=10

# ╔═╡ bd90358f-6190-4358-a25f-1a3a85a1a944
med_peaks_time = median(bench_peaks_info.times)

# ╔═╡ 30503dc1-e9e7-4876-9bbe-df8ddbb50c93
thr_df_large_times[thr_df_large_times.Threshold .== 0.0, :].MedTime / med_peaks_time

# ╔═╡ 5955727b-40f1-4eb4-b868-2ae7788c1625
md"## Для оптимизации"

# ╔═╡ 44fa6c8e-9e64-4262-9324-272a133e302e
md"Будем проводить не на полном размере трассы, а на её начале. Т.е. сперва трасса из одного участка, потом из 2-х, 3-х и т.д.

Затем будет аппроксимация сложности вычислений"

# ╔═╡ 0a4a6191-c3dd-4551-b104-dd0595dd684e
md"Проводить будем 2 эксперимента:
1. Число переменных=число сегментов (берём первые n сегментов трассы, сегменты 1:1)
2. Число переменных!=число сегментов (берём всю трассу, делим на n сегментов) "

# ╔═╡ e7c9704d-e310-4b11-aa88-f41d6692635e
md"## Короткие трассы (по числу сегментов)"

# ╔═╡ 1a956d44-14c6-41a9-8015-3691051d0cb6
md"Надо сделать:

1. Сделать benchmark на одной оптимизации
2. Обернуть в цикл и записывать результаты в датафрейм"

# ╔═╡ daa0b062-f7e4-4684-8a42-d85bc9043b71
md"### Бенчим одну оптимизацию"

# ╔═╡ f4d510db-3b37-4f06-bc27-99447fd4fefe
@bind short_track_segments confirm(NumberField(1:size(segments_peaks_pl,1), default=5))

# ╔═╡ c8117fcf-3dfa-4064-8adb-509b17e0280b
segments_short = segments_peaks_pl[1:short_track_segments,:]

# ╔═╡ 80a489d8-277d-412d-92ad-2c54dc98fbdc
track_short = track_peaks_pl[1:short_track_segments+1,:]

# ╔═╡ 2639bd8b-d88c-42f4-b34c-fd206525d3c3
start_energy_short = start_energy * last(track_short.distance) / last(track.distance)

# ╔═╡ b6d30924-90d8-4155-882e-8a72e0c7e797
function f_wrap_short_track(input_speeds)
	speeds_ms = convert_kmh_to_ms(input_speeds)
	power_use_short, solar_power_short, time_s_short = solar_trip_boundaries(
		speeds_ms, segments_short, start_datetime
	)

	energy_in_system = start_energy_short .+ solar_power_short .- power_use_short

	# energy_capacity = start_energy_short

	cost = sum(segments_short.diff_distance ./ speeds_ms) + 100 * (0. - last(energy_in_system))^2;

	# cost = last(time_s) + (
	# 	10000 * (finish_energy - last(energy_in_system))^2 +
	# 	100 * max(0, maximum(energy_in_system) - energy_capacity)
	# )
	return cost
end

# ╔═╡ 2321ad16-58d0-4332-9108-dfde9e6ec6cf
init_speeds_short = fill(30., short_track_segments)

# ╔═╡ 3ddfebdc-7b41-4a17-8cd1-9169ad4ecf07
begin
	td_short = TwiceDifferentiable(f_wrap_short_track, init_speeds_short; autodiff = :forward)
	lower_bound_short = fill(0.0, short_track_segments)
	upper_bound_short = fill(100.0, short_track_segments)
	tdc_short = TwiceDifferentiableConstraints(lower_bound_short, upper_bound_short)
end

# ╔═╡ 16887603-f096-46fe-9fa1-601d1dc2605d
@time res_short = optimize(td_short, tdc_short, init_speeds_short 
# .+ rand(vars_amount) .- 0.5
    ,
    IPNewton(),
    Optim.Options(
        x_tol = 1e-6,
        f_tol = 1e-6,
        g_tol = 1e-6
    )
)

# ╔═╡ ed06d8f1-a22e-4da2-a31c-50ae7176b31d
speeds_short = Optim.minimizer(res_short)

# ╔═╡ e5d77450-c22b-4f63-8adb-0dc1854e3b8c
simulate_run(
	speeds_short,
	track_short,
	segments_short,
	start_energy_short,
	start_datetime
)

# ╔═╡ 05599d8f-82bd-4e3c-98e1-42473e7887ad
md"Получается 80 милисекунд (для 50 участков), что-то маловато

Скорее всего что-то не то с измерением.

Либо только первый прогон длинный почему-то

Первый длинный, второй и далее - быстрые"

# ╔═╡ 9f5f1f2c-498f-460b-b6b4-a58b5c9f5804
md"1900 секунд для трассы в 1000 кусочков"

# ╔═╡ de4dee93-24e5-463b-96fd-078f1ac8f9c9
md"### Бенчим одну быструю оптимизацию"

# ╔═╡ 65fceae5-cc43-4208-b6c5-ae77fdc238ea


# ╔═╡ 08e0976a-7652-4b8b-ad07-c0b8505e086c
md"### Бенчим в цикле"

# ╔═╡ 34f706ea-918d-43fc-a653-92d311b7df18
short_segments_sizes = 10:10:500

# ╔═╡ 01ca48c3-e201-4ad5-87fe-cbec83b7a009
md"Многопоточное исполнение коверкает результаты, либо я его неправильно делаю. Поэтому тесты в однопотоке, зато с прогрессом)"

# ╔═╡ 1cc1735a-7b11-4ecf-9b11-d0aac8a79c7a
function benchmark_optimizations(segment_sizes::Vector{<:Number}, track, segments, start_energy, start_datetime)
	res_df = DataFrame(Length=Int64[], Mean=Float64[], Median=Float64[])
	# @progress @Threads.threads for size in segment_sizes
	@progress for i in 1:length(segment_sizes)
		size = segment_sizes[i]
		track_short = track[1:(size+1),:]
		segments_short = segments[1:size,:]
		start_energy_short = start_energy * last(track_short.distance) / last(track.distance)

		function f_wrap_short_track_inner(input_speeds)
			speeds_ms = convert_kmh_to_ms(input_speeds)
			power_use_short, solar_power_short, time_s_short = solar_trip_boundaries(
				speeds_ms, segments_short, start_datetime
			)
		
			energy_in_system = start_energy_short .+ solar_power_short .- power_use_short
		
			# energy_capacity = start_energy_short
		
			cost = sum(segments_short.diff_distance ./ speeds_ms) + 100 * (0. - last(energy_in_system))^2;
		
			# cost = last(time_s) + (
			# 	10000 * (finish_energy - last(energy_in_system))^2 +
			# 	100 * max(0, maximum(energy_in_system) - energy_capacity)
			# )
			return cost
		end

		init_speeds_short = fill(30., size)

		td_short = TwiceDifferentiable(f_wrap_short_track_inner, init_speeds_short; autodiff = :forward)
		lower_bound_short = fill(0.0, size)
		upper_bound_short = fill(100.0, size)
		tdc_short = TwiceDifferentiableConstraints(lower_bound_short, upper_bound_short)

		bench_info_short = @benchmark res_short = optimize($td_short, $tdc_short, $init_speeds_short 
		# .+ rand(vars_amount) .- 0.5
		    ,
		    IPNewton(),
		    Optim.Options(
		        x_tol = 1e-6,
		        f_tol = 1e-6,
		        g_tol = 1e-6
		    )
		) samples=3 evals=1 seconds=10
		
		med_time = median(bench_info_short.times) * 1e-6 # to ms
		mean_time = mean(bench_info_short.times) * 1e-6
		push!(res_df, (size, mean_time, med_time))
	end
	return res_df
end

# ╔═╡ 8a195f8f-c9ff-4f24-aca1-5b85d7421383
optim_times_df = CSV.read("optim_times.csv", DataFrame)

# ╔═╡ ec3eb92f-9232-4309-8564-d8c9b620293b
mean_fit = curve_fit(Polynomial, optim_times_df.Length, optim_times_df.Mean, 3)

# ╔═╡ 28c338dd-6a85-49ad-ac5b-e4df11210abd
mean_fit2 = curve_fit(Polynomial, optim_times_df.Length, optim_times_df.Mean, 2)

# ╔═╡ 4592706f-7f90-41c4-bf11-4254fa6450a2
median_fit = curve_fit(Polynomial, optim_times_df.Length, optim_times_df.Median, 3)

# ╔═╡ c4c9f140-086a-483e-9bb8-5af8e452181e
plot(
	optim_times_df.Length,
	[ optim_times_df.Mean/1000.  mean_fit.(optim_times_df.Length)/1000.],
	# labels=["Mean" "Median" "Mean (fit):"*text(mean_fit).str "Median (fit):"*text(median_fit).str],
	labels=["Среднее время по замеру" "Аппроксимация:"*text(mean_fit).str],
	xlabel="Количество участков",
	ylabel="Время (с)",
	seriestypes=[:scatter :path ],
	title="Время оптимизации"
)

# ╔═╡ 2fcdf706-5ff6-4e86-a6e9-839e9144655a
plot(
	optim_times_df.Length,
	[ optim_times_df.Mean/1000.  mean_fit2.(optim_times_df.Length)/1000.],
	# labels=["Mean" "Median" "Mean (fit):"*text(mean_fit).str "Median (fit):"*text(median_fit).str],
	labels=["Среднее время по замеру" "Аппроксимация:"*text(mean_fit2).str],
	xlabel="Количество участков",
	ylabel="Время (с)",
	seriestypes=[:scatter :path ],
	title="Время оптимизации"
)

# ╔═╡ 67de2dbd-3633-4728-a6e4-53e2ddfb2625
plot(
	optim_times_df.Length,
	[ optim_times_df.Mean optim_times_df.Median mean_fit.(optim_times_df.Length) median_fit.(optim_times_df.Length)],
	# labels=["Mean" "Median" "Mean (fit):"*text(mean_fit).str "Median (fit):"*text(median_fit).str],
	labels=["Mean" "Median" "Mean (fit)" "Median (fit)"],
	xlabel="Length",
	ylabel="Time(ms)",
	seriestypes=[:scatter :scatter :path :path]
)

# ╔═╡ e9821459-f1b6-4884-8d0c-3017ac2d776d
md"Выглядит как квадратичная зависимость, составляющая куба очень мала"

# ╔═╡ 5c493ed5-e368-4550-88da-bf6c7338b879
mean_fit(1000)

# ╔═╡ da45b3b3-c3c2-4fe8-a0ba-2530ec2a6e84
md"Должно быть 1900 секунд, а не 661"

# ╔═╡ 89641cc8-72b2-4a4a-9fcb-e48bdf7a778b
md"Посмотрим, что будет на полной трассе" 

# ╔═╡ f5419815-56ef-4485-b72f-ee3bf2d3461b
mean_fit(
	size(
		segments_peaks_pl,
		1
	)
) / 1000. / 3600.

# ╔═╡ 2a341070-a64e-41dc-b41b-c8e68730dc9f
md"Вышло 771 час"

# ╔═╡ f33a007e-2e91-4206-bb80-65d1997247a3
md"И для вообще полной дистанции"

# ╔═╡ 4ad15809-2bac-4c25-8139-66fb55311f9c
mean_fit(
	size(
		segments,
		1
	)
) / 1000. / 3600.

# ╔═╡ 85376195-8c5e-42a1-a307-a7f0721556b9
md"Плоховато сходится, надо больше измерений. Запустить на ночь с шагом в 100-200 до 1.5 тысяч"

# ╔═╡ d04b0d6a-96ba-4aba-b943-74cf46f8855c
md"В любом случае, какая-то зависимость есть, теперь пора делать для полной трассы"

# ╔═╡ c99b3c6a-c29d-475d-bb79-ff23d505fe0a
md"## Вся трасса с разным числом сегментов"

# ╔═╡ 7b6dd0b4-6e9c-44bf-baf7-6dcc0709086d
md"Здесь не будем даже делать без цикла, сразу в нём" 

# ╔═╡ 1938c38a-5dfa-43a2-8810-cc987a34cce4
segments_amount = 200:10:200

# ╔═╡ cdefdb31-f6fe-4f86-9f12-e01f068f1108
function benchmark_optimizations_whole(segment_sizes, track, segments, start_energy, start_datetime)
	res_df = DataFrame(Length=Int64[], Mean=Float64[], Median=Float64[])
	# @progress @Threads.threads for size in segment_sizes
	@progress for i in 1:length(segment_sizes)
		splits_amount = segment_sizes[i]
		# track_short = track[1:(size+1),:]
		# segments_short = segments[1:size,:]
		# start_energy_short = start_energy * last(track_short.distance) / last(track.distance)

		boundaries = calculate_boundaries(
			1,
			size(track, 1),
			splits_amount
		)

		function f_wrap_track_inner(input_speeds)
			speeds_ms = convert_kmh_to_ms(input_speeds)
			speed_vector = set_speeds_boundaries(speeds_ms, boundaries)
			power_use, solar_power, time_s = solar_trip_boundaries(
				speed_vector, segments, start_datetime
			)
		
			energy_in_system = start_energy .+ solar_power .- power_use
		
			# energy_capacity = start_energy_short
		
			cost = sum(segments.diff_distance ./ speed_vector) + 100 * (0. - last(energy_in_system))^2;
		
			# cost = last(time_s) + (
			# 	10000 * (finish_energy - last(energy_in_system))^2 +
			# 	100 * max(0, maximum(energy_in_system) - energy_capacity)
			# )
			return cost
		end

		init_speeds = fill(30., splits_amount)

		td = TwiceDifferentiable(f_wrap_track_inner, init_speeds; autodiff = :forward)
		lower_bound = fill(0.0, splits_amount)
		upper_bound = fill(100.0, splits_amount)
		tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)


		
		# bench_info = @benchmark res = optimize($td, $tdc, $init_speeds 
		# # .+ rand(vars_amount) .- 0.5
		#     ,
		#     IPNewton(),
		#     Optim.Options(
		#         x_tol = 1e-6,
		#         f_tol = 1e-6,
		#         g_tol = 1e-6
		#     )
		# ) samples=3 evals=1 seconds=10

		t_start = time()
		res = optimize(td, tdc, init_speeds 
		# .+ rand(vars_amount) .- 0.5
		    ,
		    IPNewton(),
		    Optim.Options(
		        x_tol = 1e-6,
		        f_tol = 1e-6,
		        g_tol = 1e-6
		    )
		)
		t_finish = time()
		exec_time = t_finish - t_start
		
		# med_time = median(bench_info.times) * 1e-6 # to ms
		# mean_time = mean(bench_info.times) * 1e-6
		med_time = exec_time
		mean_time = exec_time
		push!(res_df, (splits_amount, mean_time, med_time))
	end
	return res_df
end

# ╔═╡ aa0db854-c097-4df4-b139-058869a6be8f
optim_times_whole_df = CSV.read("optim_times_whole.csv", DataFrame)

# ╔═╡ fca157a2-743b-452f-b035-44fe7363764b
md"Очень долго считается

Значит мой speed propagation уж очеень сильно вредит оптимизации (что логично). "

# ╔═╡ 53a4439c-3714-4c95-a0db-9f4ca9d1105b
opt_time_fit = curve_fit(Polynomial, optim_times_whole_df.Length, optim_times_whole_df.OptTime, 2)

# ╔═╡ 3188f726-5bd3-46d7-bb7a-760ed1fcff05
td_time_fit = curve_fit(Polynomial, optim_times_whole_df.Length, optim_times_whole_df.TdTime, 3)

# ╔═╡ 95c2de46-2634-4ab3-a7c8-d4251cdabd32
plot(
	optim_times_whole_df.Length,
	[ optim_times_whole_df.OptTime  opt_time_fit.(optim_times_whole_df.Length)],
	# labels=["Mean" "Median" "Mean (fit):"*text(mean_fit).str "Median (fit):"*text(median_fit).str],
	labels=["Data" "Fit"],
	xlabel="Length",
	ylabel="Time(s)",
	seriestypes=[:scatter :path ]
)

# ╔═╡ 51a64f33-e0e9-4365-bab8-6d134f27517d
plot(
	optim_times_whole_df.Length,
	[ optim_times_whole_df.OptTime optim_times_whole_df.TdTime opt_time_fit.(optim_times_whole_df.Length) td_time_fit.(optim_times_whole_df.Length)],
	# labels=["Mean" "Median" "Mean (fit):"*text(mean_fit).str "Median (fit):"*text(median_fit).str],
	labels=["OptTime" "TdTime" "OptTime (fit)" "TdTime (fit)"],
	xlabel="Length",
	ylabel="Time(s)",
	seriestypes=[:scatter :scatter :path :path]
)

# ╔═╡ d154ace4-9a83-4b9a-8dce-3b2d1b569200
md"Получается очень медленно с распределением переменных"

# ╔═╡ f6b131d4-05eb-4005-b55c-5fad152f5090
md"Пробуем сколько примерно выйдет для полного размера задачи"

# ╔═╡ a8f062ed-659f-49d5-9cac-cdc47393e90a
full_time = opt_time_fit(
	size( # количество участков на трассе
		segments_peaks_pl,
		1
	)
)

# ╔═╡ 815e44dd-57cf-4104-95d9-1632ddc4ead4
full_time / 3600.

# ╔═╡ 4a69919b-38e3-4a05-9299-8e6c50b5bd25
md"Получается 1869 часов для полного размера трассы"

# ╔═╡ 80a26e90-77aa-4cb1-88de-856e8f366453
md"НО! Получается довольно быстро для нормальной оптимизации до примерно 200 переменных (но это трасса короткая)" 

# ╔═╡ 765971bf-1ff6-453e-b90c-a00126c2095c
md"Что это значит - использовать подстановку переменных только для первых итераций"

# ╔═╡ 1b4d58ea-ae61-4994-95c4-00fdfb2f4b75
md"Надо подумать, как изменить код итеративной оптимизации, чтобы работало быстрее" 

# ╔═╡ be87100b-a5c7-4050-b4ab-ad065b97091a


# ╔═╡ 0a7000be-88dd-4cba-a0b8-6f82ff1c9c07
md"# На будущее"

# ╔═╡ 308ea3b1-8145-49a7-aa16-048666312620
md"Идеи на будущее:

1. Сглаживание трассы
2. Убирание слишком коротких кусков?"

# ╔═╡ d76052a6-92a3-416e-89f7-f37e160b7d54


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CurveFit = "5a033b19-8c74-5913-a970-47c3779ef25c"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Metrics = "cb9f3049-315b-4f05-b90c-a8adaec4da78"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlotlyBase = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"

[compat]
BenchmarkTools = "~1.3.2"
CSV = "~0.10.11"
CurveFit = "~0.5.0"
DataFrames = "~1.6.1"
Metrics = "~0.1.2"
Optim = "~1.7.7"
Peaks = "~0.4.4"
PlotlyBase = "~0.8.19"
Plots = "~1.38.17"
PlutoUI = "~0.7.52"
ProgressLogging = "~0.1.4"
ProgressMeter = "~1.7.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "0d26fd02ca2063f5e9f3aa2cd9654fd748aba60f"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "44dbf560808d49041989b8a96cae4cffbeb7966a"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.11"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "d9a8f86737b665e15a9641ecbac64deef9ce6724"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.23.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fe2838a593b5f776e1597e086dcd47560d94e816"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.3"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.CurveFit]]
deps = ["LinearAlgebra", "Polynomials"]
git-tree-sha1 = "074cb8efc989bfd1ba869160889b15037560a341"
uuid = "5a033b19-8c74-5913-a970-47c3779ef25c"
version = "0.5.0"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

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
git-tree-sha1 = "f372472e8672b1d993e93dada09e23139b509f9e"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.5.0"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "d73afa4a2bb9de56077242d98cf763074ab9a970"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.9"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1596bab77f4f073a14c62424283e7ebff3072eca"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.9+1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

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
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "cb56ccdd481c0dd7f975ad2b3b62d9eda088f7e2"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.14"

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
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

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
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

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

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

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
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

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
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

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
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "5ab83e1679320064c29e8973034357655743d22d"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.25"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

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
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Metrics]]
deps = ["DataFrames", "DataStructures", "Random", "StatsBase"]
git-tree-sha1 = "6e9e77751dd230b360c29e23a10f6e6d2f4fafaf"
uuid = "cb9f3049-315b-4f05-b90c-a8adaec4da78"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

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
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "bbb5c2115d63c2f1451cb70e5ef75e8fe4707019"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.22+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "963b004d15216f8129f6c0f7d187efa136570be0"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.7"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Peaks]]
deps = ["Compat", "RecipesBase"]
git-tree-sha1 = "1627365757c8b87ad01c2c13e55a5120cbe5b548"
uuid = "18e31ff7-3703-566c-8e60-38913d67486b"
version = "0.4.4"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "9f8675a55b37a70aa23177ec110f6e3f4dd68466"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.17"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3aa2bb4982e575acd7583f01531f241af077b163"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.13"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

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

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "ee094908d720185ddbdc58dbe0c1cbe35453ec7a"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "364898e8f13f7eaaceec55fd3d08680498c0aa6e"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.4.2+3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "04bdff0b09c65ff3e06a05e3eb7b120223da3d39"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.0"

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

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

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

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

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

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "607c142139151faa591b5e80d8055a15e487095b"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.16.3"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

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

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

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
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

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
version = "1.52.0+1"

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
# ╠═cff82de0-3d15-11ee-37ff-035acc6e8e4a
# ╠═9ce6e2c1-69e4-47c6-830c-3c0eb13d784e
# ╠═7c89a14e-f6c7-4a09-8825-a81d1ba9b824
# ╠═39d80a54-3675-47ff-82bc-deb1daf03a59
# ╠═1b605aae-c293-49b0-a06d-fcff2ef3c12a
# ╠═d6fb7d82-6bd9-4da9-b183-67fb16ab9f66
# ╠═aef61af9-31df-4757-95bd-7fdfc0edbc1a
# ╠═dbe59867-4be2-4a6a-9236-35dc96d2d387
# ╠═3327e2e2-4398-4aed-83c1-eaeb7565bd48
# ╠═e53d252c-31b5-430c-a241-b41159243bc7
# ╠═164055ba-18bb-482d-816e-e278fc76a0d0
# ╠═1b6f993d-e4cc-4973-b5cf-269f01c7f409
# ╠═16405fa9-f44a-4b60-9b87-84bfd8aa8ed2
# ╠═7513698b-e1cf-4dc1-a461-ede112a76180
# ╠═57d53d47-733a-4baf-b2c4-f71d97429966
# ╠═c2da7fb6-e6f0-494c-9f6b-6c20e7127b52
# ╠═5f90dc5b-1146-4085-8ff3-4dbe33854aad
# ╠═b5c8451c-03ba-42d2-98aa-be25daa2729b
# ╠═24f87654-2894-42a8-b13a-3d2a8eacc53b
# ╠═51ae1ae6-98cb-4d76-9ce2-a0a95d88478d
# ╠═e9130295-ff5e-4cd5-92f0-844a4dfdcd8e
# ╠═418a5c92-3a0e-4640-93db-24fc0dd98c51
# ╠═7f4daafa-8a3f-46a2-81f4-fd3d6d3ade37
# ╠═1f039515-b3d7-49ae-aef5-f5951934778d
# ╠═5008c2a6-b61b-4150-a8e2-225fcbad8680
# ╠═0c7cbc87-96b7-4855-91e9-860996b2a5ba
# ╠═36d52b29-260c-427d-86f2-5d45d792e10a
# ╠═eb2f6b3b-6bdb-48b6-8c49-eb7f3a579b57
# ╠═b98aed7d-53f3-4ee4-91d6-9e5b3a9fc67b
# ╠═459d3b2c-0954-4bd8-9369-ac294e1adcb7
# ╠═04a05ccd-b5d5-4ae9-934f-b9d2290a2d58
# ╠═316d420f-b698-4d3c-988b-732e55e5089a
# ╠═77aec1eb-7e4c-4c0e-a31f-f8f7036792df
# ╠═277fcf00-4e82-4ddc-a603-a53bd6695cf7
# ╠═2f02f496-5548-46d5-9ba3-bca71b0b22b6
# ╠═3d672dd6-82dd-402f-847c-67795f87ac8a
# ╠═6081d5ae-583d-4572-8930-facd41137bb4
# ╠═b9b6d1ec-208c-4684-a6d4-3e0c89a0277d
# ╠═a23c2215-531a-47f0-8241-6e0ac06f6d46
# ╠═4d48b247-499b-42f6-8e24-085aee654de4
# ╠═611dec1f-6d0e-4b6d-9f1e-a547aec312ac
# ╠═14e1a205-8cd3-4444-957c-6902e7c41d6a
# ╠═466171d6-a749-4d18-af72-099e099dd5bd
# ╠═0131849a-9f90-4097-ab3d-9bd47449bdc2
# ╠═f7717337-d999-4150-b888-fb0644fd00d1
# ╠═5ce6ed99-e0c1-4100-a5a8-b92e2bf10fc0
# ╠═700120b8-78dd-4f11-8b91-e0bbd267be60
# ╠═64f4a5bc-7d13-45eb-a2c8-3325166cc564
# ╠═8f4fc0ed-77e4-4d30-a5f4-d2d6af16bc56
# ╠═9e2eb0bc-0ffa-4e10-be7a-df9aebc6b609
# ╠═4e2f9485-9131-40e6-afa5-e6134a638750
# ╠═aa894177-a0a5-4488-a784-3752edaeddcb
# ╠═97690c57-f630-4762-bdb8-6630cc14c3a7
# ╠═f61f4da9-0689-4752-ae25-27415a95c10d
# ╠═dadc2e8c-3aae-4944-b5e6-cc99b9c4e1f5
# ╠═0aeeb368-896e-4bc8-ad1c-d92b0c324392
# ╠═cdf8794f-67ed-41df-bd9e-b1c57c82133b
# ╠═eb2ddfbc-a8ab-4b88-ac46-37af7c563d42
# ╠═501fc5bb-8f50-4a19-b2a0-3fafe9a7738e
# ╠═322e9eb5-7da5-492f-a44d-2f7ee8c2ffcb
# ╠═fcbeb07f-3c3f-4973-ac1a-f133c7bdd4f3
# ╠═3cd8689b-5f6c-481b-9e5b-ec0356597598
# ╠═35102b6b-1a3f-4c99-b03a-86fb5ad00e14
# ╠═6950f324-b23c-4b49-b391-586bda1b872e
# ╠═a8f3f8a9-673c-497b-a311-89798970f358
# ╠═ce77f956-bd0a-48c1-a2b2-1f26a49877d8
# ╠═425a892c-b550-486e-b376-8566ae086f6a
# ╠═d9a12d12-8c19-4cd3-88ba-847e1c273850
# ╠═f84953c2-bed0-4375-a9d1-0bb646c832e9
# ╠═c26613ea-827d-46b8-9bd2-15ddf6cd885b
# ╠═b845f4f8-22c0-49cb-9580-d51c19e115bd
# ╠═d95e9d51-867e-484a-8006-e1fb2d1e0fbf
# ╠═baa774fc-d955-4b69-8a32-7bf7b6a92997
# ╠═da13556d-62b0-49ce-867f-f64b8f227599
# ╠═ecc4c5e7-7b9f-47b2-9f35-85f6c241f4ea
# ╠═691e1f9a-5d10-475f-9ede-8d687f9d07b1
# ╠═f17bfc9f-01b5-4535-b05e-a4251e6a3594
# ╠═1c54eea3-b89a-4c7c-a7ab-7759f1bd4242
# ╠═89fdba6b-0c3f-4dd2-ae5b-a26dc1a90368
# ╠═0f248c56-5ab5-4183-9260-634be0fc9736
# ╠═3e84943e-1f80-44d6-b3fc-9f7c0553c4c1
# ╠═1bb03b35-2d86-4fc4-a7ed-d2bcf7120d7d
# ╠═67ac9252-66d6-465b-9857-abd0611fb81c
# ╠═950d2611-9730-4552-af0e-72e9cfc801d3
# ╠═5f5d1144-b83b-4c67-a20e-669387b93493
# ╠═a472bc02-13d7-4ce1-b74a-2a783f0a7965
# ╠═ae90c511-3fcd-405f-b6fb-7814da8893aa
# ╠═171cc89f-05a6-49df-bd9a-deab45341778
# ╠═1a4155e3-a791-4f5c-bc46-687dd17c3742
# ╠═aa9edc8e-6357-4c3e-a8c2-f38241de7f38
# ╠═fbdbf4dc-b1b4-4eef-836f-558c276b83d1
# ╠═b031cc68-9306-4265-973b-e427b801f439
# ╠═e6cbc2b6-7c6a-4405-9899-4af0f7ba979c
# ╠═170467ab-6e26-4e24-9e2a-c244bdf08b9d
# ╠═17a4b5b0-3aa7-4f60-b843-df81823fa710
# ╠═4110b150-d44b-4e2d-9769-2edac7c1a81c
# ╠═4036eb31-d5a1-404c-aa5d-ecb26badaa70
# ╠═d7d2142b-b2e6-44e8-aaa9-c293b2a9b70c
# ╠═58436584-cddc-44ce-b19e-10a1fc7873cf
# ╠═ef6dfaa8-987f-4f79-9957-16db862aa1ff
# ╠═0aaabc2f-aac7-4f03-90e3-a42e95f08e62
# ╠═84f39121-69bb-4a04-8ca4-b89af0f6a0f8
# ╠═6c33a6a3-808a-4086-914d-6a23d749ee38
# ╠═6ef045bd-2426-41f6-9f61-4a17f136ac7a
# ╠═4642901d-e71b-41bc-8d34-5e1f01f83306
# ╠═6e697d1a-5823-4f9f-a3e0-c649adc0b449
# ╠═4d6b8b1c-081b-4db2-bede-6d67cb41fbaf
# ╠═28275f01-96c4-4988-8d2c-19648ea388b6
# ╠═9b5b321a-ed7f-4b88-8519-498f2d19fb1c
# ╠═09dddfde-100f-43d9-ad48-65eecb39923b
# ╠═91cae49c-bd36-473f-a23b-d15ecd1ac073
# ╠═dea20b86-88d4-4337-8f58-25606c7adff1
# ╠═9d390be8-b231-4fc8-827a-f561293e82b9
# ╠═1550e304-31d0-4d4c-b5a2-8e4d72561f05
# ╠═55f67256-7b17-40a3-a739-d502ac5864fb
# ╠═1e7a691c-d3cd-4d4c-a202-311e222a6320
# ╠═5da7a24e-b075-4d6f-ab20-e7bfd4a2aed9
# ╠═b62e5863-7bf9-4e96-9891-267dc4baa737
# ╠═9ab0a9a7-c071-49ba-bc0d-077b44c5c041
# ╠═26c9ff7f-2ff7-4fe8-9abe-d72f4d869c5c
# ╠═d6645aef-70a6-4af3-8cb4-4254d80a7281
# ╠═95fd37d2-f7b0-479b-9779-093efcc11066
# ╠═e42e57dd-5513-4ca8-9ddf-5816f94f5707
# ╠═205cc77c-2f01-4552-a660-08bd750d9b43
# ╠═aa4596d6-a281-4ed3-8115-f5f3e30d1f35
# ╠═9af6b676-d183-4f52-b92a-8017a0df4c85
# ╠═aa4c6f7a-d6ca-430e-8297-5e07e9e9e2b3
# ╠═81eac0f6-9c2f-4123-886f-c0317d70444b
# ╠═2790174f-7a2e-46c9-967e-065c4f9a12fd
# ╠═d36a4805-e59b-45ba-9433-54fd26733670
# ╠═8695b7a4-8127-4766-a7f5-fcc7d7f5b331
# ╠═a9bad1e6-7a9c-49c0-a0aa-33aaa6776691
# ╠═0059a826-9daa-4d9c-9213-b1acefb8b92f
# ╠═cbc854d5-d287-495c-9921-408ce23d1f30
# ╠═af6550be-7453-41ef-bcbe-dd41312eec76
# ╠═5b770bfc-64d1-4b8e-b30d-e40cf4eb5397
# ╠═f9dc21cb-d389-454d-99da-383208e268f9
# ╠═76438295-379c-4d08-a57c-ac3b341f3a90
# ╠═51bd9a66-5e73-4fe4-b7cf-47ed699d8011
# ╠═cd865e2f-c24f-4b6d-a8ed-45de9867e7b9
# ╠═ffbb006a-47be-404a-a564-d171e70c0590
# ╠═40e187fb-8f63-4b09-bb0e-110633a1ba27
# ╠═cfb3b295-61a3-4851-b511-6b9d6b162ded
# ╠═256c6718-883a-479b-838a-2e939f438c3b
# ╠═58cc06b0-14d1-476f-a5a4-ca56ef022695
# ╠═697ff289-724d-4f3b-a89d-4ec3a96567bc
# ╠═cde20eed-7581-43f1-9abc-39e6a0555478
# ╠═1133a183-3631-4ef8-b2b8-1b697a29b75f
# ╠═6f342cb3-c8e2-4358-86b9-b45ad05d2551
# ╠═151dfb5e-f19a-4f7e-8764-12f99185094b
# ╠═eede70d8-37e2-4082-a8ae-cd1491e768fa
# ╠═6b27c66e-873b-4ad9-8b9b-8987ee09c03b
# ╠═e593df22-ebbe-4096-894e-4414b4013d0b
# ╠═e46f1433-f940-4b3c-a6c5-49a70c07f671
# ╠═cc352146-3a5e-4431-919e-649995ea7fed
# ╠═e8944e65-f0e0-4169-bcb1-3d6aa2b30bf5
# ╠═6ee8254c-be59-46a2-8b9c-5e47b50237f7
# ╠═39456f2f-700d-4476-9c96-ad818c058363
# ╠═9a70e256-3fbe-4660-ae01-cd02d5136a89
# ╠═dd8cf5c2-ccb9-4e55-8b1a-7c5cde97c36a
# ╠═b2c8b3cb-7f6f-4f1e-9ccd-605bd110756a
# ╠═fc619c3a-556b-4bf7-bff9-28f9a9bc0fb2
# ╠═15bc9488-d29d-468e-8995-1ff5df9487e6
# ╠═0be83fa6-81c6-4d32-a2ff-db833750a3f1
# ╠═b7c1cde6-7d63-4497-8a3e-2df959ca6918
# ╠═9de59850-3bd2-46a9-a3f9-3b73055bdf1d
# ╠═7b07a2b8-eb69-4da7-ada5-4d44503e4274
# ╠═5d239b82-f341-446b-9aaf-64076385723f
# ╠═6a4bf9a3-fbb9-494a-b48e-000553cf5460
# ╠═4b17b7fd-3a5d-4560-bcaf-bffb1d096ab5
# ╠═aac0c842-14b6-48e3-bde1-30e32713046f
# ╠═27ffc2b1-9239-4fa4-b575-bab651196be6
# ╠═5fccfd98-c270-4cb5-b8e4-cce7e0e5f04b
# ╠═63c9378a-df62-41d1-b3a0-db5a1be3b3c5
# ╠═0902e42e-8435-4a45-ab7f-9d7969dce188
# ╠═9e6bd014-b57b-41a4-a170-9e466325d51e
# ╠═2f655bd9-ad59-4114-95b5-d4fdfb6447e6
# ╠═e8edab4f-c2cd-41dd-9aa8-2140f1641e56
# ╠═4ba9ce08-ef82-415a-bdcc-5c70b8d4a16e
# ╠═b50bbc59-dc26-4cec-b674-6d733f573421
# ╠═a7e2adf8-154b-4ad1-bc53-aed37c9846ad
# ╠═d679e948-8d30-4620-a4fb-f202c04117b8
# ╠═b018504c-bcfc-4c01-9b05-43936b43b615
# ╠═4f2abf90-3e70-422a-9a1f-a85083269a9a
# ╠═3304d896-71d8-4e81-9fb9-eb84a63d52c3
# ╠═ca2d64b3-4398-4f5d-a5c0-70d73999ae02
# ╠═2d082f91-5208-45ac-93e4-f45ce7d5fe66
# ╠═d7b60c09-d97f-4a78-9c1e-98dbf9da1372
# ╠═407f122d-01ae-4e6b-9ae0-b25b15cf23f8
# ╠═e7134856-8988-4ae9-bec2-53ab57124541
# ╠═b47c0f30-933c-47a4-8691-46c4d3d48ed7
# ╠═b6887cf2-dfd0-43f4-860f-31b18935caef
# ╠═a8b3097c-d965-435e-8712-2df777def9b6
# ╠═342e32ee-5343-4b4b-b008-9f2cdb44e9b6
# ╠═4efd0ae0-1322-4488-b80b-d7c5aa067d9d
# ╠═8cb42e3c-b5b9-4d0e-ae84-f7b9304920b0
# ╠═cdd5190e-1cb4-4cc6-96af-09d149a06f82
# ╠═c41e20bf-70b3-433d-add3-97ed29c98c84
# ╠═2ff0a042-f411-4994-b777-fa9bb2ce8c6b
# ╠═a0a42d56-1dc3-4701-80ee-ff80078e7734
# ╠═ccef9339-201e-4ea3-9eda-06972be2919a
# ╠═a8f4478d-37b8-4483-b8f1-d541f65f7b5d
# ╠═5f4819b9-3b5b-4700-8b1f-230e38a1cf28
# ╠═ed85fb2c-b35b-42fa-8933-b740ceb8f392
# ╠═5d6f15de-3c56-446b-8afa-9a688804ce4e
# ╠═092a1984-5515-4de6-b936-06344281d932
# ╠═f6bea570-3934-4c3b-8a73-d196e105479b
# ╠═d3937387-d7df-4228-8f3d-ff325a1b3d53
# ╠═76428346-cbff-4745-afd0-31a479ca2494
# ╠═ffe7d1a9-6c55-45a1-abd6-c4610d922dd0
# ╠═8c1bc363-0825-496a-bc35-d03d34a54220
# ╠═fad5c87a-2611-40ea-a6d7-a98975ddd7ed
# ╠═12ec6787-9a6b-4797-9296-fd8cbeb3c1b0
# ╠═e2342368-eaa0-41e1-be5f-0cfce2403c34
# ╠═91fc63c9-cab2-4151-be14-f856ffdeb201
# ╠═cb6f30e9-bd84-494a-856e-7dc284124a22
# ╠═66c7f41c-ba84-4522-9b9a-27b6cc3a6893
# ╠═d08e325e-5d2a-486c-87e5-cb0833e3882b
# ╠═53ce657b-1dd9-49c6-a152-fd3358d17b8f
# ╠═95f02f1c-914b-4363-a100-eb34512ee911
# ╠═df9ec094-a2ae-428d-aacb-f53a91c4eff8
# ╠═f1dd564a-1a34-4945-bba4-f142760ea533
# ╠═283b19e1-2b1b-4d7f-9e9a-3aaeb0eca290
# ╠═87d04f23-87cc-4c8d-a364-57e221ace9ce
# ╠═23dffd9a-4fd3-41ce-9645-943b6166d7f2
# ╠═71463a7f-35e0-4196-aa3c-0a0df5734696
# ╠═b59f861e-8ced-4900-9fdf-d59eef82e9e3
# ╠═5ef24d5a-e279-4ba0-83d0-4ac348df6290
# ╠═e426497b-55c1-4792-803f-012636c7559f
# ╠═20f36ba5-02a8-4883-b9ad-e9339b139b73
# ╠═84b3df87-33ae-4c0f-836d-a8b225e52a4e
# ╠═24ee0cfe-af10-4287-8265-4638ba20d042
# ╠═bd90358f-6190-4358-a25f-1a3a85a1a944
# ╠═30503dc1-e9e7-4876-9bbe-df8ddbb50c93
# ╠═5955727b-40f1-4eb4-b868-2ae7788c1625
# ╠═44fa6c8e-9e64-4262-9324-272a133e302e
# ╠═0a4a6191-c3dd-4551-b104-dd0595dd684e
# ╠═e7c9704d-e310-4b11-aa88-f41d6692635e
# ╠═1a956d44-14c6-41a9-8015-3691051d0cb6
# ╠═daa0b062-f7e4-4684-8a42-d85bc9043b71
# ╠═f4d510db-3b37-4f06-bc27-99447fd4fefe
# ╠═c8117fcf-3dfa-4064-8adb-509b17e0280b
# ╠═80a489d8-277d-412d-92ad-2c54dc98fbdc
# ╠═2639bd8b-d88c-42f4-b34c-fd206525d3c3
# ╠═b6d30924-90d8-4155-882e-8a72e0c7e797
# ╠═2321ad16-58d0-4332-9108-dfde9e6ec6cf
# ╠═3ddfebdc-7b41-4a17-8cd1-9169ad4ecf07
# ╠═16887603-f096-46fe-9fa1-601d1dc2605d
# ╠═ed06d8f1-a22e-4da2-a31c-50ae7176b31d
# ╠═e5d77450-c22b-4f63-8adb-0dc1854e3b8c
# ╠═05599d8f-82bd-4e3c-98e1-42473e7887ad
# ╠═9f5f1f2c-498f-460b-b6b4-a58b5c9f5804
# ╠═de4dee93-24e5-463b-96fd-078f1ac8f9c9
# ╠═65fceae5-cc43-4208-b6c5-ae77fdc238ea
# ╠═08e0976a-7652-4b8b-ad07-c0b8505e086c
# ╠═34f706ea-918d-43fc-a653-92d311b7df18
# ╠═01ca48c3-e201-4ad5-87fe-cbec83b7a009
# ╠═1cc1735a-7b11-4ecf-9b11-d0aac8a79c7a
# ╠═8a195f8f-c9ff-4f24-aca1-5b85d7421383
# ╠═ec3eb92f-9232-4309-8564-d8c9b620293b
# ╠═28c338dd-6a85-49ad-ac5b-e4df11210abd
# ╠═4592706f-7f90-41c4-bf11-4254fa6450a2
# ╠═c4c9f140-086a-483e-9bb8-5af8e452181e
# ╠═2fcdf706-5ff6-4e86-a6e9-839e9144655a
# ╠═67de2dbd-3633-4728-a6e4-53e2ddfb2625
# ╠═e9821459-f1b6-4884-8d0c-3017ac2d776d
# ╠═5c493ed5-e368-4550-88da-bf6c7338b879
# ╠═da45b3b3-c3c2-4fe8-a0ba-2530ec2a6e84
# ╠═89641cc8-72b2-4a4a-9fcb-e48bdf7a778b
# ╠═f5419815-56ef-4485-b72f-ee3bf2d3461b
# ╠═2a341070-a64e-41dc-b41b-c8e68730dc9f
# ╠═f33a007e-2e91-4206-bb80-65d1997247a3
# ╠═4ad15809-2bac-4c25-8139-66fb55311f9c
# ╠═85376195-8c5e-42a1-a307-a7f0721556b9
# ╠═d04b0d6a-96ba-4aba-b943-74cf46f8855c
# ╠═c99b3c6a-c29d-475d-bb79-ff23d505fe0a
# ╠═7b6dd0b4-6e9c-44bf-baf7-6dcc0709086d
# ╠═1938c38a-5dfa-43a2-8810-cc987a34cce4
# ╠═cdefdb31-f6fe-4f86-9f12-e01f068f1108
# ╠═aa0db854-c097-4df4-b139-058869a6be8f
# ╠═fca157a2-743b-452f-b035-44fe7363764b
# ╠═53a4439c-3714-4c95-a0db-9f4ca9d1105b
# ╠═3188f726-5bd3-46d7-bb7a-760ed1fcff05
# ╠═95c2de46-2634-4ab3-a7c8-d4251cdabd32
# ╠═51a64f33-e0e9-4365-bab8-6d134f27517d
# ╠═d154ace4-9a83-4b9a-8dce-3b2d1b569200
# ╠═f6b131d4-05eb-4005-b55c-5fad152f5090
# ╠═a8f062ed-659f-49d5-9cac-cdc47393e90a
# ╠═815e44dd-57cf-4104-95d9-1632ddc4ead4
# ╠═4a69919b-38e3-4a05-9299-8e6c50b5bd25
# ╠═80a26e90-77aa-4cb1-88de-856e8f366453
# ╠═765971bf-1ff6-453e-b90c-a00126c2095c
# ╠═1b4d58ea-ae61-4994-95c4-00fdfb2f4b75
# ╠═be87100b-a5c7-4050-b4ab-ad065b97091a
# ╠═0a7000be-88dd-4cba-a0b8-6f82ff1c9c07
# ╠═308ea3b1-8145-49a7-aa16-048666312620
# ╠═d76052a6-92a3-416e-89f7-f37e160b7d54
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
