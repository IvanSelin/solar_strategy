### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

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
end

# ╔═╡ 9ce6e2c1-69e4-47c6-830c-3c0eb13d784e
begin
	include("src//energy_draw.jl")
	include("src//time.jl")
	include("src//solar_radiation.jl")
	include("src//track.jl")
	include("src//utils.jl")
	include("src//strategy_calculation.jl")
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
track,segments = get_track_and_segments("data//data_australia.csv");

# ╔═╡ aef61af9-31df-4757-95bd-7fdfc0edbc1a
track

# ╔═╡ dbe59867-4be2-4a6a-9236-35dc96d2d387
plot(track.distance, track.altitude, title="Raw track data")

# ╔═╡ e53d252c-31b5-430c-a241-b41159243bc7
track_peaks, segments_peaks, points_peaks = keep_extremum_only_peaks_segments_with_points(track);

# ╔═╡ 164055ba-18bb-482d-816e-e278fc76a0d0
plot(track_peaks.distance, track_peaks.altitude, title="Peaks track data")

# ╔═╡ 1b6f993d-e4cc-4973-b5cf-269f01c7f409
md"# Parametrized merging"

# ╔═╡ 16405fa9-f44a-4b60-9b87-84bfd8aa8ed2
thresholds = 0:0.05:1

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
	income_plot = plot(track_points.distance, income, title="Power income");
	use_plot = plot(track_points.distance, use, title="Power use");
	time_plot = plot(track_points.distance, time, title="Time");
	energy_plot = plot(track_points.distance, energy, title="Energy");

	plot(income_plot, use_plot, time_plot, energy_plot, layout=(4,1), size=(1000,700), legend=false)
end

# ╔═╡ 0aeeb368-896e-4bc8-ad1c-d92b0c324392
plot_differences(test_income, test_use, test_time, test_energy, track, points[0.05])

# ╔═╡ cdf8794f-67ed-41df-bd9e-b1c57c82133b
plot(track[points[0.05],:].distance, [test_income, test_use, test_time, test_energy], title="Difference", labels=["income" "use" "time" "energy"])

# ╔═╡ b1a5d4ad-c65a-40c0-9710-18935f3f6d94
md"For some reason power use does really differs a lot!"

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

# ╔═╡ a41c1a7b-e8c5-4e4f-bf18-38e309381aaf
plot(diff(energy_diff_0_05))

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

# ╔═╡ 62e3cc8e-aa6b-45b0-b40c-3f4266d64c72
md"## Peaks"

# ╔═╡ 91ee1034-994a-49b9-a1ad-531ba8e48399
md"do it for peaks

make a comparison loop"

# ╔═╡ 41626c2c-3dbd-48e4-8b3b-9021aff1b7ec
md"make plot for energies point-by-point?"

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

# ╔═╡ c21f026f-ba8b-4e46-8f31-0b6d8b886472
last(p_energy) / start_energy

# ╔═╡ 8b71b075-0fbe-4803-8fe5-2f58ca12b550
md"Разница в 2% на peaks

А нормально ли это?

Почему разница всегда в + идёт по расходу?"

# ╔═╡ a53676bc-a9d5-4163-9d06-436cf2271e14
diff(p_use)

# ╔═╡ ea61effe-cf3c-4758-b683-5838f8698a8c
minimum(p_use)

# ╔═╡ 163eef34-8a2e-4057-9381-a798e9fdecf7
minimum(diff(p_use))

# ╔═╡ 72cc5d74-ac27-4639-b5e0-0960d60d8cd9
maximum(diff(p_use))

# ╔═╡ ee212e96-9af1-4f0f-8310-03bc84b27ce6
plot(diff(p_use))

# ╔═╡ 14d1d61e-620f-41d5-8d9f-90e2db816774
md"Редко в минус идёт, очень странно. С большими выбросами

Надо выяснять за счёт чего различия

И потом уравновесить метод, чтобы за дистанцию в ноль пришло хотя бы"

# ╔═╡ f3692e2b-2534-48a6-9f0f-be2bdd94ced1
diff_p_use = diff(p_use)

# ╔═╡ 40efb91d-70c5-4444-ac1d-cb1fcfc84447
md"буквально в 3 и 4 участке"

# ╔═╡ 3f5ea8a4-fd8c-4e6d-bcfe-4b8ee912a26e
length(diff_p_use)

# ╔═╡ 3b876e55-1b1f-43be-a18b-d958f8a88aa0
length(diff_p_use[diff_p_use .< 0.5])

# ╔═╡ 893338ee-720d-4020-be54-38fd0b50dee0
diff_p_use[diff_p_use .< 0.5]

# ╔═╡ 5a7d83ce-aed8-4d6a-b0bc-dab5a876d9dc
diff_p_use[diff_p_use .>= 0.5]

# ╔═╡ 89156cbf-8d27-4b3a-a7fa-706c94a9b772
sum(diff_p_use[diff_p_use .>= 0.5])

# ╔═╡ 9189d3f6-e9cf-4b12-839f-e2ac357fa098
length(diff_p_use[diff_p_use .< 0.25])

# ╔═╡ 9de08059-a2fb-4431-a255-523906063ae9
sum(diff_p_use[diff_p_use .>= 0.25])

# ╔═╡ f6082a6b-61ff-4f01-8a22-4143bbabb0fc
sum(diff_p_use[diff_p_use .>= 0.15])

# ╔═╡ 6e1c6655-645d-442b-bea3-566898e14ac2
sum(diff_p_use[diff_p_use .>= 0.1])

# ╔═╡ 9dd448a0-c8b0-4f3b-8dfc-e3ed1f8a4bf5
sum(diff_p_use[diff_p_use .>= 0.05])

# ╔═╡ 1c5dcf72-afa8-4d83-87da-7efd372e71b6
sum(diff_p_use[diff_p_use .>= 0.025])

# ╔═╡ e91c0b2b-727f-4391-ae77-74d3ae342f9d
Plots.histogram(diff_p_use)

# ╔═╡ fe5fcc31-ee4b-4425-bf5a-12f9c2e2a0b9
md"# Понимаем что не так"

# ╔═╡ f0c74358-79b9-41b4-8cd0-672654e62b3b
md"Сперва смотрим пиковые данные, где расхождение

А расхождение буквально в 3, 4 и 6 участке"

# ╔═╡ a0309184-5663-49c0-8285-485e16b64e5b
diff_p_use

# ╔═╡ 8e03b7ed-24af-48e9-a545-3fe8e3088f91
track_peaks[1:9,:]

# ╔═╡ 035fca65-b9ee-4b9d-9279-532dbb54c882
segments_peaks[1:8,:]

# ╔═╡ facff3dd-18ff-4b0d-b406-d15ca3204d0a
plot(track_peaks.distance[1:9], track_peaks.altitude[1:9])

# ╔═╡ 62aa9ae3-b563-43c5-885a-4bec151a002b
md"Исходная трасса"

# ╔═╡ aee37fb7-1408-499d-baaa-37cf329cb4db
track[track.distance .< 1065,:]

# ╔═╡ d94db435-491a-4799-aadd-95baf06074ad
segments[1:81,:]

# ╔═╡ 39a775fc-dc4b-4664-9fd1-dd008c45a7b8
plot(track.distance[1:82,:], track.altitude[1:82,:])

# ╔═╡ 5adf37eb-4815-489a-8602-9a14ab016593
md"Рассмотрим поплотнее участки 3, 4 и 6 из peaks"

# ╔═╡ e895ceae-d260-44df-9a40-cf46a47dd789
plot(track.distance[1:82,:], track.altitude[1:82,:], markers=:diamond)

# ╔═╡ e5769e28-8dd8-4f58-93b6-58c75c270871
md"Выяснилось что max и min значения высоты нормально учитываются, но одинаковые значения подряд (плато) не учитываются правильно

Написал новую функцию расчёта, сейчас будем пробовать"

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

# ╔═╡ 43ce04fb-8443-493c-8f9f-7abc19f01cd7
diff(pl_use)

# ╔═╡ ee1d6a44-941a-4ac5-9928-b6790329334e
md"Что-то тут не так.

Сравним с обычнми peaks"

# ╔═╡ cdf74c41-75ac-4d1d-be16-a6fff7bc2165
orig_peaks_income, orig_peaks_use, orig_peaks_time = simulate_run_income_use_time(opt_speed, segments_peaks, start_datetime)

# ╔═╡ 4fa48c6d-d4d4-4e4e-950a-f4198b046297
new_peaks_income, new_peaks_use, new_peaks_time = simulate_run_income_use_time(opt_speed, segments_peaks_pl, start_datetime)

# ╔═╡ e0bb15e1-0bca-4d10-9fb6-08e54074152d
income, use, time = simulate_run_income_use_time(opt_speed, segments, start_datetime)

# ╔═╡ aad628de-6530-4cce-8ea8-6b258716b75e
plot(track_peaks.distance, segment_data_to_track_data(orig_peaks_use, 0.), title="Исходные и новые peaks траты энергии")

# ╔═╡ 6fffcf11-0236-4e90-922b-34c6ae8796e6
plot!(track_peaks_pl.distance, segment_data_to_track_data(new_peaks_use, 0.), title="Исходные и новые peaks траты энергии")

# ╔═╡ 3e229834-bc84-4eee-8ee5-b4c4955a9553
md"на графике не видно, будем искать по-другому"

# ╔═╡ 3d9560f0-3843-4f03-9b3f-2485c7605cc6
last(new_peaks_use)

# ╔═╡ a2086714-5522-4d12-8fd7-34cb40fbbad8
last(orig_peaks_use)

# ╔═╡ 84cb6d65-fbc8-4470-a4e7-3116532cba5b
last(orig_peaks_use) - last(new_peaks_use)

# ╔═╡ 8c711de5-2bd7-492d-8f20-a68538e3f51e
md"С новыми peaks ситуация стала даже немного хуже!"

# ╔═╡ e11b95a0-a96d-4f83-a61a-f9e8698368df
md"Сравним с исходными данными"

# ╔═╡ 4d887c9a-02dc-471d-9299-f6b35ed77daf
last(use)

# ╔═╡ 21c94bfb-ee4b-40c7-ab73-40681c3f77d6
last(use) - last(new_peaks_use)

# ╔═╡ 58f2a9ab-a24b-408d-80ab-8a9ca688d572
last(use) - last(orig_peaks_use)

# ╔═╡ 9ceb2c87-80e4-4dce-8cd9-7570567d4ecb
md"Кризис пройден, стало лучше)

Но расхождение всё равно есть

И использование peaks приводит к СНИЖЕННЫМ оценкам по энергии"

# ╔═╡ 7ca5507d-6083-41ae-a127-c23062f86766
md"Но всё равно надо понять где и что идёт не так даже с новыми peaks.

Опять есть какой-то скачок в начале, который надо исследовать"

# ╔═╡ aed6dedd-5c3e-48df-aab2-c4f04b4209d5
md"## Глубокое сравнение"

# ╔═╡ 6aaa20a3-5d21-4819-b9be-b185d345a0dd
md"### Подготовка данных"

# ╔═╡ 8d1adddd-0c72-4831-9037-b9c7f9a7aa4c
track_res = deepcopy(track);

# ╔═╡ 3f4b7ad1-ee6a-4702-bc27-303e8e254e56
track_res.index = 1:size(track,1)

# ╔═╡ da6ac0a8-1d55-497d-95ab-bf8300c93c8b
track_res.use = segment_data_to_track_data(use, 0.)

# ╔═╡ d3b55c8c-0261-4a15-9667-08ed2ba43f91
track_res.income = segment_data_to_track_data(income, 0.)

# ╔═╡ 52f9b92d-f06d-410c-94b2-7475e20be328
track_res.time = segment_data_to_track_data(time, 0.)

# ╔═╡ 2f3e57a6-33e5-45b7-8998-52d6374a11d0
track_res

# ╔═╡ f167adc7-83ed-4e7c-9e93-374d0177f2fb
track_peaks_new_res = deepcopy(track_peaks_pl);

# ╔═╡ c75b72df-a899-4912-a7e5-c8900d75711c
track_peaks_new_res.use = segment_data_to_track_data(new_peaks_use, 0.)

# ╔═╡ a899b3e5-935b-43e2-8dc1-207b7bb9c657
track_peaks_new_res.income = segment_data_to_track_data(new_peaks_income, 0.)

# ╔═╡ 49d190ad-31a5-4ce9-9978-bf510e7890a1
track_peaks_new_res.time = segment_data_to_track_data(new_peaks_time, 0.)

# ╔═╡ 9f6fd121-d96f-4745-9283-f2e6b039b90a
track_peaks_new_res

# ╔═╡ ca505019-fb44-4ba2-a36d-457a5e9eccc4
track_res_points = track_res[track_peaks_new_res.index,:]

# ╔═╡ ece8624e-824c-4725-8d62-92e91cbc5e44
md"Большая разница на участке 8-14"

# ╔═╡ 9a568b18-bf33-42d4-b239-103417ce4075
Plots.histogram(diff(track_res_points.use)-diff(track_peaks_new_res.use))

# ╔═╡ 6859857f-fdd6-42b7-bde1-fce03961afb5
md"Много маленьких отличий


Смотрим подробней на участки"

# ╔═╡ 4770ab6c-da5d-4e23-b99c-6c4540579ecd
track_res[8:14,:]

# ╔═╡ 57bfd811-eec2-4840-8f7b-d9dd0cf46359
segments[8:13,:]

# ╔═╡ 8747d2fa-3630-4301-8a37-9ed7182e5b10
segments_peaks_pl[1:5,:]

# ╔═╡ dd1803c3-b9df-400e-b796-eaf4fc96d499
sum(segments.diff_distance[8:13])

# ╔═╡ 0316cd3a-202c-447f-8396-0fc78abf97d4
sum(segments.diff_altitude[8:13])

# ╔═╡ 86ce26a7-5d56-46c6-bb39-d25e6f410b31
sum(segments.altitude[8:13]) / 6.

# ╔═╡ 1801c7c6-09fc-4b39-ba96-fca8be1fec28
md"Немного неправильно считается средняя высота

Но это по идее не сильно влияет на расчёт"

# ╔═╡ 5ba67e15-0aa6-4bc8-b271-faef5dae0fe3
md"длина совпадает

Надо теперь пробовать считать энергию для обоих частей"

# ╔═╡ 18cc16be-b957-4ad7-9f1f-f2a9bce8cb2e
md"### Считаем для peaks"

# ╔═╡ b0c8fafd-8b1f-4520-a05b-ba9a1424ad5f
md"Сперва сравним, совпадает ли с ручным расчётом"

# ╔═╡ b9d3ae44-95d0-4151-9d64-fb3f9f528577
# for peaks
mech_test_peaks = mechanical_power_calculation_alloc(opt_speed / 3.6, segments_peaks_pl.slope[3], segments_peaks_pl.diff_distance[3]) / 3600.

# ╔═╡ b249b36a-d29c-416b-848d-293246ae5bfe
mech_peaks_fact = track_peaks_new_res.use[4]-track_peaks_new_res.use[3]

# ╔═╡ 88b6e352-fac5-4ba4-953b-cc1eb56b8c0c
mech_test_peaks - mech_peaks_fact

# ╔═╡ a68baf40-b94f-4f98-9bec-4f9c81d255fd
md"Совпадает

Теперь надо посчитать для обычной трассы" 

# ╔═╡ 3ca7fc23-b918-4186-b6e2-7caa610a2a5f
md"### Считаем для обычной трассы"

# ╔═╡ 8bd276cb-e8ef-4279-80c5-fafb73d389cf
md"Считаем по сегментам с точками с 8 по 14, то есть с 8 по 13 сегмент (т.к. он с 13 по 14 точку)"

# ╔═╡ 81b78f75-7bd7-4699-8d61-0d1674aec5cd
segments[8:13,:]

# ╔═╡ e4f13a18-c617-4f2c-9e17-e692cd9b74d2
md"Скармилваем эти сегменты в расчёт. Должно быть 6 сегментов (14-8=6)"

# ╔═╡ 967129f5-9f09-4974-8bc5-57a53aca2aff
mech_test = mechanical_power_calculation_alloc.(opt_speed / 3.6, segments.slope[8:13], segments.diff_distance[8:13]) / 3600.

# ╔═╡ a6e9f63d-58f5-47b5-adb4-d489996bc2af
md"Сумма по этим участкам должна быть такой же, как и в peaks трассе. Ну хотя бы примерно"

# ╔═╡ 9c8e5ce4-dfc9-4b24-8d01-c20849221e37
sum(mech_test)

# ╔═╡ f0fa67df-0d75-40e8-beca-31b5e330de7c
mech_test_peaks - sum(mech_test)

# ╔═╡ 4ad1783a-afbc-4ab7-b043-4931569c330b
md"Похоже на правду? разница в 9e-5, т.е. всё хорошо

Но откуда берётся ошибка в итоге тогда?"

# ╔═╡ 150106d9-ce77-4729-ad9c-20d48c6361b9
md"Это мы сравнили ручные расчёты.

А теперь надо сравнить с тем, что нам посчиталось по факту. Может ТАМ что-то не то? "

# ╔═╡ bdc1c037-4295-4e09-95e7-a1837b7317da
md"### Расчёт по факту для обычной трассы"

# ╔═╡ b70dda56-903f-47c0-bd42-c6ecca40c8d2
md"Смотрим результаты

В результатах сегменты мапятся на трэк. Сегментов на 1 меньше, чем точек.

Поэтому, например, в первом ряду будет 0 по энрегии, а во втором, сколько понадобилось от 1 к 2 точке."

# ╔═╡ 016af154-833c-4024-96e1-97ba3de9e487
track_res[7:14,:]

# ╔═╡ 789c8b60-06d8-45f7-9131-8710d51f5c7f


# ╔═╡ 1250bdfd-9111-494b-b823-208cdf0eb172
track_res[9:14,:]

# ╔═╡ 5f0169be-4363-42c3-9e12-315edcf37a7b
md"Но здесь результаты указаны с учётом накопления энергии.

А нам надо взять результаты без накопления. Для этого надо сделать diff.

Diff производит массив размера на 1 меньше. Поэтому возмём результаты с 8 по 14-й и от них diff"

# ╔═╡ b8f86151-a304-4d88-843d-85ba86908a72
mech_fact = diff(track_res.use[8:14])

# ╔═╡ 263b7313-2246-42bd-84ae-a032fe6e2a85
md"Смотрим сколько получилось"

# ╔═╡ be2ccdb1-b310-46d4-b319-ce51b45399a0
sum(mech_fact)

# ╔═╡ 0795f1c3-832c-4b38-a029-960b9f7099e5
md"И в сравнении с ранее посчитанным"

# ╔═╡ 4ea9933e-cef6-4638-98c9-7a77309916ac
sum(mech_fact) - sum(mech_test)

# ╔═╡ 404bd6e3-a8b0-4542-8ebe-810ef26dcb47
mech_fact - mech_test

# ╔═╡ ccc7eb5e-6275-4135-ba44-8fd6862f49e4
md"Опа-па! почему-то есть расхождения, и причём серьёзные. Где-то что-то забыто, наверное

Надо разбираться"

# ╔═╡ d4ddfacd-d430-4a7d-b503-9b9f9a260cba
md"Почему-то не совпадают ручные с не-ручными"

# ╔═╡ dcc4bcde-ca0e-490b-a56f-96997f72e667
md"### Разбираемся с расчётом для обычной трассы"

# ╔═╡ 60db43d2-9e98-4e11-a3e9-12cd653190c9
segments[8:13,:]

# ╔═╡ 72f6f45f-2e89-4a4b-872a-d47f932effd4
md"Проблемные сегменты - с 10 по 11, с 11 по 12, с 12 по 13"

# ╔═╡ a6ef9143-0fcf-4af2-a80e-83b79ac50cd2
md"Ещё раз как оно должно считаться"

# ╔═╡ b1ddcaf2-fa92-4eda-a509-0c3792316572
mech_test

# ╔═╡ 00f748d3-ad88-42c5-9b2c-6b443bd655c1
mechanical_power_calculation_alloc.(opt_speed / 3.6, segments.slope[8:13], segments.diff_distance[8:13]) / 3600.

# ╔═╡ fa896da9-ff1f-4e62-8550-84bf958cd5f9
mechanical_power_calculation_alloc.(fill(opt_speed / 3.6, 6), segments.slope[8:13], segments.diff_distance[8:13]) / 3600.

# ╔═╡ e69786b0-8641-42b4-bd64-320787dff8a1
md"### Гипотеза про segment data to track data" 

# ╔═╡ dbd57c19-c91c-4407-b798-c9b1879a3c10
md"Может проблема в segment data to track data?"

# ╔═╡ 7ba1a20e-4a3c-4e3e-be3e-9b90ebc54d37
use

# ╔═╡ 52d2c069-6a8d-4bf4-be2b-223f6eb9e72e
segment_data_to_track_data(use, 0.)

# ╔═╡ 8a079fc0-6bee-48ef-93cc-27804134b935
md"### Не учитываются электрические потери?"

# ╔═╡ 1b9f27bc-934f-48c4-86b1-37f541cb8a0c
md"не учитываются электрические потери? которые electrical power calculation"

# ╔═╡ de02cc62-f6f4-450b-a002-c5349c7f8e01
md"для peaks"

# ╔═╡ 17cbf51c-645f-49fc-be63-ae60ce86923a
mech_test_peaks_el_loss = mechanical_power_calculation_alloc(opt_speed / 3.6, segments_peaks_pl.slope[3], segments_peaks_pl.diff_distance[3]) / 3600. + electrical_power_calculation(segments_peaks_pl.diff_distance[3], opt_speed / 3.6) / 3600.

# ╔═╡ e017637e-4885-4f32-92be-6b148c3e603b
mech_test_peaks - mech_test_peaks_el_loss

# ╔═╡ bb18bce4-140b-459d-b7c3-1bf29621ae37
md"Разницы почти нет, но стоит попробовать дальше"

# ╔═╡ 88f5c027-4568-4c80-83ee-b707d6f72e0f
md"для обычной трассы"

# ╔═╡ 61c98942-8299-4dea-bafb-b728fd567f28
mech_test_el_loss = mechanical_power_calculation_alloc.(opt_speed / 3.6, segments.slope[8:13], segments.diff_distance[8:13]) / 3600. + electrical_power_calculation.(segments.diff_distance[8:13], opt_speed / 3.6) / 3600.

# ╔═╡ 9b9746c3-096a-437a-be05-fa4daba0509a
mech_test_el_loss - mech_test

# ╔═╡ 938618cc-cf5b-45ff-8ad6-69a6dc0a99f0
mech_test_el_loss - mech_fact

# ╔═╡ 6776c519-95f8-4d7f-975c-6a78b8e9e6c6
md"Вот и нашли расхождение в сравнении!

Теперь надо понять, фигурирует ли оно в сравнении peaks и не peaks" 

# ╔═╡ 8b03a23f-51a6-4a74-b6aa-1cf74c5a6677
md"У нас есть несколько расчётов:

peaks c потерями (mech_test_peaks_el_loss) и без (mech_test_peaks) руками, peaks методом (mech_peaks_fact)

обычный с потерями (mech_test_el_loss) и без (mech_test) руками, обычный методом (mech_fact)

что с чем сравнивать?"

# ╔═╡ ceaae763-3a4f-401a-a980-ce0678f11d37
md"проверяем расчёты руками без потерь, peaks и обычный"

# ╔═╡ e6b48621-eb55-4f93-80fe-8bde497512b3
sum(mech_test) - mech_test_peaks

# ╔═╡ f434f5c1-16f5-47e0-b772-74da18fefbdb
md"отлично, расхождения нет

теперь надо сравнивать потери"

# ╔═╡ d094bb28-ce6a-496f-a1b9-17e0db817085
peaks_loss = electrical_power_calculation(segments_peaks_pl.diff_distance[3], opt_speed / 3.6) / 3600.

# ╔═╡ 3a262d65-bd0d-470c-a738-f959e5e89143
track_loss = electrical_power_calculation.(segments.diff_distance[8:13], opt_speed / 3.6) / 3600.

# ╔═╡ 3c253f05-eb91-4e6f-a633-04912c36bb52
sum(track_loss) - peaks_loss

# ╔═╡ bedc9813-4204-4a9b-b784-55ea6fed5364
md"вот она разница!"

# ╔═╡ 6f83a0a8-f919-4830-90f9-5320914cba3c
md"смотрим, такая ли же в целом разница, если сравнивать с итоговыми результатами (фактическими)"

# ╔═╡ 81ccfc68-6dbc-402c-af18-871acf4022a5
sum(mech_fact) - mech_peaks_fact

# ╔═╡ eba5d743-da8d-4e2c-a908-961dacf24264
(sum(mech_fact) - mech_peaks_fact) - (sum(track_loss) - peaks_loss)

# ╔═╡ bf6b6f23-0878-4a52-bb13-b2faa5a2fa3b
md"Такая же, дело именно в электрических потерях!

Надо происследовать их!"

# ╔═╡ d85603c0-1e4d-4a92-ba00-e1a602c533da
# track
electrical_power_calculation.(segments.diff_distance[8:13], opt_speed / 3.6) / 3600.

# ╔═╡ 32356218-ecc6-4f2c-90a2-67244dcca04e
# track not broadcasted
electrical_power_calculation(segments.diff_distance[8:13], opt_speed / 3.6) / 3600.

# ╔═╡ 7af934e0-7c7c-497a-852d-343200b6e5a0
sum(electrical_power_calculation.(segments.diff_distance[8:13], opt_speed / 3.6) / 3600.)

# ╔═╡ 4194b150-1874-4ae2-b153-e3761cb724e1
segments.diff_distance[8:13]

# ╔═╡ 80342fe0-397b-4681-9178-8ab2bde90323
# peaks
electrical_power_calculation(segments_peaks_pl.diff_distance[3], opt_speed / 3.6) / 3600.

# ╔═╡ af2ed6ef-17ee-4613-9501-05f331bea410
segments_peaks_pl.diff_distance[3]

# ╔═╡ 31909481-bc95-49b3-8081-071cd888814f
function electrical_power_calculation_as_is(speed_ms, diff_distance)
    power_onboard = 40; # Wt, 0.04kWt
    return power_onboard .* diff_distance ./ speed_ms;
end

# ╔═╡ f71f9b5e-9654-49ac-803b-169652cc78e9
electrical_power_calculation_as_is(segments.diff_distance[8:13], opt_speed / 3.6) / 3600.

# ╔═╡ 94145064-52b0-4688-8947-6a1b3b34c553
function electrical_power_calculation_modified(speed_ms, diff_distance)
    power_onboard = 40; # Wt, 0.04kWt
    return power_onboard * diff_distance / speed_ms;
end

# ╔═╡ b027d115-648f-4c53-9448-0d1fbdd854fa
electrical_power_calculation_modified.(segments.diff_distance[8:13], opt_speed / 3.6) / 3600.

# ╔═╡ 304730fb-f35f-4cb6-969d-1349215ce8d7
opt_speed / 3.6

# ╔═╡ cf3fa8c1-3734-4a78-bb4b-a1d4f14c04e9
electrical_power_calculation_modified(segments_peaks_pl.diff_distance[3], opt_speed / 3.6) / 3600.

# ╔═╡ 04b5272f-3c0b-417a-bd3a-15b2716b94ec
electrical_power_calculation_modified.(segments_peaks_pl.diff_distance[3], opt_speed / 3.6) / 3600.

# ╔═╡ 886a3755-105b-42d3-97c1-2e3f4fdf8541
md"скорее всего напутал с единицами измерения

сейчас это Вт * м / (м/с) = Вт * с"

# ╔═╡ aa84cb6b-f7d1-428b-be48-bb0752fdc71e
40 * segments.diff_distance[8:13] / (opt_speed / 3.6) / 3600.

# ╔═╡ 7d72b8eb-5f99-43c9-b540-5231f9a7d40f
40. * segments_peaks_pl.diff_distance[3] / (opt_speed / 3.6) / 3600.

# ╔═╡ da48268b-5875-4d70-bb11-c57f0680218e
md"Был перепутан порядок аргументов. В функции было (speed_ms, diff_distance), а вызывалось как (diff_distance, speed_ms)"

# ╔═╡ 06a27299-6c51-42a0-9c31-98587898d895


# ╔═╡ f72e8f8f-b532-497e-892e-6708c83edab5
md"Здесь будет результат"

# ╔═╡ c3392f6c-5883-40e7-80e4-f6d752420837
md"### Отдельный кусок трассы" 

# ╔═╡ ca843464-5ba2-4d1c-ab64-5214b01fe305
md"Попробуем взять только нужный кусок трассы и прогнать симуляцию там"

# ╔═╡ 945c160f-ff97-428d-8a98-b8d59cb6afe9
track_short = deepcopy(track[8:14,:]);

# ╔═╡ 08b9835f-5600-409a-8d67-3fc599efbc03
track_short.index = 8:14;

# ╔═╡ de4ba031-9cc2-47b6-9cc9-45850ca3befd
track_short

# ╔═╡ e05889f7-d92a-4e11-b0fa-7dce89ad7f96
start_datetime

# ╔═╡ f138441b-410b-4ef3-85cf-3e052ce90351
start_datetime_short = DateTime(2023,1,1,10,0,42)

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
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Metrics = "cb9f3049-315b-4f05-b90c-a8adaec4da78"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlotlyBase = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
CSV = "~0.10.11"
DataFrames = "~1.6.1"
Metrics = "~0.1.2"
Optim = "~1.7.7"
Peaks = "~0.4.4"
PlotlyBase = "~0.8.19"
Plots = "~1.38.17"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "4fb2ef3c56a9d0920e9416faca24454229e043eb"

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
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
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
# ╠═cff82de0-3d15-11ee-37ff-035acc6e8e4a
# ╠═9ce6e2c1-69e4-47c6-830c-3c0eb13d784e
# ╠═7c89a14e-f6c7-4a09-8825-a81d1ba9b824
# ╠═39d80a54-3675-47ff-82bc-deb1daf03a59
# ╠═1b605aae-c293-49b0-a06d-fcff2ef3c12a
# ╠═d6fb7d82-6bd9-4da9-b183-67fb16ab9f66
# ╠═aef61af9-31df-4757-95bd-7fdfc0edbc1a
# ╠═dbe59867-4be2-4a6a-9236-35dc96d2d387
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
# ╠═0aeeb368-896e-4bc8-ad1c-d92b0c324392
# ╠═cdf8794f-67ed-41df-bd9e-b1c57c82133b
# ╠═b1a5d4ad-c65a-40c0-9710-18935f3f6d94
# ╠═eb2ddfbc-a8ab-4b88-ac46-37af7c563d42
# ╠═501fc5bb-8f50-4a19-b2a0-3fafe9a7738e
# ╠═322e9eb5-7da5-492f-a44d-2f7ee8c2ffcb
# ╠═fcbeb07f-3c3f-4973-ac1a-f133c7bdd4f3
# ╠═a41c1a7b-e8c5-4e4f-bf18-38e309381aaf
# ╠═3cd8689b-5f6c-481b-9e5b-ec0356597598
# ╠═35102b6b-1a3f-4c99-b03a-86fb5ad00e14
# ╠═6950f324-b23c-4b49-b391-586bda1b872e
# ╠═ce77f956-bd0a-48c1-a2b2-1f26a49877d8
# ╠═425a892c-b550-486e-b376-8566ae086f6a
# ╠═d9a12d12-8c19-4cd3-88ba-847e1c273850
# ╠═f84953c2-bed0-4375-a9d1-0bb646c832e9
# ╠═c26613ea-827d-46b8-9bd2-15ddf6cd885b
# ╠═b845f4f8-22c0-49cb-9580-d51c19e115bd
# ╠═62e3cc8e-aa6b-45b0-b40c-3f4266d64c72
# ╠═91ee1034-994a-49b9-a1ad-531ba8e48399
# ╠═41626c2c-3dbd-48e4-8b3b-9021aff1b7ec
# ╠═baa774fc-d955-4b69-8a32-7bf7b6a92997
# ╠═da13556d-62b0-49ce-867f-f64b8f227599
# ╠═ecc4c5e7-7b9f-47b2-9f35-85f6c241f4ea
# ╠═c21f026f-ba8b-4e46-8f31-0b6d8b886472
# ╠═8b71b075-0fbe-4803-8fe5-2f58ca12b550
# ╠═a53676bc-a9d5-4163-9d06-436cf2271e14
# ╠═ea61effe-cf3c-4758-b683-5838f8698a8c
# ╠═163eef34-8a2e-4057-9381-a798e9fdecf7
# ╠═72cc5d74-ac27-4639-b5e0-0960d60d8cd9
# ╠═ee212e96-9af1-4f0f-8310-03bc84b27ce6
# ╠═14d1d61e-620f-41d5-8d9f-90e2db816774
# ╠═f3692e2b-2534-48a6-9f0f-be2bdd94ced1
# ╠═40efb91d-70c5-4444-ac1d-cb1fcfc84447
# ╠═3f5ea8a4-fd8c-4e6d-bcfe-4b8ee912a26e
# ╠═3b876e55-1b1f-43be-a18b-d958f8a88aa0
# ╠═893338ee-720d-4020-be54-38fd0b50dee0
# ╠═5a7d83ce-aed8-4d6a-b0bc-dab5a876d9dc
# ╠═89156cbf-8d27-4b3a-a7fa-706c94a9b772
# ╠═9189d3f6-e9cf-4b12-839f-e2ac357fa098
# ╠═9de08059-a2fb-4431-a255-523906063ae9
# ╠═f6082a6b-61ff-4f01-8a22-4143bbabb0fc
# ╠═6e1c6655-645d-442b-bea3-566898e14ac2
# ╠═9dd448a0-c8b0-4f3b-8dfc-e3ed1f8a4bf5
# ╠═1c5dcf72-afa8-4d83-87da-7efd372e71b6
# ╠═e91c0b2b-727f-4391-ae77-74d3ae342f9d
# ╠═fe5fcc31-ee4b-4425-bf5a-12f9c2e2a0b9
# ╠═f0c74358-79b9-41b4-8cd0-672654e62b3b
# ╠═a0309184-5663-49c0-8285-485e16b64e5b
# ╠═8e03b7ed-24af-48e9-a545-3fe8e3088f91
# ╠═035fca65-b9ee-4b9d-9279-532dbb54c882
# ╠═facff3dd-18ff-4b0d-b406-d15ca3204d0a
# ╠═62aa9ae3-b563-43c5-885a-4bec151a002b
# ╠═aee37fb7-1408-499d-baaa-37cf329cb4db
# ╠═d94db435-491a-4799-aadd-95baf06074ad
# ╠═39a775fc-dc4b-4664-9fd1-dd008c45a7b8
# ╠═5adf37eb-4815-489a-8602-9a14ab016593
# ╠═e895ceae-d260-44df-9a40-cf46a47dd789
# ╠═e5769e28-8dd8-4f58-93b6-58c75c270871
# ╠═691e1f9a-5d10-475f-9ede-8d687f9d07b1
# ╠═f17bfc9f-01b5-4535-b05e-a4251e6a3594
# ╠═1c54eea3-b89a-4c7c-a7ab-7759f1bd4242
# ╠═89fdba6b-0c3f-4dd2-ae5b-a26dc1a90368
# ╠═0f248c56-5ab5-4183-9260-634be0fc9736
# ╠═3e84943e-1f80-44d6-b3fc-9f7c0553c4c1
# ╠═43ce04fb-8443-493c-8f9f-7abc19f01cd7
# ╠═ee1d6a44-941a-4ac5-9928-b6790329334e
# ╠═cdf74c41-75ac-4d1d-be16-a6fff7bc2165
# ╠═4fa48c6d-d4d4-4e4e-950a-f4198b046297
# ╠═e0bb15e1-0bca-4d10-9fb6-08e54074152d
# ╠═aad628de-6530-4cce-8ea8-6b258716b75e
# ╠═6fffcf11-0236-4e90-922b-34c6ae8796e6
# ╠═3e229834-bc84-4eee-8ee5-b4c4955a9553
# ╠═3d9560f0-3843-4f03-9b3f-2485c7605cc6
# ╠═a2086714-5522-4d12-8fd7-34cb40fbbad8
# ╠═84cb6d65-fbc8-4470-a4e7-3116532cba5b
# ╠═8c711de5-2bd7-492d-8f20-a68538e3f51e
# ╠═6d0d758e-1df3-4dec-b6cf-fb6cc453a6c1
# ╠═3a73d079-aa41-4b61-90b8-a93b94d8159c
# ╠═4c1c1493-f48a-49ff-8e49-bd850819d7f8
# ╠═7c68e07a-29bf-48d2-b6c4-595df52abdf7
# ╠═f85e9a18-df29-449a-bc7e-34a3024d03e0
# ╠═e11b95a0-a96d-4f83-a61a-f9e8698368df
# ╠═4d887c9a-02dc-471d-9299-f6b35ed77daf
# ╠═21c94bfb-ee4b-40c7-ab73-40681c3f77d6
# ╠═58f2a9ab-a24b-408d-80ab-8a9ca688d572
# ╠═9ceb2c87-80e4-4dce-8cd9-7570567d4ecb
# ╠═7ca5507d-6083-41ae-a127-c23062f86766
# ╠═aed6dedd-5c3e-48df-aab2-c4f04b4209d5
# ╠═6aaa20a3-5d21-4819-b9be-b185d345a0dd
# ╠═8d1adddd-0c72-4831-9037-b9c7f9a7aa4c
# ╠═3f4b7ad1-ee6a-4702-bc27-303e8e254e56
# ╠═da6ac0a8-1d55-497d-95ab-bf8300c93c8b
# ╠═d3b55c8c-0261-4a15-9667-08ed2ba43f91
# ╠═52f9b92d-f06d-410c-94b2-7475e20be328
# ╠═2f3e57a6-33e5-45b7-8998-52d6374a11d0
# ╠═f167adc7-83ed-4e7c-9e93-374d0177f2fb
# ╠═c75b72df-a899-4912-a7e5-c8900d75711c
# ╠═a899b3e5-935b-43e2-8dc1-207b7bb9c657
# ╠═49d190ad-31a5-4ce9-9978-bf510e7890a1
# ╠═9f6fd121-d96f-4745-9283-f2e6b039b90a
# ╠═ca505019-fb44-4ba2-a36d-457a5e9eccc4
# ╠═ece8624e-824c-4725-8d62-92e91cbc5e44
# ╠═9a568b18-bf33-42d4-b239-103417ce4075
# ╠═6859857f-fdd6-42b7-bde1-fce03961afb5
# ╠═4770ab6c-da5d-4e23-b99c-6c4540579ecd
# ╠═57bfd811-eec2-4840-8f7b-d9dd0cf46359
# ╠═8747d2fa-3630-4301-8a37-9ed7182e5b10
# ╠═dd1803c3-b9df-400e-b796-eaf4fc96d499
# ╠═0316cd3a-202c-447f-8396-0fc78abf97d4
# ╠═86ce26a7-5d56-46c6-bb39-d25e6f410b31
# ╠═1801c7c6-09fc-4b39-ba96-fca8be1fec28
# ╠═5ba67e15-0aa6-4bc8-b271-faef5dae0fe3
# ╠═18cc16be-b957-4ad7-9f1f-f2a9bce8cb2e
# ╠═b0c8fafd-8b1f-4520-a05b-ba9a1424ad5f
# ╠═b9d3ae44-95d0-4151-9d64-fb3f9f528577
# ╠═b249b36a-d29c-416b-848d-293246ae5bfe
# ╠═88b6e352-fac5-4ba4-953b-cc1eb56b8c0c
# ╠═a68baf40-b94f-4f98-9bec-4f9c81d255fd
# ╠═3ca7fc23-b918-4186-b6e2-7caa610a2a5f
# ╠═8bd276cb-e8ef-4279-80c5-fafb73d389cf
# ╠═81b78f75-7bd7-4699-8d61-0d1674aec5cd
# ╠═e4f13a18-c617-4f2c-9e17-e692cd9b74d2
# ╠═967129f5-9f09-4974-8bc5-57a53aca2aff
# ╠═a6e9f63d-58f5-47b5-adb4-d489996bc2af
# ╠═9c8e5ce4-dfc9-4b24-8d01-c20849221e37
# ╠═f0fa67df-0d75-40e8-beca-31b5e330de7c
# ╠═4ad1783a-afbc-4ab7-b043-4931569c330b
# ╠═150106d9-ce77-4729-ad9c-20d48c6361b9
# ╠═bdc1c037-4295-4e09-95e7-a1837b7317da
# ╠═b70dda56-903f-47c0-bd42-c6ecca40c8d2
# ╠═016af154-833c-4024-96e1-97ba3de9e487
# ╠═789c8b60-06d8-45f7-9131-8710d51f5c7f
# ╠═1250bdfd-9111-494b-b823-208cdf0eb172
# ╠═5f0169be-4363-42c3-9e12-315edcf37a7b
# ╠═b8f86151-a304-4d88-843d-85ba86908a72
# ╠═263b7313-2246-42bd-84ae-a032fe6e2a85
# ╠═be2ccdb1-b310-46d4-b319-ce51b45399a0
# ╠═0795f1c3-832c-4b38-a029-960b9f7099e5
# ╠═4ea9933e-cef6-4638-98c9-7a77309916ac
# ╠═404bd6e3-a8b0-4542-8ebe-810ef26dcb47
# ╠═ccc7eb5e-6275-4135-ba44-8fd6862f49e4
# ╠═d4ddfacd-d430-4a7d-b503-9b9f9a260cba
# ╠═dcc4bcde-ca0e-490b-a56f-96997f72e667
# ╠═60db43d2-9e98-4e11-a3e9-12cd653190c9
# ╠═72f6f45f-2e89-4a4b-872a-d47f932effd4
# ╠═a6ef9143-0fcf-4af2-a80e-83b79ac50cd2
# ╠═b1ddcaf2-fa92-4eda-a509-0c3792316572
# ╠═00f748d3-ad88-42c5-9b2c-6b443bd655c1
# ╠═fa896da9-ff1f-4e62-8550-84bf958cd5f9
# ╠═e69786b0-8641-42b4-bd64-320787dff8a1
# ╠═dbd57c19-c91c-4407-b798-c9b1879a3c10
# ╠═7ba1a20e-4a3c-4e3e-be3e-9b90ebc54d37
# ╠═52d2c069-6a8d-4bf4-be2b-223f6eb9e72e
# ╠═8a079fc0-6bee-48ef-93cc-27804134b935
# ╠═1b9f27bc-934f-48c4-86b1-37f541cb8a0c
# ╠═de02cc62-f6f4-450b-a002-c5349c7f8e01
# ╠═17cbf51c-645f-49fc-be63-ae60ce86923a
# ╠═e017637e-4885-4f32-92be-6b148c3e603b
# ╠═bb18bce4-140b-459d-b7c3-1bf29621ae37
# ╠═88f5c027-4568-4c80-83ee-b707d6f72e0f
# ╠═61c98942-8299-4dea-bafb-b728fd567f28
# ╠═9b9746c3-096a-437a-be05-fa4daba0509a
# ╠═938618cc-cf5b-45ff-8ad6-69a6dc0a99f0
# ╠═6776c519-95f8-4d7f-975c-6a78b8e9e6c6
# ╠═8b03a23f-51a6-4a74-b6aa-1cf74c5a6677
# ╠═ceaae763-3a4f-401a-a980-ce0678f11d37
# ╠═e6b48621-eb55-4f93-80fe-8bde497512b3
# ╠═f434f5c1-16f5-47e0-b772-74da18fefbdb
# ╠═d094bb28-ce6a-496f-a1b9-17e0db817085
# ╠═3a262d65-bd0d-470c-a738-f959e5e89143
# ╠═3c253f05-eb91-4e6f-a633-04912c36bb52
# ╠═bedc9813-4204-4a9b-b784-55ea6fed5364
# ╠═6f83a0a8-f919-4830-90f9-5320914cba3c
# ╠═81ccfc68-6dbc-402c-af18-871acf4022a5
# ╠═eba5d743-da8d-4e2c-a908-961dacf24264
# ╠═bf6b6f23-0878-4a52-bb13-b2faa5a2fa3b
# ╠═d85603c0-1e4d-4a92-ba00-e1a602c533da
# ╠═32356218-ecc6-4f2c-90a2-67244dcca04e
# ╠═7af934e0-7c7c-497a-852d-343200b6e5a0
# ╠═4194b150-1874-4ae2-b153-e3761cb724e1
# ╠═80342fe0-397b-4681-9178-8ab2bde90323
# ╠═af2ed6ef-17ee-4613-9501-05f331bea410
# ╠═31909481-bc95-49b3-8081-071cd888814f
# ╠═f71f9b5e-9654-49ac-803b-169652cc78e9
# ╠═94145064-52b0-4688-8947-6a1b3b34c553
# ╠═b027d115-648f-4c53-9448-0d1fbdd854fa
# ╠═304730fb-f35f-4cb6-969d-1349215ce8d7
# ╠═cf3fa8c1-3734-4a78-bb4b-a1d4f14c04e9
# ╠═04b5272f-3c0b-417a-bd3a-15b2716b94ec
# ╠═886a3755-105b-42d3-97c1-2e3f4fdf8541
# ╠═aa84cb6b-f7d1-428b-be48-bb0752fdc71e
# ╠═7d72b8eb-5f99-43c9-b540-5231f9a7d40f
# ╠═da48268b-5875-4d70-bb11-c57f0680218e
# ╠═06a27299-6c51-42a0-9c31-98587898d895
# ╠═f72e8f8f-b532-497e-892e-6708c83edab5
# ╠═c3392f6c-5883-40e7-80e4-f6d752420837
# ╠═ca843464-5ba2-4d1c-ab64-5214b01fe305
# ╠═945c160f-ff97-428d-8a98-b8d59cb6afe9
# ╠═08b9835f-5600-409a-8d67-3fc599efbc03
# ╠═de4ba031-9cc2-47b6-9cc9-45850ca3befd
# ╠═e05889f7-d92a-4e11-b0fa-7dce89ad7f96
# ╠═f138441b-410b-4ef3-85cf-3e052ce90351
# ╠═0a7000be-88dd-4cba-a0b8-6f82ff1c9c07
# ╠═308ea3b1-8145-49a7-aa16-048666312620
# ╠═d76052a6-92a3-416e-89f7-f37e160b7d54
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
