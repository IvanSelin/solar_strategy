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

# ╔═╡ c6a2816d-e81e-449e-af8c-93676c3fd077
begin
	using DataFrames
	using CSV
	using Plots # default
	import PlotlyJS as plotjs
	# using PlotlyBase
	# using PlotlySave
	using TimeZones
	using Dates
	using Optim
	using LineSearches
	using PlutoUI
	using Peaks
	using WebIO
	using StatsBase
	using Distributions
	# using LinearAlgebra
	# using Ranges
	using Random
	using ProgressMeter
end

# ╔═╡ 3e46a49b-7bb7-4889-8c94-b842977899e4
begin
	include("../src/energy_draw.jl")
	include("../src/time.jl")
	include("../src/solar_radiation.jl")
	include("../src/track.jl")
	include("../src/utils.jl")
	include("../src/strategy_calculation.jl")
	include("../src/weather.jl")
end

# ╔═╡ 94da1f00-efe7-11ed-2d94-f9a905085f40
Threads.nthreads()

# ╔═╡ 7282988b-0667-409d-9634-d874e7767d16
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

# ╔═╡ b20a1cbd-705e-4b49-b671-d042d1511afe
PlutoUI.TableOfContents()

# ╔═╡ d5213897-9cb6-453a-bc4e-a7015909c886
rng = MersenneTwister(1234)

# ╔═╡ 769d7af1-e95c-426f-aa97-85b219c5b65e
plotlyjs()

# ╔═╡ 64ecb4d9-b67c-4a41-a79e-300affd0440f
md"# План

План довольно прост. Надо доказать, что нам есть зачем выбирать разные скорости на разных участках трассы.

Шаги для осуществления плана:
1. Взять ровную трассу, подобрать одну скорость на все участки (отлично работает)
2. Взять ровную трассу, подобрать разные скорости на все участки (работает, но избыточно)
3. Взять трассу с холмом, подобрать одну скорость на все участки (плохо работает, нужен другой подход)
4. Взять трассу с холмом, подобрать разные скорости для каждого участка (должно норм работать)
5. Взять трассу ощутимой длины (участков 300-500) с холмами, чтобы оно считалось достаточно долго. Подобрать разные скорости обычной оптимизацией и сказать что долго выходит (по идее норм посчитает, но долго)
6. Взять эту же трассу с холмами, подобрать скорости моим методом (в идеале должно посчитать примерно так же, но быстрее)
"

# ╔═╡ 9529a535-3655-496c-8077-30c148341fb9
md"## Необходимые приготовления (функции и пр.)"

# ╔═╡ e3fa515b-10bd-4ae3-9c47-d8b5e365c4ea
md"### Переменная на все участки"

# ╔═╡ 1029bac5-c3ae-4dd9-b992-db077cf69c18
function single_optim(track, segments, start_energy, start_datetime)
	segments_length = size(segments, 1)
	
	function f_wrap_single(input_speed)
		# speed_vector = fill(first(input_speeds) / 3.6, segments_length)
		speed_vector = fill(input_speed / 3.6, segments_length)
		power_use_f, solar_power_f, time_s_f = solar_trip_boundaries(
			speed_vector, segments, start_datetime
		)
		energy_in_system_f = start_energy .+ solar_power_f .- power_use_f
		pushfirst!(energy_in_system_f, start_energy)
	
		# cost = sum(segments_clean.diff_distance ./ input_speeds) + 100 * (0. - last(energy_in_system_clean))^2;
		
		cost = sum(segments.diff_distance ./ speed_vector) + 1000 * abs(minimum(energy_in_system_f))^2 + 10 * (0. - last(energy_in_system_f))^2;
	
		# cost = sum(segments.diff_distance ./ input_speeds) + 100 * abs(minimum(energy_in_system_f))
		return cost
	end

	# td = TwiceDifferentiable(f_wrap_single, [speed]; autodiff = :forward)
	# lower_bound = fill(0.0, var_num)
	# upper_bound = fill(100.0, var_num)
	# tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)


	# result = optimize(td, tdc, speeds 
	# # .+ rand(vars_amount) .- 0.5
	#     ,
	#     IPNewton(),
	#     Optim.Options(
	#         x_tol = 1e-6,
	#         f_tol = 1e-6,
	#         g_tol = 1e-6
	#     )
	# )

	result = optimize(f_wrap_single, 0., 150.)
	minimized_speeds = fill(Optim.minimizer(result), segments_length)
	minimized_speeds_ms = minimized_speeds / 3.6

	power_use, solar_power, time_s = solar_trip_boundaries(
		minimized_speeds_ms, segments, start_datetime
	)
	# track points, not segments, that's why it is size is +1 
	energy_in_system_new = start_energy .+ solar_power .- power_use
	lowest_energy = minimum(energy_in_system_new)
	last_energy = last(energy_in_system_new)
	pushfirst!(energy_in_system_new, start_energy)

	track_plot = plot(track.distance, track.altitude, title="Track, $(segments_length) segments",
		color=:green,
		ylabel="altitude(m)")

	speed_plot = plot(
		get_mean_data(track.distance),
		minimized_speeds,
		seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed plot for 1 var, total time $(round(last(time_s), digits=3)) sec",
		ylabel="speed(kmh)"
	)

	low_energy_red = fill(0., size(track.distance, 1))

	energy_plot = plot(
		track.distance,
		[energy_in_system_new low_energy_red],
		linewidth=[1 3],
		title="Energy for 1 var, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)"
	)
	
	plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(650,700), legend=false)
	
end

# ╔═╡ 4b3430b1-619f-4123-9a5b-6db7e7083ee9
md"### Переменная на каждый участок"

# ╔═╡ 018eeb89-8463-414d-badb-f64e1105ffdd
# подбираем скорости на всех участках
function regular_optim(track, segments, speeds, start_energy, start_datetime, upper_speed_bound=150.)

	var_num = size(segments,1)
	function f_wrap_reg(input_speeds)
		speed_vector = input_speeds / 3.6;
		power_use_f, solar_power_f, time_s_f = solar_trip_boundaries(
			speed_vector, segments, start_datetime
		)
		energy_in_system_f = start_energy .+ solar_power_f .- power_use_f
		pushfirst!(energy_in_system_f, start_energy)
	
		# cost = sum(segments_clean.diff_distance ./ input_speeds) + 100 * (0. - last(energy_in_system_clean))^2;
		
		cost = sum(segments.diff_distance ./ speed_vector) + 1000 * abs(minimum(energy_in_system_f))^2 + 10 * (0. - last(energy_in_system_f))^2;
	
		# cost = sum(segments.diff_distance ./ input_speeds) + 100 * abs(minimum(energy_in_system_f))
		return cost
	end
	
	td = TwiceDifferentiable(f_wrap_reg, speeds; autodiff = :forward)
	lower_bound = fill(0.0, var_num)
	upper_bound = fill(upper_speed_bound, var_num)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)

	result = optimize(td, tdc, speeds 
	# .+ rand(vars_amount) .- 0.5
	    ,
	    IPNewton(),
	    Optim.Options(
	        x_tol = 1e-8,
	        f_tol = 1e-8,
	        g_tol = 1e-8
	    )
	)

	minimized_speeds = Optim.minimizer(result)
	minimized_speeds_ms = minimized_speeds / 3.6

	power_use, solar_power, time_s = solar_trip_boundaries(
		minimized_speeds_ms, segments, start_datetime
	)
	# track points, not segments, that's why it is size is +1 
	energy_in_system_new = start_energy .+ solar_power .- power_use
	lowest_energy = minimum(energy_in_system_new)
	last_energy = last(energy_in_system_new)
	pushfirst!(energy_in_system_new, start_energy)

	track_plot = plot(track.distance, track.altitude, title="Track, $(var_num) segments",
		color=:green,
		ylabel="altitude(m)")

	speed_plot = plot(
		get_mean_data(track.distance),
		minimized_speeds,
		seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed plot for $(var_num) vars, total time $(round(last(time_s), digits=3)) sec",
		ylabel="speed(kmh)"
	)

	low_energy_red = fill(0., size(track.distance, 1))
	
	energy_plot = plot(
		track.distance,
		[energy_in_system_new low_energy_red],
		linewidth=[1 3],
		title="Energy for $(var_num) vars, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)"
	)
	
	plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(650,700), legend=false)

end

# ╔═╡ 2274f3b1-6a63-4206-a9ab-fb38c7a638f3
function simulate_run(speeds, track, segments, start_energy, start_datetime)
	minimized_speeds_ms = speeds / 3.6
	
	power_use, solar_power, time_s = solar_trip_boundaries(
		minimized_speeds_ms, segments, start_datetime
	)
	# track points, not segments, that's why it is size is +1 
	energy_in_system_new = start_energy .+ solar_power .- power_use
	lowest_energy = minimum(energy_in_system_new)
	last_energy = last(energy_in_system_new)
	pushfirst!(energy_in_system_new, start_energy)

	track_plot = plot(track.distance, track.altitude, title="Track",
		color=:green,
		ylabel="altitude(m)", label="track")

	speed_plot = plot(
		get_mean_data(track.distance),
		speeds,
		seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed plot with manually set speeds, total time $(round(last(time_s), digits=3)) sec",
		ylabel="speed(kmh)",
		label="speed"
	)

	low_energy_red = fill(0., size(track.distance, 1))
	
	energy_plot = plot(
		track.distance,
		[energy_in_system_new low_energy_red],
		linewidth=[1 3],
		title="Energy with manually set speeds, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)",
		label="energy"
	)
	
	plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(750,700), legend=false)
end

# ╔═╡ df5b2134-2bd7-4e07-91ff-e5ec736f2dfa
md"# Ровная трасса"

# ╔═╡ e9a6a82a-ecde-4baf-ad2c-f0a858bb2ed0
track_flat, segments_flat = get_track_and_segments("../data/data_test_flat.csv");

# ╔═╡ 9e47be95-5f9a-43d6-876e-abb33a67de0a
track_flat

# ╔═╡ 9de8094a-5d24-4b08-8879-c7853cf6954b
segments_flat

# ╔═╡ b965ae63-8e6b-4c0a-91b6-c29bae610590
single_optim(
	track_flat,
	segments_flat,
	75.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ c892b78b-58dd-4a61-b731-6f4c1432ea1a
md"Получили одну скорость в 79.35 км/ч, на которой проезжаем всю дистанцию за 226.831 секунд"

# ╔═╡ c9daaf07-3dc2-4bef-ae23-b69acd2cb513
regular_optim(
	track_flat,
	segments_flat,
	fill(75., size(segments_flat,1)),
	75.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ a4ebc847-08fd-4976-87a7-2c35ca7920be
md"Получили 5 примерно одинаковых скоростей около тех же самых 79 км/ч, на которых проезжаем всю дистанцию опять за 226 секунд

Особой разницы нет, и в такой ситуации рассматривать каждый участок трассы как отдельный нет смысла"

# ╔═╡ 8f1d25de-f4dd-481a-a7f5-b8691586e705
md"# Трасса с холмом"

# ╔═╡ 46a01c40-c920-483c-aa93-1297eb0cddc8
md"А теперь попробуем более сложную конфигурацию. Допустим, нам надо преодолеть холм, т.е. на трассе будет подъём и спуск"

# ╔═╡ 7f3807d7-e3ef-4c8f-8a21-8abdd4f25259
track_hill, segments_hill = get_track_and_segments("../data/data_test_hill.csv");

# ╔═╡ 9f151595-8f5a-4337-90fb-a8ef10f453b3
track_hill

# ╔═╡ 67e01c26-832a-4eda-863d-e1842aae13dd
segments_hill

# ╔═╡ 074e1eeb-98db-4728-8bfe-4b4ed91b554f
single_optim(
	track_hill,
	segments_hill,
	75.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 85d293fc-2992-4a38-a080-93179ff4a9f1
md"При попытке поиска оптимального режима движения для холмистого варианта, скорость стала гораздо ниже, всего лишь 25 км/ч и итоговое время составило 704 секунды. Но почему?

Потому что на подъём тратится большое количество энергии. И если поехать хоть немного быстрее, то количество энергии в системе упадёт ниже 0, чего быть не может.

Решение получается далёким от оптимального. Это очевидно, т.к. в конце трассы остаётся много неиспользованной энергии. А это значит, что можно ехать быстрее

Попробуем изменять скорость на каждом участке по отдельности"

# ╔═╡ b63f2571-079b-4ef3-9f96-6a465c255733
regular_optim(
	track_hill,
	segments_hill,
	fill(45., size(segments_hill,1)),
	75.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 48291326-f37d-4397-aeb0-fee0855970d8
md"Стало гораздо лучше!

Теперь трасса проезжается за 400+ секунд. Скорости при этом разнятся от ~15 км/ч в самом начале, до ~150 км/ч в конце. При этом фундаментальные ограничения не нарушаются."

# ╔═╡ 3509ebf4-b758-4ce3-a69f-fb20c16615ef
md"Однако, в предложенном плане движения очень высокие скорости в конце, которые далеко не всегда могут быть реализованы. Имеет смысл ограничить сверху максимальную скорость."

# ╔═╡ 72773944-4541-4220-8f0a-151ccd67bcaa
regular_optim(
	track_hill,
	segments_hill,
	fill(45., size(segments_hill,1)),
	75.,
	DateTime(2022,1,1,10,0,0),
	80.
)

# ╔═╡ d6921b5f-d9ae-4886-a40e-98f9ad052ed2
md"Как видно, ограничив максимальную скорость 80 км/ч, итоговое время не сильно пострадало, но план получился гораздо более реалистичным."

# ╔═╡ 49c58d5c-11e2-4aff-a408-cd4aabff7312
md"### Лирическое отступление"

# ╔═╡ 53c9b844-5d1e-4444-abe2-ec16f7cb1ee9
md"Доделать:

1. более унифицированные трассы с и без холма - сделано
2. энергия до холма и после осталась примерно такой же. что это означает? Это значит, что на значениях энергии не около нуля, на это всё равно. теперь надо проверить, можно ли проехать с определённой скоростью, чтобы было такое же время и энергия прохождения дистанции - нужно только для ворлд солар челленджа. пока без этого"

# ╔═╡ e6540edc-e971-48e9-8679-689a3ed73e7b
simulate_run(
	fill(47.81, size(segments_hill, 1)),
	track_hill,
	segments_hill,
	75.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 0b1564a2-8129-43ae-8831-77e18959be67
md"Получается быстрее за то же самое количество энергии

Выходит, просто не стоит париться за нижний порог энергии?

Нужно, если брать как одну большую задачу"

# ╔═╡ 53e5f786-6e6a-4892-8252-aa3389682a80
md"# Длинная трасса с повторяющимися холмами"

# ╔═╡ afbea56b-194a-407e-9f74-8a3c9d0cbbed
track_hills, segments_hills = get_track_and_segments("../data/data_test_hills.csv");

# ╔═╡ db0b009d-1abd-42c0-bbad-b6ba53af5fb5
track_hills

# ╔═╡ f664a938-2e64-4e37-9404-578720b4bd57
segments_hills

# ╔═╡ ec45bd26-942d-4180-81f5-a9db7e227c4f
single_optim(
	track_hills,
	segments_hills,
	250.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ d82662d3-5533-49cf-967f-5bfa1fb65b13
regular_optim(
	track_hills,
	segments_hills,
	fill(75., size(segments_hills,1)),
	250.,
	DateTime(2022,1,1,10,0,0),
	150.
)

# ╔═╡ d907fa8d-f398-48ee-bd2d-e1fea7149fb1
md"# Длинная трасса с разными холмами"

# ╔═╡ de68d838-4e10-49fa-a92b-a5524ef378d7
track_hills2, segments_hills2 = get_track_and_segments("../data/data_test_hills2.csv");

# ╔═╡ fc78cf3e-c12c-43cb-9ce1-07665134b2cc
single_optim(
	track_hills2,
	segments_hills2,
	500.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ eb3bd221-532a-4299-88bc-3a19ccb4c7e3
regular_optim(
	track_hills2,
	segments_hills2,
	fill(50., size(segments_hills2,1)),
	500.,
	DateTime(2022,1,1,10,0,0),
	150.
)

# ╔═╡ 0c1b0b93-23e0-422a-a6ce-30b1c6b1ab51
md"# Длинная трасса с разными холмами (финиш выше старта)"

# ╔═╡ f7a72d85-ec30-4406-9a9e-a148d2244cf2
track_hills3, segments_hills3 = get_track_and_segments("../data/data_test_hills3.csv");

# ╔═╡ bb34bc3a-4232-4c52-8358-08c571566055
single_optim(
	track_hills3,
	segments_hills3,
	500.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ cd45b86d-7231-418f-a43c-b78b4fa2d879
regular_optim(
	track_hills3,
	segments_hills3,
	fill(50., size(segments_hills3,1)),
	500.,
	DateTime(2022,1,1,10,0,0),
	200.
)

# ╔═╡ 69bf6ffc-8c79-4573-bfb2-d264c66a7040
track_hills3a, segments_hills3a = get_track_and_segments("../data/data_test_hills3a.csv");

# ╔═╡ e410277c-87ef-4f18-89cb-7ce8817920b5
single_optim(
	track_hills3a,
	segments_hills3a,
	500.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 6e74e3be-fb37-4bdd-a059-87b513f8382a
regular_optim(
	track_hills3a,
	segments_hills3a,
	fill(50., size(segments_hills3a,1)),
	500.,
	DateTime(2022,1,1,10,0,0),
	200.
)

# ╔═╡ 08027dc4-0f87-4c31-8c30-703314a0c512
md"Пока очень похоже на то, что скорости подбираются следующим образом:

Подобрать максимально возможную скорость одной переменной на всю дистанцию, чтобы минимум энергии не был меньше 0.

Зафиксировать скорости до точки, где энергия стала равной нулю. Далее рассматривать только участки после этой точки.

Подбирать такую скорость одной переменной на всю оставшуюся дистанцию, чтобы минимум энергии не был меньше 0. Зафиксировать. Далее рассматривать только после этой точки.

И так до тех пор, пока точка не будет финишем
"

# ╔═╡ b4e85535-d731-4e90-bf33-24ba6ee69aa4
md"# Австралия с сокращённым количеством точек"

# ╔═╡ 940a6c60-a5e5-4bab-940b-4a3b027ecc65
begin
	track_aus, segments_aus = get_track_and_segments("../data/data_australia_random.csv");
	track_aus.altitude = track_aus.altitude * 10;
	segments_aus = get_segments_for_track(track_aus);
end

# ╔═╡ 7adbc40f-8f0f-4434-96da-3962f7107ce9
single_optim(
	track_aus,
	segments_aus,
	5100.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 52c423ab-248a-4398-af19-6f80acd6164f
regular_optim(
	track_aus,
	segments_aus,
	fill(40., size(segments_aus,1)),
	5100.,
	DateTime(2022,1,1,10,0,0),
	150.
)

# ╔═╡ f8d9184d-e88e-430c-80e9-e2a660fce6ba
md"Похоже что на длинных дистанциях работает нормально. Чем больше перепад высот, тем лучше работает"

# ╔═╡ aae20b83-351f-44ef-8cc3-79404167fdc4
md"Идея 2: выигрыш будет тем больше, чем раньше будет самый большой холм. Что это значит: можно подбирать скорость до первого самого большого холма (где энергия 0), потом подбирать её опять, но уже начиная с момента нулевой энергии. До следующего нулевого холма)"

# ╔═╡ 9c686f0c-831a-4d5c-b18c-d5a459fe871a
md"Похоже что пора делать осадки, ну или увеличивать количество получаемой энергии, чтобы была разница в режимах движения"

# ╔═╡ 0eec4f71-31a7-435d-9ed5-1c86c2751f5b
md"# Начинаем делать осадки"

# ╔═╡ 3b0cd1f5-b2f2-4178-8f0e-156f2e6ffd10
mapbox_style="open-street-map"
# mapbox_style="stamen-terrain"
# "open-street-map", "carto-positron", "carto-darkmatter", "stamen-terrain", "stamen-toner" or "stamen-watercolor" yield maps composed of raster tiles from various public tile servers which do not require signups or access tokens

# ╔═╡ 9ad47da6-2b71-4354-9bd5-7113f98139d7
md"Начнём с отображения дистанции, для того чтобы моделировать осадки"

# ╔═╡ 6bc90ca1-dc8d-41e2-a9df-b54ccac74f79
plotjs.Plot(
	plotjs.scattergeo(
		lat=track_aus.latitude,
		lon=track_aus.longitude
	),
	plotjs.Layout(
		width=700,
		height=600,
		# geo_scope="world",
		# geo_resolution=50,
		geo_fitbounds="locations", # geo_fitbounds="geojson",
		mapbox_style="stamen-terrain"
	)
)

# ╔═╡ 57f47d65-5a3e-4c86-be15-587d94fbf677
md"## Создание данных осадков

Осадки у нас локализованы по:
1. Времени
2. Месту

Это значит что мне надо генерировать некоторую пространственную сетку с разными значениями в узлах (место), и иметь кучу таких сеток на каждый момент времени (время)

План таков:
1. Сделать статичную карту осадков, отобразить её
2. Модифицировать код, чтобы он её использовал
3. Помножить на локализацию по времени
4. Опять модифицировать код, чтобы эти данные использовались"

# ╔═╡ ac174dbd-11d5-4859-b8f5-751af94f9ac6
md"## Генерируем карту сетку и накладываем осадки"

# ╔═╡ e5088f3d-87f1-41e8-b8e6-e5a98b1c6e83
md"Пробуем сперва даже без сетки!"

# ╔═╡ 64628a15-fb96-4e65-ac4d-d07d9dab7fcb
md"## Функция генерация осадков на сетке"

# ╔═╡ 5709b86c-6e68-4dac-b7d8-a880e9eca00d
function generate_clouds(
	lat_from,
	lon_from,
	lat_to,
	lon_to,
	lat_peak,
	lon_peak,
	lat_std,
	lon_std,
	ndims,
	coef
)
	lat_distr = rand(Normal(lat_peak, 2), 10000 * ndims)
	lon_distr = rand(Normal(lon_peak, 2), 10000 * ndims)

	edges_lat = range(min(lat_from, lat_to), max(lat_from, lat_to), ndims + 1) #lat_from:step:lat_to
	edges_lon = range(min(lon_from, lon_to), max(lon_from, lon_to), ndims + 1)#lon_from:step:lon_to

	hist = fit(
		Histogram,
		(lat_distr, lon_distr),
		(edges_lat, edges_lon)
	)

	# return hist
	normed_weights = hist.weights / maximum(hist.weights) * coef#, edges_lat, edges_lon

	edges_lat_collected = collect(edges_lat)
	edges_lon_collected = collect(edges_lon)
	
	if lat_from > lat_to
		reverse!(normed_weights, dims=1)
		reverse!(edges_lat_collected)
	end

	if lon_from > lon_to
		reverse!(normed_weights, dims=2)
		reverse!(edges_lon_collected)
	end

	return normed_weights, edges_lat_collected, edges_lon_collected
	
	#hist.weights, hist.edges[1], hist.egdes[2]
end

# ╔═╡ 8ade2b34-98e9-495c-b9e0-faa655652cef
ndims = 15

# ╔═╡ 6cab74b6-650e-4d6f-b6a9-0cbc0282f104
w, elat, elon = generate_clouds(
	-10,
	130,
	-18,
	135,
	-12.5,
	131.2,
	0.5,
	0.5,
	ndims,
	0.75
);

# ╔═╡ 5fca02b2-24ae-4081-98b0-e7f640e65405
function generate_density_mapbox(w, edges_lat, edges_lon)
	ndims = size(w,1)
	w_arr = collect(Iterators.flatten(w))
	edges_lat_rep = repeat(get_mean_data(edges_lat), outer=ndims)
	edges_lon_rep = repeat(get_mean_data(edges_lon), inner=ndims)

	df = DataFrame(lat=edges_lat_rep, lon=edges_lon_rep, z=w_arr)
	traces_vector::AbstractVector{plotjs.AbstractTrace} = [];

	push!(
		traces_vector,
		plotjs.densitymapbox(
			lat=df.lat,
			lon=df.lon,
			z=df.z,
			opacity=0.5
		)
	)
	
	push!(
		traces_vector,
		plotjs.scattermapbox(
			lat=track_aus.latitude,
			lon=track_aus.longitude,
			marker_color="red",
			marker_size=1,
			mode="lines"
		)
	)
	
	plotjs.Plot(
		traces_vector,
		plotjs.Layout(
			width=650,
			height=600,
			geo_fitbounds="locations",
			autosize=true,
			# mapbox_style="stamen-terrain"
			mapbox_style=mapbox_style,
			mapbox_center_lat=-25.0,
			mapbox_center_lon=132.0,
			mapbox_zoom=3
		)
	)
end

# ╔═╡ e4188c54-b620-48cd-b1ff-448ac4a15950
generate_density_mapbox(w, elat, elon)

# ╔═╡ 7bf6b5ed-ed8b-42c9-8772-38f4e2b3a0f9
md"Выглядит правдиво, осталось только нормально настроить для норм внешнего вида"

# ╔═╡ 7bd787da-997c-48bc-8c5b-f181148ac964
md"## Heat map (self-made)"

# ╔═╡ 8e8ca19d-a381-4055-aef7-f27292aac611
function generate_heatmap_traces(w, edges_lat, edges_lon)
	traces_vector::AbstractVector{plotjs.AbstractTrace} = [];
	push!(
		traces_vector,
		plotjs.scattermapbox(
			lat=track_aus.latitude,
			lon=track_aus.longitude,
			marker_color="red",
			marker_size=1,
			mode="lines"
		)
	)
	for i=1:length(edges_lat)-1
		for j=1:length(edges_lon)-1
			# println("lat: $(edges_lat[i]), lon: $(edges_lon[j]), w: $(w[i,j])")
			# println("lat+1: $(edges_lat[i+1]), lon+1: $(edges_lon[j+1]), w+1: $(w[i,j])")
			trace = plotjs.scattermapbox(
				fill="toself",
				lat = [
					edges_lat[i],
					edges_lat[i],
					edges_lat[i+1],
					edges_lat[i+1],
					edges_lat[i]
				],
				lon = [
					edges_lon[j],
					edges_lon[j+1],
					edges_lon[j+1],
					edges_lon[j],
					edges_lon[j]
				],
				# lat=[-11.,-11.,-12.,-12.5,-12., -12.],
				# lon=[130.,131.,131.,130.5,130., 130.],
				marker_size=1,
				# marker_color="orange",
				opacity=0.5,
				showlegend=false,
				# marker_colorscale=w[i,j]
				# marker_colorscale="Viridis",
				# marker_color=w[i,j]
				marker_color="rgb($(w[i,j]*255),$((1-w[i,j])*255),0)",
				name="$(w[i,j])"
			);
			push!(traces_vector, trace)
		end
		# println()
	end
	plotjs.Plot(
		traces_vector,
		plotjs.Layout(
			width=700,
			height=600,
			geo_fitbounds="locations",
			mapbox_style=mapbox_style,
			autosize=true,
			mapbox_center_lat=-25.0,
			mapbox_center_lon=132.0,
			mapbox_zoom=3
		)
	)
end

# ╔═╡ 75f0759b-285c-402b-afdf-243cf5d81bbe
generate_heatmap_traces(w, elat, elon)

# ╔═╡ e7575b4b-49bf-4498-b74b-69bf9d1762cb
begin
	w_test, edges_lat_test, edges_lon_test = generate_clouds(
		-10,
		125,
		-35,
		145,
		-20.,
		134.,
		0.5,
		0.5,
		15,
		0.75
	);
	generate_heatmap_traces(
		w_test, edges_lat_test, edges_lon_test
	)
end

# ╔═╡ 45183999-58ff-4bdd-a4b0-225b1e50a5ee
md"### Регулировка пика"

# ╔═╡ 00973819-68b4-4d7e-8c2a-bebfce2a820a
@bind lat_peak Slider(-35:0.1:-10, default=-20)

# ╔═╡ 3498e807-a860-4f73-9d97-68c6f46d7efe
lat_peak

# ╔═╡ b1449601-ce61-4097-8b32-345275f859b8
@bind lon_peak Slider(125:0.1:145, default=134)

# ╔═╡ aba59435-42d9-4ac8-b768-0f017591d788
lon_peak

# ╔═╡ 63afb4e4-dc19-4106-8bb5-b51d0881675d
begin
	w_test2, edges_lat_test2, edges_lon_test2 = generate_clouds(
		-10,
		125,
		-35,
		145,
		lat_peak,
		lon_peak,
		0.5,
		0.5,
		50,
		0.75
	);
	generate_density_mapbox(
		w_test2, edges_lat_test2, edges_lon_test2
	)
end

# ╔═╡ 642a882f-de23-4414-a579-b61088e98bd9
w_test2

# ╔═╡ 78327f2d-7cf2-49e3-98fa-d90aeb479c45
md"Density гораздо быстрее, его лучше использовать для анимаций"

# ╔═╡ 760e7514-e3c4-4366-93fc-e1379f96c2ad
md"Можно приступать к использованию сгенерированных облаков"

# ╔═╡ e24f1d14-e0d9-45e8-a415-172e56e7cc9c
md"А потом и к распространению во времени"

# ╔═╡ a803b045-f7ce-4a14-b363-1c166e34fbe7
md"## Обновлённый функционал со статичными осадками"

# ╔═╡ e07dc08f-9cd8-46b6-8e92-062372b955a3
weather_segm = calculate_weather_weights_for_segments(
	w_test2,
	edges_lat_test2,
	edges_lon_test2,
	segments_aus
)

# ╔═╡ 39e09aba-c4c5-46b1-8ef9-63403b0225eb


# ╔═╡ f7d29865-473d-47f5-a40c-08f48b531401
function solar_trip_weather(input_speed, segments, start_datetime,
	weather_weights, weather_edges_lat, weather_edges_lon
)
    # input speed in m/s
	# @debug "func solar_trip_calculation_bounds input_speed size is $(size(input_speed, 1)), track size is $(size(track.distance, 1)))"

    # calculating time needed to spend to travel across distance
    # time_df = calculate_travel_time_datetime(input_speed, segments, start_datetime)

    #### calculcations
    # mechanical calculations are now in separate file
    mechanical_power = mechanical_power_calculation_alloc.(input_speed, segments.slope, segments.diff_distance)

    # electical losses
    electrical_power = electrical_power_calculation(segments.diff_distance, input_speed)
    # converting mechanical work to elecctrical power and then power use
    # power_use = calculate_power_use(mechanical_power, electrical_power)
    power_use_accumulated_wt_h = mechanical_power + electrical_power
	cumsum!(power_use_accumulated_wt_h, power_use_accumulated_wt_h)
	power_use_accumulated_wt_h = power_use_accumulated_wt_h / 3600.

    # get solar energy income
	# @debug "track size is $(size(track.latitude, 1))"

	time_seconds = calculate_travel_time_seconds(input_speed, segments)
	mean_seconds = get_mean_data(time_seconds)
	pushfirst!(mean_seconds, 0.)
	milliseconds = round.(mean_seconds .* 1000)
	mean_segment_utc = start_datetime .+ Dates.Millisecond.(milliseconds)
	
    solar_power = solar_power_income_alloc.(
		segments.latitude,
		segments.longitude, 
		segments.altitude, 
		# time_df.utc_time,
		mean_segment_utc,
		segments.diff_distance,
		input_speed
	)
	weather_coef = zeros(size(segments.latitude,1))
	for i in eachindex(weather_coef)
		_, lat_index = findmin(abs.(weather_edges_lat .- segments.latitude[i]))
		_, lon_index = findmin(abs.(weather_edges_lon .- segments.longitude[i]))
		weather_coef[i] = 1 .- weather_weights[lat_index, lon_index]
		# println("lat $(segments.latitude[i]), lon $(segments.longitude[i]), w is $(weather_coef[i])")
	end
	
	solar_power_adjusted = solar_power .* weather_coef
	# println()
	# print(weather_coef)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power_adjusted)

    # TODO: calculate night charging - do it later since it is not critical as of right now
    time_seconds = calculate_travel_time_seconds(input_speed, segments)
    return power_use_accumulated_wt_h, solar_power_accumulated, time_seconds
end

# ╔═╡ 525a59f1-e377-4bb4-a922-51fb940ad688
pua_weather,spa_weather,ts_weather = solar_trip_weather(
	fill(55., size(segments_aus.longitude, 1)),
	segments_aus,
	DateTime(2022,1,1,10,0,0),
	w_test2,
	edges_lat_test2,
	edges_lon_test2
)

# ╔═╡ 1e9d17e9-d9fb-4b25-8233-4fc2fa6ce23e
plot(spa_weather, title="Накопленная солнечная энергия с учётом осадков")

# ╔═╡ 63bf3f77-6abe-4280-ab2f-520beb93a353
pua_no_weather,spa_no_weather,ts_no_weather = solar_trip_boundaries(
	fill(55., size(segments_aus.longitude, 1)),
	segments_aus,
	DateTime(2022,1,1,10,0,0)
);

# ╔═╡ ac2a21e4-491b-4b92-826c-5d4d74cf9ef1
plot(spa_no_weather, title="Накопленная солнечная энергия без учёта осадков")

# ╔═╡ 8d5be066-4625-4b60-a789-557ff2bd4660
plot(diff(spa_weather), title="Cолнечная энергия (абсолютная) с учётом осадков")

# ╔═╡ 8ae2eaaa-4a63-4c01-b0f7-0050e5e41b04
plot(diff(spa_no_weather), title="Солнечная энергия (абсолютная) без учёта осадков")

# ╔═╡ 7cba5ab0-8844-454d-8f7d-1fe5552f201c
plot(diff(spa_no_weather) .- diff(spa_weather),
	title="Разница в количестве солнечной энергии"
)

# ╔═╡ 6ae05363-c7b9-492d-91f0-ebc7935e40db
md"Получили меньше энергии на 3-10 участке

А это

от $(segments_aus.latitude[3]); $(segments_aus.longitude[3])

до $(segments_aus.latitude[10]); $(segments_aus.longitude[10]) "

# ╔═╡ b1a41a59-3842-468f-af9c-84387129ca45
weather_difference = last(spa_no_weather)-last(spa_weather)

# ╔═╡ b184d70e-9cdf-4531-87fc-6b606c54a8e8
md"Всего разница в энергии:

С осадками: $(last(spa_weather))

Без осадков: $(last(spa_no_weather))

Разница: $(weather_difference), т.е. $(weather_difference / last(spa_no_weather) * 100)%"

# ╔═╡ 52c09620-90cd-4e67-b84c-c7458d67c2c2
md"### Функция с погодой, 1 переменная на все"

# ╔═╡ 036ba7a9-8695-4022-a97b-dec4b3a08921
function single_optim_weather(track, segments, start_energy, start_datetime,
weather_weights, weather_edges_lat, weather_edges_lon)
	segments_length = size(segments, 1)
	
	function f_wrap_single(input_speed)
		# speed_vector = fill(first(input_speeds) / 3.6, segments_length)
		speed_vector = fill(input_speed / 3.6, segments_length)
		power_use_f, solar_power_f, time_s_f = solar_trip_weather(
			speed_vector, segments, start_datetime,
			weather_weights, weather_edges_lat, weather_edges_lon
		)
		energy_in_system_f = start_energy .+ solar_power_f .- power_use_f
		pushfirst!(energy_in_system_f, start_energy)
	
		# cost = sum(segments_clean.diff_distance ./ input_speeds) + 100 * (0. - last(energy_in_system_clean))^2;
		
		cost = sum(segments.diff_distance ./ speed_vector) + 1000 * abs(minimum(energy_in_system_f))^2 + 10 * (0. - last(energy_in_system_f))^2;
	
		# cost = sum(segments.diff_distance ./ input_speeds) + 100 * abs(minimum(energy_in_system_f))
		return cost
	end

	# td = TwiceDifferentiable(f_wrap_single, [speed]; autodiff = :forward)
	# lower_bound = fill(0.0, var_num)
	# upper_bound = fill(100.0, var_num)
	# tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)


	# result = optimize(td, tdc, speeds 
	# # .+ rand(vars_amount) .- 0.5
	#     ,
	#     IPNewton(),
	#     Optim.Options(
	#         x_tol = 1e-6,
	#         f_tol = 1e-6,
	#         g_tol = 1e-6
	#     )
	# )

	result = optimize(f_wrap_single, 0., 150.)
	minimized_speeds = fill(Optim.minimizer(result), segments_length)
	minimized_speeds_ms = minimized_speeds / 3.6

	power_use, solar_power, time_s = solar_trip_weather(
		minimized_speeds_ms, segments, start_datetime,
		weather_weights, weather_edges_lat, weather_edges_lon
	)
	# track points, not segments, that's why it is size is +1 
	energy_in_system_new = start_energy .+ solar_power .- power_use
	lowest_energy = minimum(energy_in_system_new)
	last_energy = last(energy_in_system_new)
	pushfirst!(energy_in_system_new, start_energy)

	track_plot = plot(track.distance, track.altitude, title="Track, $(segments_length) segments",
		color=:green,
		ylabel="altitude(m)")

	speed_plot = plot(
		get_mean_data(track.distance),
		minimized_speeds,
		seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed plot for 1 var, total time $(round(last(time_s), digits=3)) sec",
		ylabel="speed(kmh)"
	)

	low_energy_red = fill(0., size(track.distance, 1))

	energy_plot = plot(
		track.distance,
		[energy_in_system_new low_energy_red],
		linewidth=[1 3],
		title="Energy for 1 var, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)"
	)
	
	plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(650,700), legend=false)
	
end

# ╔═╡ c93555bf-a560-4d52-9545-c18af2749faa
md"### Функция с погодой, переменная на каждый участок"

# ╔═╡ bd36b76f-51bb-43c3-a94f-eb675fa3ac75
# подбираем скорости на всех участках
function regular_optim_weather(track, segments, speeds, start_energy, start_datetime, 
	weather_weights, weather_edges_lat, weather_edges_lon,
	upper_speed_bound=150.
)

	var_num = size(segments,1)
	function f_wrap_reg(input_speeds)
		speed_vector = input_speeds / 3.6;
		power_use_f, solar_power_f, time_s_f = solar_trip_weather(
			speed_vector, segments, start_datetime,
			weather_weights, weather_edges_lat, weather_edges_lon
		)
		energy_in_system_f = start_energy .+ solar_power_f .- power_use_f
		pushfirst!(energy_in_system_f, start_energy)
	
		# cost = sum(segments_clean.diff_distance ./ input_speeds) + 100 * (0. - last(energy_in_system_clean))^2;
		
		cost = sum(segments.diff_distance ./ speed_vector) + 1000 * abs(minimum(energy_in_system_f))^2 + 10 * (0. - last(energy_in_system_f))^2;
	
		# cost = sum(segments.diff_distance ./ input_speeds) + 100 * abs(minimum(energy_in_system_f))
		return cost
	end
	
	td = TwiceDifferentiable(f_wrap_reg, speeds; autodiff = :forward)
	lower_bound = fill(0.0, var_num)
	upper_bound = fill(upper_speed_bound, var_num)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)

	result = optimize(td, tdc, speeds 
	# .+ rand(vars_amount) .- 0.5
	    ,
	    IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10
	    )
	)

	minimized_speeds = Optim.minimizer(result)
	minimized_speeds_ms = minimized_speeds / 3.6

	power_use, solar_power, time_s = solar_trip_weather(
		minimized_speeds_ms, segments, start_datetime,
		weather_weights, weather_edges_lat, weather_edges_lon
	)
	# track points, not segments, that's why it is size is +1 
	energy_in_system_new = start_energy .+ solar_power .- power_use
	lowest_energy = minimum(energy_in_system_new)
	last_energy = last(energy_in_system_new)
	pushfirst!(energy_in_system_new, start_energy)

	track_plot = plot(track.distance, track.altitude, title="Track, $(var_num) segments",
		color=:green,
		ylabel="altitude(m)")

	speed_plot = plot(
		get_mean_data(track.distance),
		minimized_speeds,
		seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed plot for $(var_num) vars, total time $(round(last(time_s), digits=3)) sec",
		ylabel="speed(kmh)"
	)

	low_energy_red = fill(0., size(track.distance, 1))

	energy_plot = plot(
		track.distance,
		[energy_in_system_new low_energy_red],
		linewidth=[1 3],
		title="Energy for $(var_num) vars, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)"
	)
	println(time_s)
	println(energy_in_system_new)
	plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(650,700), legend=false)

end

# ╔═╡ 941b2f24-83a7-410c-aba3-51083c99c43e
function regular_optim_weather_from(track_input, segments_input, start_index, speeds, start_energy, start_datetime, 
	weather_weights, weather_edges_lat, weather_edges_lon,
	upper_speed_bound=150.
)
	orig_size = size(track_input.distance, 1)
	track = track_input[start_index:orig_size,:]
	segments = segments_input[start_index:orig_size-1,:]
	
	var_num = size(segments,1)
	
	function f_wrap_reg(input_speeds)
		speed_vector = input_speeds / 3.6;
		power_use_f, solar_power_f, time_s_f = solar_trip_weather(
			speed_vector, segments, start_datetime,
			weather_weights, weather_edges_lat, weather_edges_lon
		)
		energy_in_system_f = start_energy .+ solar_power_f .- power_use_f
		pushfirst!(energy_in_system_f, start_energy)
	
		# cost = sum(segments_clean.diff_distance ./ input_speeds) + 100 * (0. - last(energy_in_system_clean))^2;
		
		cost = sum(segments.diff_distance ./ speed_vector) + 1000 * abs(minimum(energy_in_system_f))^2 + 10 * (0. - last(energy_in_system_f))^2;
	
		# cost = sum(segments.diff_distance ./ input_speeds) + 100 * abs(minimum(energy_in_system_f))
		return cost
	end
	
	td = TwiceDifferentiable(f_wrap_reg, speeds; autodiff = :forward)
	lower_bound = fill(0.0, var_num)
	upper_bound = fill(upper_speed_bound, var_num)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)

	result = optimize(td, tdc, speeds 
	# .+ rand(vars_amount) .- 0.5
	    ,
	    IPNewton(),
	    Optim.Options(
	        x_tol = 1e-10,
	        f_tol = 1e-10,
	        g_tol = 1e-10
	    )
	)

	minimized_speeds = Optim.minimizer(result)
	minimized_speeds_ms = minimized_speeds / 3.6

	power_use, solar_power, time_s = solar_trip_weather(
		minimized_speeds_ms, segments, start_datetime,
		weather_weights, weather_edges_lat, weather_edges_lon
	)
	# track points, not segments, that's why it is size is +1 
	energy_in_system_new = start_energy .+ solar_power .- power_use
	lowest_energy = minimum(energy_in_system_new)
	last_energy = last(energy_in_system_new)
	pushfirst!(energy_in_system_new, start_energy)
	time_utc_res = travel_time_to_datetime(time_s, start_datetime)
	projected_finish_time = last(time_utc_res)
	# projected_finish_time = start_datetime + Dates.Millisecond(round(last(time_s) * 1000))

	point_times = start_datetime .+ Dates.Millisecond.(round.(time_s*1000))
	pushfirst!(point_times, start_datetime)

	for i in 1:var_num+1
		println("$i: $(point_times[i]), $(energy_in_system_new[i])")
	end

	track_plot = plot(track.distance, track.altitude, title="Track, $(var_num) segments",
		color=:green,
		ylabel="altitude(m)",
		xlimits=(-20, last(track.distance)+1000)
	)

	speed_plot = plot(
		get_mean_data(track.distance),
		minimized_speeds,
		seriestype=:bar,
		bar_width=segments.diff_distance,
		title="Speed from $(start_index) point, proj. finish at $(projected_finish_time)",
		ylabel="speed(kmh)",
		xlimits=(-20, last(track.distance)+1000)
	)

	low_energy_red = fill(0., size(track.distance, 1))

	energy_plot = plot(
		track.distance,
		[energy_in_system_new low_energy_red],
		linewidth=[1 3],
		title="Energy from $(start_index) point, lowest $(round(lowest_energy, digits=2)), last $(round(last_energy, digits=2))",
		xlabel="distance(m)", ylabel="energy(w*h)",
		xlimits=(-20, last(track.distance)+1000)
	)
	
	plot(track_plot, speed_plot, energy_plot, layout=(3,1), size=(700,700), legend=false)

end

# ╔═╡ 9147e247-e7bf-4f16-b1a9-a270bae7c7e4
md"## Простая австралия, опять"

# ╔═╡ 8b34704f-7d4f-4929-94cf-6fa6e550d6af
md"Ещё раз простая австралия с одной переменной"

# ╔═╡ 9d014721-b6dc-4e2e-9ec8-814f5ece3d6a
single_optim(
	track_aus,
	segments_aus,
	5100.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 89bffc0d-79af-4011-8867-3d85000978e1
md"Одна переменная с облаками"

# ╔═╡ d36c07b3-87fc-47a4-8b64-40ec984de403
single_optim_weather(
	track_aus,
	segments_aus,
	5100.,
	DateTime(2022,1,1,10,0,0),
	w_test2,
	edges_lat_test2,
	edges_lon_test2
)

# ╔═╡ 6efe9e33-2346-4176-bf77-cd89e7691f9e
md"Со всеми переменными"

# ╔═╡ d14b8b83-26ef-4530-b6b2-d700e674075d
regular_optim(
	track_aus,
	segments_aus,
	fill(50., size(segments_aus,1)),
	5100.,
	DateTime(2022,1,1,10,0,0),
	200.
)

# ╔═╡ 08f3f82d-f4b5-47da-a645-5da40cde61a5
md"Все переменные, с погодой"

# ╔═╡ 0350f096-b3cb-4141-9c58-c221748fd12f
regular_optim_weather(
	track_aus,
	segments_aus,
	fill(55., size(segments_aus,1)),
	5100.,
	DateTime(2022,1,1,10,0,0),
	w_test2,
	edges_lat_test2,
	edges_lon_test2,
	150.
)

# ╔═╡ 9dc3a519-3672-4956-8237-00d16da2ea26
md"Ну, это уже успех!

Почему? Потому что в самом проблемном месте (где тучка), соптимизировалась наибольшая скорость.

То есть надо как можно быстрее проезжать тучные места.

И стало чуточку быстрее.

Но у нас слишком много энрегии накопилось (надо что-то делать с ограничениями, или как-то скейлить штрафы в зависимости от размера задачи"

# ╔═╡ 23612b01-0bc3-4620-96f5-416895ce55cd
regular_optim_weather(
	track_aus,
	segments_aus,
	fill(55., size(segments_aus,1)),
	5100.,
	DateTime(2022,1,1,10,0,0),
	w_test2,
	edges_lat_test2,
	edges_lon_test2,
	150.
)

# ╔═╡ 6fb38907-fbec-4f2b-81b9-63adb1a1b5d2
md"Вот я поменял всего лишь максимальную скорость, и всё рассыпалось. Без черри-пикинга не работает"

# ╔═╡ c94bc74e-f68a-4522-a7fa-bdcfbd8ec413
md"Итого, что надо делать то?

вариант 1 - Переделать модель на JuMP, может станет лучше, т.к. правильней будут учитываться ограничения

вариант 2 - Попробовать посчитать это итеративно? То есть сужаем пространство поиска, т.к. оптимизатор застрявает в локальном минимуме

вариант 3 - попробовать backward diff - не факт что получится.

вариант 4 - тупая стратегия про одна скорость до холма?
"

# ╔═╡ ff6c92b6-783b-445e-a3b7-f70b689df4ad
md"# План №2 (исправленный)

План довольно прост. Надо доказать, что нам есть зачем выбирать разные скорости на разных участках трассы.

Шаги для осуществления плана:
1. Взять ровную трассу, подобрать одну скорость на все участки (отлично работает) +
2. Взять ровную трассу, подобрать разные скорости на все участки (работает, но избыточно) + 
3. Взять трассу с холмом, подобрать одну скорость на все участки (плохо работает, нужен другой подход) +
4. Взять трассу с холмом, подобрать разные скорости для каждого участка (должно норм работать) +
5. Взять трассу ощутимой длины (участков 300-500) с холмами, чтобы оно считалось достаточно долго. Подобрать разные скорости обычной оптимизацией и сказать что долго выходит (по идее норм посчитает, но долго) - (работает иногда нормально, иногда не сходится)
6. Взять эту же трассу с холмами, подобрать скорости моим методом (в идеале должно посчитать примерно так же, но быстрее) - до этого даже не дошл
7. Сделали статические осадки. Холм с осадками +
8. Холм с осадками разные скорости +, есть улучшение


ДОПИСАТЬ!!"

# ╔═╡ 2bb84981-520c-4e6c-9e19-568dda484881
md"# Проверяем накопление энергии"

# ╔═╡ f181c3e1-b5c4-407c-8003-b1e624e65761
md"Есть сомнение, что правильно считается накопление энергии в ночное время"

# ╔═╡ 5296b7de-93f9-415d-b268-a61eb4390dfa
Dates.dayofyear(DateTime(2022,5,1,0,0))

# ╔═╡ 2c95d568-eb5d-447d-a21e-786df523a1b4
begin
	# through 1 24-hour period, whole day (day+night)
	seconds = 0:60*60*25
	date = DateTime(2023,6,20)
	utc_times = date .+ Dates.Second.(seconds)
	# utc_times = DateTime(2023,1,1,0,0) .+ Dates.Second.(seconds)

	array_len = length(utc_times)

	lat = fill(-23.7, array_len)
	lon = fill(133.87, array_len)
	alt = fill(2000., array_len)

	lat_spb = fill(59.9386300, array_len)
	lon_spb = fill(30.3141300, array_len)
	alt_spb = fill(3., array_len)
	
	irradiance = solar_radiation_alloc.(
		lat,
		lon,
		alt,
		utc_times
	)
	irradiance_spb = solar_radiation_alloc.(
		lat_spb,
		lon_spb,
		alt_spb,
		utc_times
	)
	plot(
		utc_times,
		[irradiance irradiance_spb],
		title="Дневная солнечная радиация на $(Date(date))",
		# legend=false,
		legend=:right,
		labels=["Alice Springs, Австралия" "Санкт-Петербург, Россия"],
		ylimits=(-20,1300)
	)
end

# ╔═╡ 35a3b594-89c9-4d96-9eae-e5d20cca5476
md"Выглядит нормально"

# ╔═╡ 175f6fca-6c59-4814-9460-4dc9b1e2cf79
md"# Генерируем картинки"

# ╔═╡ 6ad29152-b8ad-47fd-bfa9-df87b64396d0
simulate_run(
	fill(10., size(segments_hill, 1)),
	track_hill,
	segments_hill,
	75.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 172784cf-2df6-4c90-910c-e6362847999a
simulate_run(
	fill(60., size(segments_hill, 1)),
	track_hill,
	segments_hill,
	75.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 12bbb3a5-e742-40c4-b89e-4e94574153cb
simulate_run(
	[55,45,60,40,30,40,50,55,60,40,65,35,45,60,55,60,50,45],
	track_aus,
	segments_aus,
	5100.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 7ba245e1-861c-47bd-b6b6-f7bd253cacb9
simulate_run(
	[109.17170893064288, 62.994584397266614, 108.33077015054835, 147.885468922807, 51.21062154781847, 98.67897603416483, 129.04537257465546, 133.39540229034336, 44.171273543445515, 41.24716187838297, 62.79237229568634, 1.8570579263272768, 51.58712986415076, 131.66717962274188, 104.66828180485733, 57.06231084905412, 89.68752269233823, 6.966200571105846],
	track_aus,
	segments_aus,
	5100.,
	DateTime(2022,1,1,10,0,0)
)

# ╔═╡ 14a6ee0b-ccef-4d6c-9008-cc67e8c486c9
println(rand(18) .* 150.)

# ╔═╡ 469732ca-7e68-43f2-a3f2-234c7fad425a
Random.seed!(1234)

# ╔═╡ 0859e307-6ab7-492a-8c57-802c44b1d814
rand(5)

# ╔═╡ 7e44936d-a96f-4fc5-b807-0bf41b298180
plot(
	[10,30,50,70,90],
	[10,20,15,25,10],
	seriestype=:bar,
	bar_width=fill(20, 5),
	title="Пример разбиения переменных",
	ylabel="Скорость(kmh)",
	xlabel="Номер участка",
	xticks=[0,20,40,60,80,100],
	legend=false
)

# ╔═╡ ff95e97c-360f-4a2a-bfa3-75fcf886aacd
md"## Картинки для плана энергосбережения"

# ╔═╡ e33c012d-ea0d-4a48-bc67-50dd45c9a967
plotjs.Plot(
	plotjs.scattermapbox(
		lat=track_aus.latitude,
		lon=track_aus.longitude,
		marker_color="red",
		marker_size=1,
		mode="lines"
	),
	plotjs.Layout(
		width=700,
		height=600,
		geo_fitbounds="locations",
		mapbox_style=mapbox_style,
		autosize=true,
		mapbox_center_lat=-25.0,
		mapbox_center_lon=132.0,
		mapbox_zoom=3,
		title="Полная дистанция"
	)
)

# ╔═╡ 53b029d2-7ad7-49be-a911-b9f6e54fc891
plotjs.Plot(
	plotjs.scattermapbox(
		lat=track_aus.latitude[8:end],
		lon=track_aus.longitude[8:end],
		marker_color="red",
		marker_size=1,
		mode="lines"
	),
	plotjs.Layout(
		width=700,
		height=600,
		geo_fitbounds="locations",
		mapbox_style=mapbox_style,
		autosize=true,
		mapbox_center_lat=-25.0,
		mapbox_center_lon=132.0,
		mapbox_zoom=3,
		title="Оставшаяся часть дистанции"
	)
)

# ╔═╡ 6c705457-ffd6-46d4-bc0a-6a6b60493fa9
function generate_plan(track, segments, start_index, start_datetime)
	# нужна карта, управляющие воздействия и ожидаемое время финиша

	traces_vector::AbstractVector{plotjs.AbstractTrace} = [];
	push!(
		traces_vector,
		plotjs.scattermapbox(
			lat=track.latitude[start_index:end],
			lon=track.longitude[start_index:end],
			marker_color="red",
			marker_size=1,
			mode="lines"
		)
	)
	map_plot = plotjs.Plot(
		plotjs.scattermapbox(
			lat=track.latitude[start_index:end],
			lon=track.longitude[start_index:end],
			marker_color="red",
			marker_size=1,
			mode="lines"
		),
		plotjs.Layout(
			width=700,
			height=600,
			geo_fitbounds="locations",
			mapbox_style=mapbox_style,
			autosize=true,
			mapbox_center_lat=-25.0,
			mapbox_center_lon=132.0,
			mapbox_zoom=3,
			title="Оставшаяся часть дистанции"
		)
	)
	track_plot = plotjs.Plot(
		track.distance[start_index:end],
		track.altitude[start_index:end],
		plotjs.Layout(
			width=600,
			height=300
		)
	)
	plots = [map_plot; track_plot]

	plots.layout["showlegend"] = false
    plots.layout["width"] = 700
    plots.layout["height"] = 600
    return plotjs.Plot(plots)

end

# ╔═╡ 478c720b-7c13-4990-b14b-d1e5c3f80953
generate_plan(track_aus, segments_aus, 1, DateTime(2023,1,1,10,0,0))

# ╔═╡ fd7dc6b5-100a-438b-a7d3-662b9ca32e66
plotjs.Plot(
	track_aus.distance,
	track_aus.altitude,
	plotjs.Layout(
		width=600,
		height=300
	)
)

# ╔═╡ b1e3829c-93b8-4546-9f4a-a0bcd0de120b
regular_optim_weather(
	track_aus,
	segments_aus,
	fill(55., size(segments_aus,1)),
	5100.,
	DateTime(2022,1,1,10,0,0),
	w_test2,
	edges_lat_test2,
	edges_lon_test2,
	200.
)

# ╔═╡ 17be5e6f-de9e-46a0-ade2-8d7f59d4e334
md"## План с момента"

# ╔═╡ 01f41a82-1783-4595-9fc0-d7507819bc85
regular_optim_weather_from(
	track_aus,
	segments_aus,
	1,
	fill(55., size(segments_aus,1)-1+1),
	5100.,
	DateTime(2022,1,1,10,0,0),
	w_test2,
	edges_lat_test2,
	edges_lon_test2,
	150.
)

# ╔═╡ 5f56b21f-8df4-472b-a7da-cb50e597aa99
regular_optim_weather_from(
	track_aus,
	segments_aus,
	8,
	fill(45., size(segments_aus,1)-8+1),
	126.451,
	DateTime(2022,1,2,15,14,59),
	w_test2,
	edges_lat_test2,
	edges_lon_test2,
	150.
)

# ╔═╡ bfdede64-c065-44ad-847f-59ec460b0591


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
TimeZones = "f269a46b-ccf7-5d73-abea-4c690281aa53"
WebIO = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"

[compat]
CSV = "~0.10.10"
DataFrames = "~1.5.0"
Distributions = "~0.25.90"
LineSearches = "~7.2.0"
Optim = "~1.7.5"
Peaks = "~0.4.3"
PlotlyJS = "~0.18.10"
Plots = "~1.38.11"
PlutoUI = "~0.7.51"
ProgressMeter = "~1.7.2"
StatsBase = "~0.33.21"
TimeZones = "~1.9.2"
WebIO = "~0.8.20"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "53f1e4e9f960959315ecfa6b208c19e1fd9ae3e1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38911c7737e123b28182d89027f4216cfc8a9da7"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.3"

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

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Blink]]
deps = ["Base64", "Distributed", "HTTP", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Pkg", "Reexport", "Sockets", "WebIO"]
git-tree-sha1 = "88616b94aa805689cf12f74b2509410135c00f43"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "ed28c86cbde3dc3f53cf76643c2e9bc11d56acc7"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.10"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

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
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"

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
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

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
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "eead66061583b6807652281c0fbf291d7a9dc497"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.90"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

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
git-tree-sha1 = "fc86b4fd3eff76c3ce4f5e96e2fdfa6282722885"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.0.0"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "6604e18a0220650dbbea7854938768f15955dd8e"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.20.0"

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
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

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

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

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
git-tree-sha1 = "efaac003187ccc71ace6c755b197284cd4811bfe"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.4"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4486ff47de4c18cb511a0da420efebb314556316"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.4+0"

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
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "877b7bc42729aa2c90bbbf5cb0d4294bd6d42e5a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "84204eae2dd237500835990bcade263e27674a93"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.16"

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
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "b413a73785b98474d8af24fd4c8a975e31df3658"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.4"

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

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43032da5832754f58d14a91ffbe86d5f176acda9"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.2.1+0"

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
git-tree-sha1 = "099e356f267354f46ba65087981a77da23a279b7"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.0"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "782e258e80d68a73d8c916e55f8ced1de00c2cea"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.6"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "87c371d27dbf2449a5685652ab322be163269df0"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.15"

[[deps.Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "MbedTLS", "Pkg", "Sockets"]
git-tree-sha1 = "0bdaa479939d2a1f85e2f93e38fbccfcb73175a5"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "1.0.1"

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

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "a89b11f0f354f06099e4001c151dffad7ebab015"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.5"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Peaks]]
deps = ["Compat", "RecipesBase"]
git-tree-sha1 = "ca47b866754525ede84e5dec84a104c45f92afb6"
uuid = "18e31ff7-3703-566c-8e60-38913d67486b"
version = "0.4.3"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

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

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "7452869933cd5af22f59557390674e8679ab2338"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.10"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "6c7f47fd112001fc95ea1569c2757dffd9e81328"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.11"

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
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

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
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "213579618ec1f42dea7dd637a42785a608b1ea9c"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

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

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

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
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

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
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "c262c8e978048c2b095be1672c9bee55b4619521"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.24"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

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

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

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

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Scratch", "Unicode"]
git-tree-sha1 = "a5404eddfee0cf451cabb8ea8846413323712e25"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.9.2"

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
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "976d0738247f155d0dcd77607edea644f069e1e9"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.20"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "4162e95e05e79922e44b9952ccbc262832e4ad07"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.6.0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

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
# ╠═94da1f00-efe7-11ed-2d94-f9a905085f40
# ╠═7282988b-0667-409d-9634-d874e7767d16
# ╠═c6a2816d-e81e-449e-af8c-93676c3fd077
# ╠═b20a1cbd-705e-4b49-b671-d042d1511afe
# ╠═d5213897-9cb6-453a-bc4e-a7015909c886
# ╠═769d7af1-e95c-426f-aa97-85b219c5b65e
# ╠═3e46a49b-7bb7-4889-8c94-b842977899e4
# ╠═64ecb4d9-b67c-4a41-a79e-300affd0440f
# ╠═9529a535-3655-496c-8077-30c148341fb9
# ╠═e3fa515b-10bd-4ae3-9c47-d8b5e365c4ea
# ╠═1029bac5-c3ae-4dd9-b992-db077cf69c18
# ╠═4b3430b1-619f-4123-9a5b-6db7e7083ee9
# ╠═018eeb89-8463-414d-badb-f64e1105ffdd
# ╠═2274f3b1-6a63-4206-a9ab-fb38c7a638f3
# ╠═df5b2134-2bd7-4e07-91ff-e5ec736f2dfa
# ╠═e9a6a82a-ecde-4baf-ad2c-f0a858bb2ed0
# ╠═9e47be95-5f9a-43d6-876e-abb33a67de0a
# ╠═9de8094a-5d24-4b08-8879-c7853cf6954b
# ╠═b965ae63-8e6b-4c0a-91b6-c29bae610590
# ╟─c892b78b-58dd-4a61-b731-6f4c1432ea1a
# ╠═c9daaf07-3dc2-4bef-ae23-b69acd2cb513
# ╟─a4ebc847-08fd-4976-87a7-2c35ca7920be
# ╟─8f1d25de-f4dd-481a-a7f5-b8691586e705
# ╟─46a01c40-c920-483c-aa93-1297eb0cddc8
# ╠═7f3807d7-e3ef-4c8f-8a21-8abdd4f25259
# ╠═9f151595-8f5a-4337-90fb-a8ef10f453b3
# ╠═67e01c26-832a-4eda-863d-e1842aae13dd
# ╟─074e1eeb-98db-4728-8bfe-4b4ed91b554f
# ╟─85d293fc-2992-4a38-a080-93179ff4a9f1
# ╠═b63f2571-079b-4ef3-9f96-6a465c255733
# ╟─48291326-f37d-4397-aeb0-fee0855970d8
# ╟─3509ebf4-b758-4ce3-a69f-fb20c16615ef
# ╠═72773944-4541-4220-8f0a-151ccd67bcaa
# ╟─d6921b5f-d9ae-4886-a40e-98f9ad052ed2
# ╟─49c58d5c-11e2-4aff-a408-cd4aabff7312
# ╟─53c9b844-5d1e-4444-abe2-ec16f7cb1ee9
# ╠═e6540edc-e971-48e9-8679-689a3ed73e7b
# ╠═0b1564a2-8129-43ae-8831-77e18959be67
# ╠═53e5f786-6e6a-4892-8252-aa3389682a80
# ╠═afbea56b-194a-407e-9f74-8a3c9d0cbbed
# ╠═db0b009d-1abd-42c0-bbad-b6ba53af5fb5
# ╠═f664a938-2e64-4e37-9404-578720b4bd57
# ╠═ec45bd26-942d-4180-81f5-a9db7e227c4f
# ╠═d82662d3-5533-49cf-967f-5bfa1fb65b13
# ╠═d907fa8d-f398-48ee-bd2d-e1fea7149fb1
# ╠═de68d838-4e10-49fa-a92b-a5524ef378d7
# ╠═fc78cf3e-c12c-43cb-9ce1-07665134b2cc
# ╠═eb3bd221-532a-4299-88bc-3a19ccb4c7e3
# ╠═0c1b0b93-23e0-422a-a6ce-30b1c6b1ab51
# ╠═f7a72d85-ec30-4406-9a9e-a148d2244cf2
# ╠═bb34bc3a-4232-4c52-8358-08c571566055
# ╠═cd45b86d-7231-418f-a43c-b78b4fa2d879
# ╠═69bf6ffc-8c79-4573-bfb2-d264c66a7040
# ╠═e410277c-87ef-4f18-89cb-7ce8817920b5
# ╠═6e74e3be-fb37-4bdd-a059-87b513f8382a
# ╠═08027dc4-0f87-4c31-8c30-703314a0c512
# ╠═b4e85535-d731-4e90-bf33-24ba6ee69aa4
# ╠═940a6c60-a5e5-4bab-940b-4a3b027ecc65
# ╠═7adbc40f-8f0f-4434-96da-3962f7107ce9
# ╠═52c423ab-248a-4398-af19-6f80acd6164f
# ╠═f8d9184d-e88e-430c-80e9-e2a660fce6ba
# ╠═aae20b83-351f-44ef-8cc3-79404167fdc4
# ╠═9c686f0c-831a-4d5c-b18c-d5a459fe871a
# ╠═0eec4f71-31a7-435d-9ed5-1c86c2751f5b
# ╠═3b0cd1f5-b2f2-4178-8f0e-156f2e6ffd10
# ╠═9ad47da6-2b71-4354-9bd5-7113f98139d7
# ╠═6bc90ca1-dc8d-41e2-a9df-b54ccac74f79
# ╠═57f47d65-5a3e-4c86-be15-587d94fbf677
# ╠═ac174dbd-11d5-4859-b8f5-751af94f9ac6
# ╠═e5088f3d-87f1-41e8-b8e6-e5a98b1c6e83
# ╠═64628a15-fb96-4e65-ac4d-d07d9dab7fcb
# ╠═5709b86c-6e68-4dac-b7d8-a880e9eca00d
# ╠═8ade2b34-98e9-495c-b9e0-faa655652cef
# ╠═6cab74b6-650e-4d6f-b6a9-0cbc0282f104
# ╠═5fca02b2-24ae-4081-98b0-e7f640e65405
# ╠═e4188c54-b620-48cd-b1ff-448ac4a15950
# ╠═7bf6b5ed-ed8b-42c9-8772-38f4e2b3a0f9
# ╠═7bd787da-997c-48bc-8c5b-f181148ac964
# ╠═8e8ca19d-a381-4055-aef7-f27292aac611
# ╠═75f0759b-285c-402b-afdf-243cf5d81bbe
# ╠═e7575b4b-49bf-4498-b74b-69bf9d1762cb
# ╠═45183999-58ff-4bdd-a4b0-225b1e50a5ee
# ╠═00973819-68b4-4d7e-8c2a-bebfce2a820a
# ╠═3498e807-a860-4f73-9d97-68c6f46d7efe
# ╠═b1449601-ce61-4097-8b32-345275f859b8
# ╠═aba59435-42d9-4ac8-b768-0f017591d788
# ╠═63afb4e4-dc19-4106-8bb5-b51d0881675d
# ╠═642a882f-de23-4414-a579-b61088e98bd9
# ╠═78327f2d-7cf2-49e3-98fa-d90aeb479c45
# ╠═760e7514-e3c4-4366-93fc-e1379f96c2ad
# ╠═e24f1d14-e0d9-45e8-a415-172e56e7cc9c
# ╠═a803b045-f7ce-4a14-b363-1c166e34fbe7
# ╠═e07dc08f-9cd8-46b6-8e92-062372b955a3
# ╠═39e09aba-c4c5-46b1-8ef9-63403b0225eb
# ╠═f7d29865-473d-47f5-a40c-08f48b531401
# ╠═525a59f1-e377-4bb4-a922-51fb940ad688
# ╠═1e9d17e9-d9fb-4b25-8233-4fc2fa6ce23e
# ╠═63bf3f77-6abe-4280-ab2f-520beb93a353
# ╠═ac2a21e4-491b-4b92-826c-5d4d74cf9ef1
# ╠═8d5be066-4625-4b60-a789-557ff2bd4660
# ╠═8ae2eaaa-4a63-4c01-b0f7-0050e5e41b04
# ╠═7cba5ab0-8844-454d-8f7d-1fe5552f201c
# ╠═6ae05363-c7b9-492d-91f0-ebc7935e40db
# ╠═b1a41a59-3842-468f-af9c-84387129ca45
# ╠═b184d70e-9cdf-4531-87fc-6b606c54a8e8
# ╠═52c09620-90cd-4e67-b84c-c7458d67c2c2
# ╠═036ba7a9-8695-4022-a97b-dec4b3a08921
# ╠═c93555bf-a560-4d52-9545-c18af2749faa
# ╠═bd36b76f-51bb-43c3-a94f-eb675fa3ac75
# ╠═941b2f24-83a7-410c-aba3-51083c99c43e
# ╠═9147e247-e7bf-4f16-b1a9-a270bae7c7e4
# ╠═8b34704f-7d4f-4929-94cf-6fa6e550d6af
# ╠═9d014721-b6dc-4e2e-9ec8-814f5ece3d6a
# ╠═89bffc0d-79af-4011-8867-3d85000978e1
# ╠═d36c07b3-87fc-47a4-8b64-40ec984de403
# ╠═6efe9e33-2346-4176-bf77-cd89e7691f9e
# ╠═d14b8b83-26ef-4530-b6b2-d700e674075d
# ╠═08f3f82d-f4b5-47da-a645-5da40cde61a5
# ╠═0350f096-b3cb-4141-9c58-c221748fd12f
# ╠═9dc3a519-3672-4956-8237-00d16da2ea26
# ╠═23612b01-0bc3-4620-96f5-416895ce55cd
# ╠═6fb38907-fbec-4f2b-81b9-63adb1a1b5d2
# ╠═c94bc74e-f68a-4522-a7fa-bdcfbd8ec413
# ╠═ff6c92b6-783b-445e-a3b7-f70b689df4ad
# ╠═2bb84981-520c-4e6c-9e19-568dda484881
# ╠═f181c3e1-b5c4-407c-8003-b1e624e65761
# ╠═5296b7de-93f9-415d-b268-a61eb4390dfa
# ╠═2c95d568-eb5d-447d-a21e-786df523a1b4
# ╠═35a3b594-89c9-4d96-9eae-e5d20cca5476
# ╠═175f6fca-6c59-4814-9460-4dc9b1e2cf79
# ╠═6ad29152-b8ad-47fd-bfa9-df87b64396d0
# ╠═172784cf-2df6-4c90-910c-e6362847999a
# ╠═12bbb3a5-e742-40c4-b89e-4e94574153cb
# ╠═7ba245e1-861c-47bd-b6b6-f7bd253cacb9
# ╠═14a6ee0b-ccef-4d6c-9008-cc67e8c486c9
# ╠═469732ca-7e68-43f2-a3f2-234c7fad425a
# ╠═0859e307-6ab7-492a-8c57-802c44b1d814
# ╠═7e44936d-a96f-4fc5-b807-0bf41b298180
# ╠═ff95e97c-360f-4a2a-bfa3-75fcf886aacd
# ╠═e33c012d-ea0d-4a48-bc67-50dd45c9a967
# ╠═53b029d2-7ad7-49be-a911-b9f6e54fc891
# ╠═6c705457-ffd6-46d4-bc0a-6a6b60493fa9
# ╠═478c720b-7c13-4990-b14b-d1e5c3f80953
# ╠═fd7dc6b5-100a-438b-a7d3-662b9ca32e66
# ╠═b1e3829c-93b8-4546-9f4a-a0bcd0de120b
# ╠═17be5e6f-de9e-46a0-ade2-8d7f59d4e334
# ╠═01f41a82-1783-4595-9fc0-d7507819bc85
# ╠═5f56b21f-8df4-472b-a7da-cb50e597aa99
# ╠═bfdede64-c065-44ad-847f-59ec460b0591
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
