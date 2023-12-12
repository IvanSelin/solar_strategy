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

# ╔═╡ 621782b4-0289-41b3-8702-02b96ab5a37d
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
	using NLSolversBase
	using Printf
end

# ╔═╡ 95d48aa1-d500-4ed0-85da-65184ccf3364
begin
	include("../src/energy_draw.jl")
	include("../src/time.jl")
	include("../src/solar_radiation.jl")
	include("../src/track.jl")
	include("../src/utils.jl")
	include("../src/strategy_calculation.jl")
	include("../src/weather.jl")
end

# ╔═╡ cef231a9-6584-440f-926c-6d33d7df1250
md"В этом ноутбуке мы будем получать результаты научной новизны"

# ╔═╡ 82888790-7eff-11ee-1ea8-69148277dca2
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

# ╔═╡ 651b3d97-232a-4d9a-8534-f35828db5380
PlutoUI.TableOfContents()

# ╔═╡ 8e4c579d-d9cf-47be-8987-6e35f7aff9ff
plotlyjs()

# ╔═╡ 309ee8a3-f136-46ab-bb05-9fae649ce940
begin
	track_full, segments_full = get_track_and_segments("../data/data_australia.csv")
	track_peaks, segments_peaks, points_peaks = keep_extremum_only_peaks_segments_with_points(track_full)
	plot(track_peaks.distance, track_peaks.altitude, title="Track extremum only data built w/ Peaks.jl")
end

# ╔═╡ a2d44ee3-a967-4c3e-8c28-782c58384870
begin
	dimensions=15;
	w, elat, elon = generate_clouds(
		-10,
		130,
		-18,
		135,
		-12.5,
		131.2,
		0.5,
		0.5,
		dimensions,
		0.75
	);
	weather_coeff = calculate_weather_weights_for_segments(
	    w,
	    elat,
	    elon,
	    segments_peaks
	);
	segments_peaks.weather_coeff = weather_coeff
end

# ╔═╡ 3fb0db03-1cbe-4711-9ff7-2ccf8b0ee08a
begin
	weather_full_coeff = calculate_weather_weights_for_segments(
		w,
		elat,
		elon,
		segments_full
	)
	segments_full.weather_coeff = weather_full_coeff
end

# ╔═╡ a48f9dc8-654b-4072-a8ea-27b0e812324e
md"# 1 - применение AD"

# ╔═╡ e5dfb032-07c2-4d2e-bdcb-370c6d9f7c82
md"В этом разделе мы будем сравнивать оптимизацию с применением автоматического дифференцирования, и без него

Что надо сделать:
 - Взять трек длинной где-то в 500 участков
 - Оптимизировать скорости nelder-mead
 - методом 1-го порядка (BFGS)
 - IPNewton
 - Сравнить результаты (время оптимизации)
 - повторить для других размерностей
"

# ╔═╡ 146ede7c-64dd-443e-922e-e5fddc788d55
@bind short_track_len confirm(NumberField(1:size(track_full, 1), default=2))

# ╔═╡ 5a7b2e00-3819-4d48-8165-af5f0ab69011
md"## N=$(short_track_len)"

# ╔═╡ e39e5375-4d8f-4e24-b0d9-97d743428f0f


# ╔═╡ 6e727497-bec2-4b77-93ff-2344512795aa
variables_num = short_track_len -1

# ╔═╡ bba824d1-3ee6-4ccb-8a0b-3ebf73a08099
mutable struct OptimResultsStub
	time_run::Real
	minimizer::Vector{Real}
	OptimResultsStub() = new(1e5, fill(1e-5, variables_num))
end

# ╔═╡ c5ef2c9b-fc8f-450b-9edf-e543c85723c2
track_short = track_peaks[1:+short_track_len,:];

# ╔═╡ 8f0207f2-f7b9-4e88-acfa-bf142e50e6a9
segments_short = segments_peaks[1:variables_num,:];

# ╔═╡ 6a3c26ea-5de0-49e4-9067-af429173aed4
max_energy = 5100.

# ╔═╡ 94f81fb6-92ad-4b92-99cc-8453fe5b1589
start_energy = last(track_short.distance)/last(track_peaks.distance)*max_energy

# ╔═╡ aef683c8-fd5b-479d-83b2-e047ba0f8713
finish_energy = 0.

# ╔═╡ 6a9ac47d-fee4-4609-b3ed-fc5598b88762
start_datetime = DateTime(2023,1,1,10,0,0)

# ╔═╡ 35950f4c-a4f3-4e94-9ca9-b3e3c1550d09
function f_wrap_regular(input_speeds :: Vector{<: Real})
	speed_vector :: Vector{<: Real} = convert_kmh_to_ms_typed(input_speeds)
	power_use, solar_power, time_s = solar_trip_boundaries_typed(
		speed_vector, segments_short, start_datetime
	)

	last_energy = last(solar_power) - last(power_use) + start_energy
	solar_power -= power_use
	min_penalty = abs(minimum(solar_power) + start_energy)

	cost = sum(segments_short.diff_distance ./ abs.(speed_vector)) + 150000 * min_penalty^2 + 10000 * (finish_energy - last_energy)^2;

	return cost
end

# ╔═╡ 7a30a132-9b75-4aef-9752-c5fcefbca922
init_speeds = fill(40., variables_num)

# ╔═╡ 6f54da25-f910-441e-8851-b4bf3e847a4e
begin
	lower = fill(0., variables_num)
	upper = fill(100., variables_num)
end

# ╔═╡ e2737024-4cb2-4eec-90a1-d1137bd5e487
md"## 0-order. Nelder Mead"

# ╔═╡ a82feccb-57c8-4ab8-a9de-ec6117f3a38c
nm_optimizer = NelderMead()

# ╔═╡ 8608e4d1-8ff6-4c57-80eb-45e8b721c4b3
begin
	if short_track_len < 155
		results_nm = optimize(f_wrap_regular, lower, upper, init_speeds, Fminbox(nm_optimizer), Optim.Options(iterations=10000))
	else
		results_nm = OptimResultsStub()
	end
end

# ╔═╡ 9df91593-d0dc-42ed-a88a-b22c92f91486
# results_nm_unconstrained = optimize(f_wrap_regular, init_speeds, nm_optimizer, Optim.Options(iterations=1000000))

begin
	if short_track_len < 505
		results_nm_unconstrained = optimize(f_wrap_regular, init_speeds, nm_optimizer, Optim.Options(iterations=1000000))
	else
		results_nm_unconstrained = OptimResultsStub()
	end
end

# ╔═╡ dfdf787e-5fd4-430f-875a-b538182b085e
minimized_speeds_nm = Optim.minimizer(results_nm)

# ╔═╡ fbd63615-69b1-4f1b-bfae-3a295c46c406
minimized_speeds_nm_unconstrained = Optim.minimizer(results_nm_unconstrained)

# ╔═╡ ea338a3a-620a-4971-80bf-cc47b576ca63
md"Sanity check"

# ╔═╡ 50bdbde7-ed0f-437f-a563-d20eb4314f60
simulate_run(minimized_speeds_nm, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 068f2adf-63da-4cd9-9c17-aefc92e5f3b0
md"That's sane!"

# ╔═╡ 84eaf793-abb4-4dd9-8708-2ccbf55806dc
simulate_run(minimized_speeds_nm_unconstrained, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 494a4eec-af31-47a3-b552-7faf371d81d2
md"## 1-order. (L-)BFGS"

# ╔═╡ 18cf4375-dff5-4912-a34b-d2b929475443
bfgs_optimizer = BFGS()

# ╔═╡ 08313d1e-7822-4139-9495-1641836716fe
begin
	od = OnceDifferentiable(f_wrap_regular, init_speeds; autodiff = :forward)
	odc = OnceDifferentiableConstraints(lower, upper)
end

# ╔═╡ 1b6df0b2-572b-4255-8dd5-ee42d1a9a08d
begin
	if short_track_len < 505
		results_bfgs_no_diff = optimize(f_wrap_regular, init_speeds, bfgs_optimizer, Optim.Options(iterations=1000))
	else
		results_bfgs_no_diff = OptimResultsStub()
	end
end

# ╔═╡ c58758ea-4fd4-42ce-8a82-0d1a09f4f07b
minimized_speeds_bfgs_no_diff = Optim.minimizer(results_bfgs_no_diff)

# ╔═╡ 4ab6e5b4-d16c-4016-a4d8-d2b7c8db5b9b
begin
	if short_track_len < 20
		results_bfgs = optimize(f_wrap_regular, lower, upper, init_speeds, Fminbox(bfgs_optimizer), Optim.Options(iterations=1000))
	else
		results_bfgs = OptimResultsStub()
	end
end

# ╔═╡ 49bf8d7f-605d-4dbd-8077-8238573aea92
minimized_speeds_bfgs = Optim.minimizer(results_bfgs)

# ╔═╡ 4ccac8b1-070c-43cd-aeaf-69955d3ce2c3
results_bfgs_diff = optimize(f_wrap_regular, lower, upper, init_speeds, Fminbox(bfgs_optimizer), Optim.Options(iterations=10000); autodiff= :forward)

# ╔═╡ 54ff073c-c751-4a13-b6be-05d744537291
minimized_speeds_bfgs_diff = Optim.minimizer(results_bfgs_diff)

# ╔═╡ f405aaa5-e520-4d0c-bd63-ea4276a1be93
# results_bfgs_od = optimize(od, odc, init_speeds, BFGS(), Optim.Options(iterations=10000)) # not working

# ╔═╡ 978810c4-8f70-47c3-b1e3-0608ca855b1a
results_bfgs_diff_unconstrained = optimize(f_wrap_regular, init_speeds, bfgs_optimizer, Optim.Options(iterations=10000); autodiff= :forward)

# ╔═╡ 2039e2ff-a56d-4a44-89ec-9073ae98480b
minimized_speeds_bfgs_diff_unconstrained = Optim.minimizer(results_bfgs_diff_unconstrained)

# ╔═╡ f04149e4-a1f3-4a1d-bc3c-7c4e520ba976
simulate_run(minimized_speeds_bfgs_no_diff, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 2106d6cb-9a70-4ccb-907d-3f157cfe1c46
simulate_run(minimized_speeds_bfgs, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 82deaf73-4205-4b43-8e7e-4f6500012834
simulate_run(minimized_speeds_bfgs_diff, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 8bab99e6-b6de-4fa1-b5ae-0a31af1cf59a
simulate_run(minimized_speeds_bfgs_diff_unconstrained, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ cedabc1e-0cb6-4920-9bed-200f85b861c7
md"## 1-order. Conjugate Gradient"

# ╔═╡ e0a006aa-cb6b-474d-91af-be792a03deaa
cg_optimizer = ConjugateGradient()

# ╔═╡ a77ac935-1c8f-43af-b0ff-62bb39204562
results_cg = optimize(f_wrap_regular, lower, upper, init_speeds, Fminbox(cg_optimizer), Optim.Options(iterations=10000))

# ╔═╡ ebb6de9a-eedb-4352-af99-ba2f95c057d3
minimized_speeds_cg = Optim.minimizer(results_cg)

# ╔═╡ cefe8ebf-cbec-40b1-a906-b4ca714538ae
results_cg_diff = optimize(f_wrap_regular, lower, upper, init_speeds, Fminbox(cg_optimizer), Optim.Options(iterations=10000); autodiff = :forward)

# ╔═╡ 18b859fa-b004-4406-bbf3-81249e85be40
minimized_speeds_cg_diff = Optim.minimizer(results_cg_diff)

# ╔═╡ aad06b97-aade-47ae-8196-58bb8a0f9c75
results_cg_diff_unconstrained = optimize(f_wrap_regular, init_speeds, cg_optimizer, Optim.Options(iterations=10000); autodiff= :forward)

# ╔═╡ 70c353c1-c103-45ee-9ab0-b355326f96a6
minimized_speeds_cg_diff_unconstrained = Optim.minimizer(results_cg_diff_unconstrained)

# ╔═╡ 18d5923d-eb7f-48d1-80ea-5d147d99d8a5
results_cg_no_diff = optimize(f_wrap_regular, init_speeds, cg_optimizer, Optim.Options(iterations=10000))

# ╔═╡ 190d79b4-fcbc-4473-8e30-8316f2cded96
minimized_speeds_cg_no_diff = Optim.minimizer(results_cg_no_diff)

# ╔═╡ 3b499e5e-5235-4133-b35d-123f2ca18ada
simulate_run(minimized_speeds_cg, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 08493508-869b-46a2-aeb3-098b1ca99380
simulate_run(minimized_speeds_cg_diff, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 674597af-873d-44f6-bab7-8960bfa19ee3
simulate_run(minimized_speeds_cg_diff_unconstrained, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 4f298aa8-de16-4414-a794-e5a532c733b5
simulate_run(minimized_speeds_cg_no_diff, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ a5688643-48a7-459d-8ee4-f32e3e86211c
md"## 2-order. IPNewton"

# ╔═╡ ab2f7e0b-8aa9-4ddd-aec5-16c3e5de70b3
ipnewton_optimizer = IPNewton()

# ╔═╡ 8295f55e-69ec-4a17-8cb6-5e5441a632ce
begin
	td = TwiceDifferentiable(f_wrap_regular, init_speeds; autodiff = :forward)
	lower_bound = fill(0.0, variables_num)
	upper_bound = fill(100.0, variables_num)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
end

# ╔═╡ c0b97d5e-4a54-41ad-83d3-0d8923290d70
results_ipnewton = optimize(td, tdc, init_speeds 
# .+ rand(vars_amount) .- 0.5
    ,
    IPNewton(),
    Optim.Options(
        x_tol = 1e-6,
        f_tol = 1e-6,
        g_tol = 1e-6,
		allow_f_increases = true
    )
)

# ╔═╡ a477f042-55ac-4a8c-9c7d-fb3f6a85bc1a
minimized_speeds_ipnewton = Optim.minimizer(results_ipnewton)

# ╔═╡ b0685b3a-42eb-4312-86a2-c9ab9432b565
tdc_unconstrained = TwiceDifferentiableConstraints(
	fill(-Inf, variables_num),
	fill(Inf, variables_num)
)

# ╔═╡ 3d0e6774-90ad-4026-abf1-6feaa7579f30
clear!(td)

# ╔═╡ da104d8b-67f8-406f-b6b7-095a40596df8
results_ipnewton_unconstrained = optimize(td, tdc_unconstrained, init_speeds, ipnewton_optimizer, Optim.Options(iterations=10000))

# ╔═╡ 95d3b5a0-a4e6-4cc1-8b23-48e1d18961cf
minimized_speeds_ipnewton_unconstrained = Optim.minimizer(results_ipnewton_unconstrained)

# ╔═╡ 9d5fc300-afdc-4389-be3e-2f548fcafa76
simulate_run(minimized_speeds_ipnewton, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ e2409a9e-580b-4fc4-bb1a-3d07e53dceb4
simulate_run(minimized_speeds_ipnewton_unconstrained, track_short, segments_short, start_energy, start_datetime)

# ╔═╡ 530ec0eb-87d9-408b-9ca0-4f76ba29307d
md"## Результаты"

# ╔═╡ 7c7f9b0b-3648-4ae7-8686-0546294337fe
typeof(5)

# ╔═╡ 5b45416b-ab8e-4127-a637-99dbd76cef98
typeof(OptimResultsStub())

# ╔═╡ 439882df-7137-4157-89cf-d8ac92a6a8e6
typeof(results_nm) <: OptimResultsStub

# ╔═╡ ad5d9933-b68d-4eb4-816f-b083306189ef
function sim_time(optim_result)
	# if typeof(results_nm) <: OptimResultsStub
	# 	speeds = optim
	# speeds = Optim.minimizer(optim_result)
	speeds = optim_result.minimizer
	return sum(segments_short.diff_distance ./ (speeds ./ 3.6))
end

# ╔═╡ bf28df12-d9bf-4310-951a-bd22eab4b12b
function get_optim_res_string(optim_result)
	return "$(@sprintf("%10.3f",optim_result.time_run)) c, $(@sprintf("%10.3f",sim_time(optim_result))) c"#, $(optim_result.f_converged)"
end

# ╔═╡ d50301b9-de09-4b91-9104-af1a149ca6c0
function get_optim_res_string_excel(optim_result)
	return [@sprintf("%.3f",optim_result.time_run),@sprintf("%.3f",sim_time(optim_result))]
end

# ╔═╡ 690226f5-a3df-4aad-8cc7-c1f16c0ab033
get_optim_res_string_excel(results_nm)

# ╔═╡ c6798503-84ac-4d29-b1f2-f3a2037ced11


# ╔═╡ 89203abf-3941-41f1-8361-aae96f867c58
# print("N=$(short_track_len)
# run_time, sim_time 

# Nelder Mead Fminbox:
# $(get_optim_res_string(results_nm))

# Nelder Mead:
# $(get_optim_res_string(results_nm_unconstrained))

# BFGS Fminbox:
# $(get_optim_res_string(results_bfgs))

# BFGS Fminbox autodiff:
# $(get_optim_res_string(results_bfgs_diff))

# BFGS autodiff:
# $(get_optim_res_string(results_bfgs_diff_unconstrained))

# Conjugate Gradient Fminbox:
# $(get_optim_res_string(results_cg))

# Conjugate Gradient Fminbox autodiff:
# $(get_optim_res_string(results_cg_diff))

# Conjugate Gradient autodiff:
# $(get_optim_res_string(results_cg_diff_unconstrained))

# IPNewton:
# $(get_optim_res_string(results_ipnewton))

# IPNewton unconstrained:
# $(get_optim_res_string(results_ipnewton_unconstrained))
# ")

# ╔═╡ 3680d006-9322-41b9-93c6-51498d56889f
function generate_excel_string(optim_results_list)
	elems = []
	for optim_result in optim_results_list
		append!(elems, get_optim_res_string_excel(optim_result))
	end
	return join(elems,",")
end

# ╔═╡ fde27fab-aaa2-45bf-8c71-e92bb53dc25c
md"## Строка для копирования"

# ╔═╡ 6fcc5460-c724-4f6e-bb04-8785cddb7f79
generate_excel_string([results_nm, results_nm_unconstrained, results_bfgs_no_diff, results_bfgs, results_bfgs_diff, results_bfgs_diff_unconstrained, results_cg_no_diff,results_cg, results_cg_diff, results_cg_diff_unconstrained, results_ipnewton, results_ipnewton_unconstrained])

# ╔═╡ 23906645-c779-4eb1-9201-6302fec045ea
md"Это не способ повышения производительности, а повышения применимости?" 

# ╔═╡ 9ffbda77-5314-493e-89f6-bcaaff48a737
md"# 2 - Сокращение трассы"

# ╔═╡ 62499f26-9863-47f4-9aab-37d2f6f1b255
plot(
	track_peaks.distance / 1000.,
	track_peaks.altitude,
	title="Машрут",
	xlabel="Дистанция (км)",
	ylabel="Высота над уровнем моря (м)",
	legend=false
	# margin=10Plots.px,
)

# ╔═╡ 35b780d1-bf92-42bc-ac72-0a95e2585f56
md"Здесь мы будем убеждать всех, что мы хорошо сокращаем трассу"

# ╔═╡ 955df201-d9d4-4da2-b604-16619c92dd85
md"В принципе, уже всё сделано в ноутбуке про упрощение трассы

Надо только прирост скорости оптимизации сделать" 

# ╔═╡ 127ca813-a652-4cc7-9995-588a0ab863cb
md"Надо взять трассу до какого-то момента в peaks и параметрического (желательно одну и ту же точку)

И взять оригинальную трассу

Измерить время оптимизации

Затем взять ещё несколько таких точек"

# ╔═╡ b53eec1b-c644-4ef9-8086-6d90eb89935b
md"Ещё нужно сравнение на разных трассах для параметрического (желательно)"

# ╔═╡ f497ffca-5630-4d30-811f-bd1e15e8d34b
k = 1.75

# ╔═╡ 1246d429-eb90-41df-bb14-df94f67a152f
track_k, points_k = parametrized_track_simplification(track_full, k);

# ╔═╡ 83934f3d-bf9c-4cd2-976c-0d50a91cbbb8
segments_k = get_segments_for_track(track_k);

# ╔═╡ d7ab2d22-c546-4b6c-b2d8-7aad0af65261
md"Выясним, где одинаковые точки у peaks и k"

# ╔═╡ e269d75f-bc01-4621-b193-555f54746ead
intersect_points = intersect(points_peaks, points_k)

# ╔═╡ 93b53e5d-69e1-47d5-a6bc-4d1756467b2d
length(points_k)

# ╔═╡ 990ac9a0-fc6f-4afb-aed2-e679ddcadfbf
intersect_df = DataFrame((intersect=intersect_points, distance=track_full[intersect_points, :].distance))

# ╔═╡ 429fcd9c-300f-437b-ae15-16ff5040f3d3
peaks_bit_vector = in.(track_peaks.distance, Ref(intersect_df.distance))

# ╔═╡ b11ef547-15b2-4425-bc83-fa6b8071e0a9
peaks_indexes = findall(peaks_bit_vector)

# ╔═╡ bf08d4cd-d34b-4a42-82df-edbbb8d02437
intersect_df.peaks = peaks_indexes;

# ╔═╡ 0c17cb82-91e8-423a-8764-87b768e71a46
k_bit_vector = in.(track_k.distance, Ref(intersect_df.distance))

# ╔═╡ 251b48a8-3a82-4859-a1b7-95306e484375
k_indexes = findall(k_bit_vector)

# ╔═╡ 4807e032-e508-4795-9ec8-e830d9abfcc5
intersect_df.k = k_indexes;

# ╔═╡ cf34f9ff-744a-498b-9835-b478a2c2d6c5
md"Выберем некоторые из точек"

# ╔═╡ 9d1b383c-5dee-45b7-bfce-f1f4bd9ef1e1
# poi = [28, 54, 109, 202, 400, 518, 749, 1007, 1490, 2065, 2534] # k=1.0
poi = [5,28,54,109, 150, 281, 327, 518, 742, 1207, 1560, 2046, 2619] # k=1.75
# poi = [109, 150, 281]

# ╔═╡ 468f734c-5f23-420a-8b87-e068e44cc977
poi_df = intersect_df[in.(intersect_df.intersect, Ref(poi)), :]

# ╔═╡ db538803-f8c5-48f3-9051-adcfa1c4160d
md"sanity check"

# ╔═╡ cc8085a3-8ce6-45b0-a353-2bc3eddc4381
track_peaks.distance[42]

# ╔═╡ 9482aafe-37b6-482d-ae11-a9a64b183bfd
track_k.distance[22]

# ╔═╡ 63989b8f-fa2f-4151-b2cf-394e79b80474
md"sanity check passed

Теперь можно и посравнивать)"

# ╔═╡ 362f6a0b-b2e2-4bf0-a74d-59c1ec38578a
md"## Сравнение по выбранным точкам"

# ╔═╡ 09032708-dbde-43b1-9d18-70a881d65b7d
begin
	for row in eachrow(poi_df)
		println(row.intersect)
	end
end

# ╔═╡ 215e45a7-e1f3-45d5-af34-ec7e89fb96fb
begin
	track_optimization_comparison_df = DataFrame(
		points_regular=Int[], points_peaks=Int[], points_k=Int[],
		opt_time_regular=Float64[], opt_time_peaks=Float64[], opt_time_k=Float64[],
		res_time_regular=Float64[], res_time_peaks=Float64[], res_time_k=Float64[],
		min_energy_regular=Float64[], min_energy_peaks=Float64[], min_energy_k=Float64[],
		finish_energy_regular=Float64[], finish_energy_peaks=Float64[], finish_energy_k=Float64[],
		regular_speeds = Vector{Float64}[], peaks_speeds = Vector{Float64}[], k_speeds = Vector{Float64}[]
	)
	for row in eachrow(poi_df)
		# println(row)
		track_short_comparison = track_full[1:row.intersect,:]
		segments_short_comparison = get_segments_for_track(track_short_comparison);
		track_peaks_comparison = track_peaks[1:row.peaks,:]
		segments_peaks_comparison = get_segments_for_track(track_peaks_comparison);
		track_k_comparison = track_k[1:row.k,:]
		segments_k_comparison = get_segments_for_track(track_k_comparison);
		# println("$(size(segments_short_comparison,1)), $(size(segments_peaks_comparison,1)), $(size(segments_k_comparison,1))")

		start_energy_comparison = row.distance/last(track_full.distance) * max_energy
		init_speeds_regular = fill(40. , size(segments_short_comparison, 1) )
		init_speeds_peaks = fill(40. , size(segments_peaks_comparison, 1) )
		init_speeds_k = fill(40. , size(segments_k_comparison, 1) )
		# now it's time to optimize

		function f_wrap_regular_comparison(input_speeds :: Vector{<: Real})
			speed_vector :: Vector{<: Real} = convert_kmh_to_ms_typed(input_speeds)
			power_use, solar_power, time_s = solar_trip_boundaries_typed(
				speed_vector, segments_short_comparison, start_datetime
			)
		
			last_energy = last(solar_power) - last(power_use) + start_energy_comparison
			solar_power -= power_use
			min_penalty = abs(minimum(solar_power) + start_energy_comparison)
		
			cost = sum(segments_short_comparison.diff_distance ./ abs.(speed_vector)) + 150000 * min_penalty^2 + 10000 * (finish_energy - last_energy)^2;
		
			return cost
		end

		function f_wrap_peaks_comparison(input_speeds :: Vector{<: Real})
			speed_vector :: Vector{<: Real} = convert_kmh_to_ms_typed(input_speeds)
			power_use, solar_power, time_s = solar_trip_boundaries_typed(
				speed_vector, segments_peaks_comparison, start_datetime
			)
		
			last_energy = last(solar_power) - last(power_use) + start_energy_comparison
			solar_power -= power_use
			min_penalty = abs(minimum(solar_power) + start_energy_comparison)
		
			cost = sum(segments_peaks_comparison.diff_distance ./ abs.(speed_vector)) + 150000 * min_penalty^2 + 10000 * (finish_energy - last_energy)^2;
		
			return cost
		end

		function f_wrap_k_comparison(input_speeds :: Vector{<: Real})
			speed_vector :: Vector{<: Real} = convert_kmh_to_ms_typed(input_speeds)
			power_use, solar_power, time_s = solar_trip_boundaries_typed(
				speed_vector, segments_k_comparison, start_datetime
			)
		
			last_energy = last(solar_power) - last(power_use) + start_energy_comparison
			solar_power -= power_use
			min_penalty = abs(minimum(solar_power) + start_energy_comparison)
		
			cost = sum(segments_k_comparison.diff_distance ./ abs.(speed_vector)) + 150000 * min_penalty^2 + 10000 * (finish_energy - last_energy)^2;
		
			return cost
		end

		# или может лучше использовать NelderMead? чтобы больше была видна разница
		results_regular = optimize(f_wrap_regular_comparison, init_speeds_regular,
			NelderMead(), Optim.Options(iterations=10000))

		results_peaks = optimize(f_wrap_peaks_comparison, init_speeds_peaks, 
			NelderMead(), Optim.Options(iterations=10000))

		results_k = optimize(f_wrap_k_comparison, init_speeds_k, 
			NelderMead(), Optim.Options(iterations=10000))

		# simulate energies and run time
		regular_speeds = results_regular.minimizer ./ 3.6
		peaks_speeds = results_peaks.minimizer ./ 3.6
		k_speeds = results_k.minimizer ./ 3.6

		pu_regular, sp_regular, time_regular = solar_trip_boundaries_typed(
			regular_speeds, segments_short_comparison, start_datetime
		)
		pu_peaks, sp_peaks, time_peaks = solar_trip_boundaries_typed(
			peaks_speeds, segments_peaks_comparison, start_datetime
		)
		pu_k, sp_k, time_k = solar_trip_boundaries_typed(
			k_speeds, segments_k_comparison, start_datetime
		)

		time_sum_regular = last(time_regular)
		time_sum_peaks = last(time_peaks)
		time_sum_k = last(time_k)
		
		energy_regular = start_energy_comparison .+ sp_regular .- pu_regular
		low_e_regular = minimum(energy_regular)

		energy_peaks= start_energy_comparison .+ sp_peaks .- pu_peaks
		low_e_peaks = minimum(energy_peaks)

		energy_k = start_energy_comparison .+ sp_k .- pu_k
		low_e_k = minimum(energy_k)

		
		# minimized_speeds_ms = speeds / 3.6
	
		# power_use, solar_power, time_s = solar_trip_boundaries(
		# 	minimized_speeds_ms, segments, start_datetime
		# )
		# # track points, not segments, that's why it is size is +1 
		# energy_in_system_new = start_energy .+ solar_power .- power_use
		# lowest_energy = minimum(energy_in_system_new)
		# last_energy = last(energy_in_system_new)


		push!(
			track_optimization_comparison_df,
			(
				row.intersect, row.peaks, row.k,
				results_regular.time_run, results_peaks.time_run, results_k.time_run,
				time_sum_regular, time_sum_peaks, time_sum_k,
				low_e_regular, low_e_peaks, low_e_k,
				last(energy_regular), last(energy_peaks), last(energy_k),
				results_regular.minimizer, results_peaks.minimizer, results_k.minimizer
			)
		)
		# push!(res_df, (thr, last_diff, mae_val, mse_val, rmse_val, r2_val, number_of_segments, exec_time))
		
		# [@sprintf("%.3f",optim_result.time_run),@sprintf("%.3f",sim_time(optim_result))]
		
		
	end
	track_optimization_comparison_df
end

# ╔═╡ 88694580-4e7f-4e35-8c33-e90bf5679c7e
simulate_run(
	track_optimization_comparison_df.regular_speeds[6],
	track_full[1:poi_df.intersect[6],:],
	get_segments_for_track(track_full[1:poi_df.intersect[6],:]),
	poi_df.distance[6]/last(track_full.distance) * max_energy,
	start_datetime
)
# regular for 1st poi

# ╔═╡ 480e5187-19d7-486e-ad08-3a340bd7a760
simulate_run(
	track_optimization_comparison_df.peaks_speeds[6],
	track_peaks[1:poi_df.peaks[6],:],
	get_segments_for_track(track_peaks[1:poi_df.peaks[6],:]),
	poi_df.distance[6]/last(track_full.distance) * max_energy,
	start_datetime
)
# peaks for 1st poi

# ╔═╡ e3f80864-8fd2-48bd-acb8-8c91aba28e76
simulate_run(
	track_optimization_comparison_df.k_speeds[8],
	track_k[1:poi_df.k[8],:],
	get_segments_for_track(track_k[1:poi_df.k[8],:]),
	poi_df.distance[8]/last(track_full.distance) * max_energy,
	start_datetime
)
# k for 1st poi

# ╔═╡ 985db9d8-b795-45c7-991f-8b68d81feee9
track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_k

# ╔═╡ b32881c8-59c1-4806-91e9-f6447f2aa22b
plot(
	track_optimization_comparison_df.points_regular,
	track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_k,
	label="Прирост производительности",
	xlabel="Участков",
	legend=:topleft
)

# ╔═╡ d1146d90-664a-4235-95b0-8570e70713fc
plot(
	track_optimization_comparison_df.points_regular,
	track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_k,
	ylabel="Performance gain",
	xlabel="Number of segments",
	label="k=1.75",
	title="k=1.75 performance",
	# legend=:topleft,
	legend=false,
	size=(400,300),
	# marginbottom=5Plots.mm
	bottom_margin = 10Plots.px
	# , layout=(2,1)
)

# ╔═╡ f8780df8-5868-49da-a5a8-a556afdb3f18
plot(
	track_optimization_comparison_df.points_regular,
	track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_peaks,
	ylabel="Performance gain",
	xlabel="Number of segments",
	label="Peaks",
	title="Peaks performance",
	# legend=:topleft,
	legend=false,
	size=(400,300),
	# marginbottom=5Plots.mm
	bottom_margin = 10Plots.px
	# , layout=(2,1)
)

# ╔═╡ 75a0a6e3-4462-4f1d-baa6-3364d211a4b5
plot(
	track_optimization_comparison_df.points_regular,
	[
		track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_peaks track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_k
	],
	ylabel="Performance gain",
	xlabel="Number of segments",
	labels=["Extremum points" "Parametric, k=1.75"],
	# title="Peaks performance",
	markers=[:diamond :circle],
	legend=:topleft,
	# legend=false,
	# size=(400,300),
	# marginbottom=5Plots.mm
	bottom_margin = 10Plots.px
	# , layout=(2,1)
)

# ╔═╡ a0775a71-0cca-43a6-835a-d66850605a45
good_peaks_index = abs.(track_optimization_comparison_df.min_energy_peaks) .< 5.0

# ╔═╡ fd6a8bb0-7691-4a03-a928-999ddf6e3076
good_k_index = abs.(track_optimization_comparison_df.min_energy_k) .< 5.0

# ╔═╡ f9d7c72f-c666-4a4a-b31a-c92b37507da9
plot(
	track_optimization_comparison_df.points_regular[good_peaks_index],
	track_optimization_comparison_df.opt_time_regular[good_peaks_index] ./ track_optimization_comparison_df.opt_time_peaks[good_peaks_index] ,
	ylabel="Performance gain",
	xlabel="Number of segments",
	labels=["Extremum points" ""],
	# title="Peaks performance",
	marker=:diamond,
	legend=:topleft,
	# legend=false,
	# size=(400,300),
	# marginbottom=5Plots.mm
	bottom_margin = 10Plots.px
	# , layout=(2,1)
)

# ╔═╡ 62043829-4ca6-4dd7-bf4f-85c7343ee177
plot!(
	track_optimization_comparison_df.points_regular[good_k_index],
	track_optimization_comparison_df.opt_time_regular[good_k_index] ./ track_optimization_comparison_df.opt_time_k[good_k_index] ,
	ylabel="Performance gain",
	xlabel="Number of segments",
	labels=["k=$(k)" "2"],
	# title="Peaks performance",
	marker=:circle,
	legend=:topleft,
	# legend=false,
	# size=(400,300),
	# marginbottom=5Plots.mm
	bottom_margin = 10Plots.px
	# , layout=(2,1)
)

# ╔═╡ 46c90566-c373-40c1-b4ca-1ae1c9872a87
good_orig_index = abs.(track_optimization_comparison_df.min_energy_regular) .< 5.0

# ╔═╡ dcdf64c2-01df-4d73-a365-1b6381c9c3e1
until_num = findlast(good_orig_index)

# ╔═╡ 3a6eb336-4ec3-400a-b58b-418c7429ca2d
sum((track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_k)[1:until_num]) / until_num

# ╔═╡ 2f0d3995-cbfd-420f-8936-380fc5d5da88
sum((track_optimization_comparison_df.opt_time_regular ./ track_optimization_comparison_df.opt_time_peaks)[1:until_num]) / until_num

# ╔═╡ a694736a-a239-4b3c-8670-9982aaf0204e
md"### Сравнение всех трёх "

# ╔═╡ 714c393b-bb42-4622-848f-1dca1868169d
function draw_comparison_plot(index)

	reg_plot = simulate_run(
		track_optimization_comparison_df.regular_speeds[index],
		track_full[1:poi_df.intersect[index],:],
		get_segments_for_track(track_full[1:poi_df.intersect[index],:]),
		poi_df.distance[index]/last(track_full.distance) * max_energy,
		start_datetime
	)


	peaks_plot = simulate_run(
		track_optimization_comparison_df.peaks_speeds[index],
		track_peaks[1:poi_df.peaks[index],:],
		get_segments_for_track(track_peaks[1:poi_df.peaks[index],:]),
		poi_df.distance[index]/last(track_full.distance) * max_energy,
		start_datetime
	)

	k_plot = simulate_run(
		track_optimization_comparison_df.k_speeds[index],
		track_k[1:poi_df.k[index],:],
		get_segments_for_track(track_k[1:poi_df.k[index],:]),
		poi_df.distance[index]/last(track_full.distance) * max_energy,
		start_datetime
	)

	plot(
		reg_plot, peaks_plot, k_plot,
		layout=(3,1),
		size=(750,1500),
		legend=false
	)
	
end

# ╔═╡ a13bb6c3-0bed-4c2f-868f-5acded90aada
draw_comparison_plot(5)

# ╔═╡ 75983ca9-e54a-419d-895c-892916135484
md"Сокращённые управляющие воздействия надо раскрывать на всю трассу, а не гонять по сокращённой"

# ╔═╡ 75ec8028-db00-4e56-936c-dddf5c286735
points_peaks[1:poi_df.peaks[5]]

# ╔═╡ d105fc37-e079-4c5d-aa6b-288e4828c993
function peaks_speeds_to_regular(index)
	peaks_speeds = track_optimization_comparison_df.peaks_speeds[index]
	# peaks_track = track_peaks[1:poi_df.peaks[index],:]
	# peaks_segments = get_segments_for_track(track_full[1:poi_df.intersect[index],:])
	# start_energy = poi_df.distance[index]/last(track_full.distance) * max_energy

	# orig_track = track_full[1:poi_df.intersect[index],:]

	points_to_analyze = points_peaks[1:poi_df.peaks[index]]
	new_speed_vector_for_regular = []
	for i=2:size(points_to_analyze,1)
		append!(
			new_speed_vector_for_regular,
			fill(
				peaks_speeds[i-1],
				points_to_analyze[i] - points_to_analyze[i-1]
			)
		)
	end
	# println("new size $(size(new_speed_vector_for_regular,1))")
	# println("old size $(poi_df.intersect[index]-1)")
	# @assert size(new_speed_vector_for_regular,1)==poi_df.intersect-1
	return new_speed_vector_for_regular
end

# ╔═╡ 2e35d429-cbb1-4e63-a03b-eaa9662229ab
peaks_speeds_to_regular(2)

# ╔═╡ 2d8a9498-aa45-46d8-947e-94da8a07c27f
function k_speeds_to_regular(index)
	k_speeds = track_optimization_comparison_df.k_speeds[index]
	# k_track = track_peaks[1:poi_df.k[index],:]
	# k_segments = get_segments_for_track(track_k[1:poi_df.intersect[index],:])
	# start_energy = poi_df.distance[index]/last(track_full.distance) * max_energy

	# orig_track = track_full[1:poi_df.intersect[index],:]

	points_to_analyze = points_k[1:poi_df.k[index]]
	new_speed_vector_for_regular = []
	for i=2:size(points_to_analyze,1)
		append!(
			new_speed_vector_for_regular,
			fill(
				k_speeds[i-1],
				points_to_analyze[i] - points_to_analyze[i-1]
			)
		)
	end
	# println("new size $(size(new_speed_vector_for_regular,1))")
	# println("old size $(poi_df.intersect[index]-1)")
	# @assert size(new_speed_vector_for_regular,1)==poi_df.intersect-1
	return new_speed_vector_for_regular
end

# ╔═╡ 31ac9c67-77f5-422c-b57a-6cf0c1021f9b
function draw_comparison_plot_peaks_speeds_on_reg(index)

	reg_plot = simulate_run(
		track_optimization_comparison_df.regular_speeds[index],
		track_full[1:poi_df.intersect[index],:],
		get_segments_for_track(track_full[1:poi_df.intersect[index],:]),
		poi_df.distance[index]/last(track_full.distance) * max_energy,
		start_datetime
	)


	peaks_plot = simulate_run(
		peaks_speeds_to_regular(index),
		track_full[1:poi_df.intersect[index],:],
		get_segments_for_track(track_full[1:poi_df.intersect[index],:]),
		poi_df.distance[index]/last(track_full.distance) * max_energy,
		start_datetime
	)

	k_plot = simulate_run(
		k_speeds_to_regular(index),
		track_full[1:poi_df.intersect[index],:],
		get_segments_for_track(track_full[1:poi_df.intersect[index],:]),
		poi_df.distance[index]/last(track_full.distance) * max_energy,
		start_datetime
	)

	plot(
		reg_plot, peaks_plot, k_plot,
		layout=(3,1),
		size=(750,1200),
		legend=false
	)
	
end

# ╔═╡ caa3f733-8766-4674-8d1e-f3b67fdb84a8
draw_comparison_plot_peaks_speeds_on_reg(4)

# ╔═╡ b3437d4d-eb4b-418b-a139-1bc5bf4e45a7
md"Похоже что упрощение по k работает отлично (при 1.75)

Peaks работает хуже чем полное представление трассы

После 6-го индекса оригинальное решение начинает расходиться. (но результат ещё лучше peaks)

После индекса 8 и peaks начинает работать лучше

На индексе 10 peaks начинает расходиться

k работает всегда"

# ╔═╡ 5cf54f86-1377-4342-a7b8-700af58bbe0a
md"Итого, k работает всегда и вообще пусечка

Peaks в малом интевале значений между 8 и 10 индексом"

# ╔═╡ 4f73ca8a-97d7-47af-bda8-003496ce38a4
md"### Рисуем новую картинку"

# ╔═╡ 800979ad-68c4-48ea-8bac-7bb0e6cb10e5
plot!(
	track_optimization_comparison_df.points_regular[good_k_index],
	track_optimization_comparison_df.opt_time_regular[good_k_index] ./ track_optimization_comparison_df.opt_time_k[good_k_index] ,
	ylabel="Performance gain",
	xlabel="Number of segments",
	labels=["k=$(k)" "2"],
	# title="Peaks performance",
	marker=:circle,
	legend=:topleft,
	# legend=false,
	# size=(400,300),
	# marginbottom=5Plots.mm
	bottom_margin = 10Plots.px
	# , layout=(2,1)
)

# ╔═╡ c7c4ad0e-600b-47f3-ad26-9f76cecf2931
md"# 3 - Алгоритм разбиения на подзадачи" 

# ╔═╡ 1eecf941-1a74-45f8-a104-1db2d6c8e377
md"Скорее всего в ноутбуку" 

# ╔═╡ b1effadf-79ca-422a-91e3-08f72b0239ad

function process_subtask_novelty!(subtask::Subtask, scaling_coef_variables:: Real, segments::DataFrame, iteration_num::Integer)
	# 3. each subtask comes with its own chunks (variables)			
	# split each task on parts
	subtask.variables_boundaries = calculate_boundaries(
		subtask.subtask_boundaries.from,
		subtask.subtask_boundaries.to,
		scaling_coef_variables
	)

	is_track_divisible_further = false
	if (subtask.subtask_boundaries.size) >= scaling_coef_variables
		# at least one subtask can be divided, so, there should be next iteration
		is_track_divisible_further = true
	end

	# TODO: how to calculate amount of speeds?
	vars_amount = size(subtask.variables_boundaries, 1)

	# TODO: change here
	# write a function that repeats array values to get another array of bigger size
	# if size(subtask.problem.initial_speeds, 1) == 2
	prev_iter_speeds = fill_array(
		subtask.problem.initial_speeds, 
		vars_amount
	)
	# else
	# 	prev_iter_speeds = fill_array(
	# 	35., 
	# 	vars_amount
	# 	)
	# end

	subtask_segments = get_segments_interval_typed(
		segments,
		subtask.subtask_boundaries.from,
		subtask.subtask_boundaries.to
	)

	function f_iter(input_speeds :: Vector{<: Real})
		return solar_partial_trip_wrapper_iter_typed(
		# return solar_partial_trip_wrapper_iter_with_low_energy(
			input_speeds, subtask_segments, subtask.variables_boundaries,
			subtask.problem.start_energy, subtask.problem.finish_energy,
			subtask.problem.start_datetime
		)
	end

	function f_iter_low_energy(input_speeds :: Vector{<: Real})
		# return solar_partial_trip_wrapper_iter(
		return solar_partial_trip_wrapper_iter_with_low_energy_typed(
			input_speeds, subtask_segments, subtask.variables_boundaries,
			subtask.problem.start_energy, subtask.problem.finish_energy,
			subtask.problem.start_datetime
		)
	end

	if size(subtask.problem.initial_speeds, 1) == 2
		td = TwiceDifferentiable(f_iter_low_energy, prev_iter_speeds; autodiff = :forward)
		random_term = fill(0., vars_amount)
		optim_func = f_iter_low_energy
	else
		td = TwiceDifferentiable(f_iter, prev_iter_speeds; autodiff = :forward)
		# random_term = rand(vars_amount)
		random_term = fill(1., vars_amount)
		optim_func = f_iter
	end

	# lower_bound = fill(0.0, vars_amount)
	# upper_bound = fill(100.0, vars_amount)
	# tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
	# result = optimize(td, tdc, prev_iter_speeds 
	# .+ random_term .- 0.5
	# 	,
	# 	IPNewton(),
	# 	Optim.Options(
	# 		x_tol = 1e-12,
	# 		f_tol = 1e-12,
	# 		g_tol = 1e-12
	# 	)
	# )

	result = optimize(
		optim_func,
		prev_iter_speeds 
		# .+ random_term .- 0.5
		,
		NelderMead()
		# autodiff = :forward
	)

	minimized_speeds = Optim.minimizer(result)
	# println(Optim.minimizer(result))
	# println(Optim.minimum(result))
	subtask.solution = minimized_speeds
	return is_track_divisible_further
end

# ╔═╡ 1fb8c4d4-1cdb-4ee8-860c-dd82cc0e7043
function iterative_optimization_novelty(
		track :: DataFrame,
		segments :: DataFrame,
		scaling_coef_subtasks :: Integer,
		scaling_coef_subtask_input_speeds :: Integer,
		start_energy :: Real,
		start_datetime=DateTime(2022,7,1,0,0,0)::DateTime
	)
	# general algorithm:
	# 0. data setup
	# 1. exit loop check
	# 2. 	split the track into subtasks
	# 3. 		each subtask comes with its own chunks (variables)
	# 4. 		solve optimization problem for every subtask with its chunks (variables)
	# 5. 	tie everything together (collect speeds to one array)
	# 6. 	make full simulation with said speeds
	# 7. 	prepare data for next iteration. should be subtask's chunks as subtasks
	# 8. 	go to 2:
	# 9. final calculations?

	# 0. data setup
	track_size = size(track,1)
	scaling_coef_variables = scaling_coef_subtasks * scaling_coef_subtask_input_speeds

	start_speeds = minimize_single_speed(
		track,
		segments,
		start_energy,
		start_datetime,
		31.
	)

	start_n_speeds = minimize_n_speeds(
		track,
		segments,
		2,
		start_energy,
		start_datetime,
		first(start_speeds)
	)

	# boundaries = calculate_boundaries(1, size(track, 1), scaling_coef_subtask_input_speeds)



	# solar_trip_boundaries(
	# 		convert_kmh_to_ms(iteration_speeds),
	# 		segments,
	# 		start_datetime
	# 	)


	# for i in eachindex(start_n_speeds)
	# 	subtask = Subtask(
	# 		boundaries[i],
	# 		# calculate_boundaries(1, track_size, scaling_coef),
	# 		[],
	# 		SubtaskProblem(
	# 			start_energy,
	# 			0.,
	# 			start_speeds,
	# 			start_datetime
	# 		),
	# 		[]
	# 	)
	# end

	zero_subtask = Subtask(
		Boundaries(1, track_size),
		# calculate_boundaries(1, track_size, scaling_coef),
		[],
		SubtaskProblem(
			start_energy,
			0.,
			# start_speeds,
			start_n_speeds,
			start_datetime
		),
		[]
	)

	iteration_1 = Iteration(
		[ zero_subtask ],
		1,
		IterationSolution(
			[],
			[],
			[],
			[]
		)
	);

	iterations::Vector{Iteration} = []
	push!(iterations, iteration_1)

	iteration_num = 1;
	# 1. exit loop check
	is_track_divisible_further = true
	# while iteration_num <= 2
	while is_track_divisible_further # && iteration_num <= 2
	# while is_track_divisible_further && iteration_num <= 3

		iteration = iterations[iteration_num]
		println("Iteration $(iteration.number)")
		
		# 2. split the track into subtasks
		# or grab the result of previous division for subtasks
		# подумать где у меня точки, где участки, а где трассы по результатам деления
		# и ещё подумать как получаются подзадачи и их переменные на разных итерациях
		# № итерации, подзадачи, количество переменных
		# 1			1			scale
		# 2			scale		scale^2
		# ...
		# n			scale^(n-1)	scale^n - но не совсем, где как места хватит

		# по идее достаточно делить трассу, когда определяемся с участками для подзадач
		# потом в конце итерации после всех оптимизаций склеивать единый массив из точек разделения
		# и на следующей итерации использовать уже его


		is_track_divisible_further = false
		# # is_there_single_subtask_where_track_is_divisible = false
		# # Threads.@threads for subtask_index in eachindex(iteration.subtasks)
		# # Threads.@threads for subtask in iteration.subtasks
		# for subtask in iteration.subtasks
		# # @floop for subtask in iteration.subtasks
		# # for subtask_index in eachindex(iteration.subtasks)
		# # for subtask_index in eachindex(subtasks_splits_iteration)
		# 	# println("Subtask $subtask_index")
		# 	# subtask = iteration.subtasks[subtask_index]

		# 	# println("Analyzing subtask from $(subtask.subtask_boundaries.from) to $(subtask.subtask_boundaries.to)")

		# 	# 3. each subtask comes with its own chunks (variables)			
		# 	# split each task on parts
		# 	subtask.variables_boundaries = calculate_boundaries(
		# 		subtask.subtask_boundaries.from,
		# 		subtask.subtask_boundaries.to,
		# 		scaling_coef_variables
		# 	)

		# 	if (subtask.subtask_boundaries.size) >= scaling_coef_variables
		# 		# at least one subtask can be divided, so, there should be next iteration
		# 		is_track_divisible_further = true
		# 	end

		# 	#######################

		# 	# TODO: how to calculate amount of speeds?
		# 	vars_amount = size(subtask.variables_boundaries, 1)

		# 	# TODO: change here
		# 	# write a function that repeats array values to get another array of bigger size
		# 	prev_iter_speeds = fill_array(
		# 		subtask.problem.initial_speeds, 
		# 		vars_amount
		# 	)

		# 	# prev_iter_speeds = fill(
		# 	# 	first(subtask.problem.initial_speeds), 
		# 	# 	vars_amount
		# 	# )

		# 	subtask_segments = get_segments_interval(
		# 		segments,
		# 		subtask.subtask_boundaries.from,
		# 		subtask.subtask_boundaries.to
		# 	)

		# 	# 4. solve optimization problem for every subtask with its chunks (variables)
		# 	# function f_iter(input_speeds)
		# 	# 	if iteration_num == 1
		# 	# 		return solar_partial_trip_wrapper_iter_with_low_energy(
		# 	# 		# return solar_partial_trip_wrapper_iter_with_low_energy(
		# 	# 			input_speeds, subtask_segments, subtask.variables_boundaries,
		# 	# 			subtask.problem.start_energy, subtask.problem.finish_energy,
		# 	# 			subtask.problem.start_datetime
		# 	# 		)
		# 	# 	else
		# 	# 		return solar_partial_trip_wrapper_iter(
		# 	# 		# return solar_partial_trip_wrapper_iter_with_low_energy(
		# 	# 			input_speeds, subtask_segments, subtask.variables_boundaries,
		# 	# 			subtask.problem.start_energy, subtask.problem.finish_energy,
		# 	# 			subtask.problem.start_datetime
		# 	# 		)
		# 	# 	end
		# 	# end

		# 	function f_iter(input_speeds)
		# 		return solar_partial_trip_wrapper_iter(
		# 		# return solar_partial_trip_wrapper_iter_with_low_energy(
		# 			input_speeds, subtask_segments, subtask.variables_boundaries,
		# 			subtask.problem.start_energy, subtask.problem.finish_energy,
		# 			subtask.problem.start_datetime
		# 		)
		# 	end

		# 	function f_iter_low_energy(input_speeds)
		# 		# return solar_partial_trip_wrapper_iter(
		# 		return solar_partial_trip_wrapper_iter_with_low_energy(
		# 			input_speeds, subtask_segments, subtask.variables_boundaries,
		# 			subtask.problem.start_energy, subtask.problem.finish_energy,
		# 			subtask.problem.start_datetime
		# 		)
		# 	end

		# 	if iteration_num == 1
		# 		td = TwiceDifferentiable(f_iter_low_energy, prev_iter_speeds; autodiff = :forward)
		# 	else
		# 		td = TwiceDifferentiable(f_iter, prev_iter_speeds; autodiff = :forward)	
		# 	end

		# 	# td = TwiceDifferentiable(f_iter, prev_iter_speeds; autodiff = :forward)
		# 	lower_bound = fill(0.0, vars_amount)
		# 	upper_bound = fill(100.0, vars_amount)
		# 	# upper_bound = fill(200.0, vars_amount)
		# 	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
		# 	# line_search = LineSearches.BackTracking();
		# 	# result = optimize(td, fill(speed, vars_amount),
		# 		#Newton(; linesearch = line_search),
		# 	result = optimize(td, tdc, prev_iter_speeds 
		# 	.+ rand(vars_amount) .- 0.5
		# 		,
		# 		IPNewton(),
		# 		Optim.Options(
		# 			x_tol = 1e-12,
		# 			f_tol = 1e-12,
		# 			g_tol = 1e-12
		# 		)
		# 	)
		# 	minimized_speeds = Optim.minimizer(result)
		# 	# println(Optim.minimizer(result))
		# 	# println(Optim.minimum(result))
		# 	subtask.solution = minimized_speeds

		# 	# TODO: check optimization procedure in compliance with article
		# 	# TODO: save result somewhere - in subtask struct
		# 	# OR, in subtaskResult struct

		# end
		# is_track_divisible_further = !is_there_single_subtask_where_track_is_divisible
		# push!(subtasks_splits_general, variables_split_iteration)

		# @floop for subtask in iteration.subtasks
		start_time_next_subtask = start_datetime
		# @floop for subtask in iteration.subtasks
		@showprogress for subtask in iteration.subtasks
			subtask.problem.start_datetime = start_time_next_subtask
			is_divisible = process_subtask_novelty!(subtask, scaling_coef_variables, segments, iteration_num)
			# @reduce(is_track_divisible_further |= is_divisible)
			is_track_divisible_further |= is_divisible
			# calculate start time for next optimization
			subtask_segments = get_segments_interval_typed(
				segments,
				subtask.subtask_boundaries.from,
				subtask.subtask_boundaries.to
			)
			subtask_speeds = convert_kmh_to_ms(
				set_speeds_boundaries(subtask.solution, subtask.variables_boundaries)
			)
			subtask_times = subtask_segments.diff_distance ./ subtask_speeds
			subtask_time = sum(subtask_times)
			start_time_next_subtask = subtask.problem.start_datetime .+ Dates.Millisecond.(round.(subtask_time .* 1000))
			# start_datetime .+ Dates.Millisecond.(round.(time_seconds .* 1000))
		end

		println()
		# 5. tie everything together (collect speeds to one array)
		iteration_speeds = []
		for subtask in iteration.subtasks
			speed_vector = set_speeds_boundaries(subtask.solution, subtask.variables_boundaries)
			append!(iteration_speeds, speed_vector)
		end

		# 6. make full simulation with said speeds
		power_use, solar_power, time_seconds = solar_trip_boundaries(
			convert_kmh_to_ms(iteration_speeds),
			segments,
			start_datetime
		)

		println("solar sum $(sum(solar_power))")
		println("use sum $(sum(power_use))")
		

		energy_in_system = []
		push!(energy_in_system, start_energy)
		total_energy = start_energy .+ solar_power .- power_use
		append!(energy_in_system, total_energy)
		println("min energy $(minimum(energy_in_system))")

		times = start_datetime .+ Dates.Millisecond.(round.(time_seconds .* 1000))
		pushfirst!(times, start_datetime)
		# times = travel_time_to_datetime(time_seconds, start_datetime)
		iteration.solution = IterationSolution(
			iteration_speeds,
			energy_in_system,
			time_seconds,
			times
		)
		# 7. 	prepare data for next iteration. should be subtask's chunks as subtasks

		# TODO: re-split on new subtasks
		# we can't use variables boundaries anymore

		if is_track_divisible_further
			iteration_num += 1;
			next_iteration = Iteration(
				[],
				iteration_num,
				IterationSolution(
					[],
					[],
					[],
					[]
				)
			)
			for subtask in iteration.subtasks
				new_subtasks_boundaries = calculate_boundaries(
					subtask.subtask_boundaries.from,
					subtask.subtask_boundaries.to,
					scaling_coef_subtasks
				)
				solution_index_counter = 0
				for new_subtask_index in eachindex(new_subtasks_boundaries)
					new_subtask_boundaries = new_subtasks_boundaries[new_subtask_index]
					variables_boundaries = calculate_boundaries(
						new_subtask_boundaries.from,
						new_subtask_boundaries.to,
						scaling_coef_subtask_input_speeds
					)
					# println(subtask.solution)
					# TODO: make proper speed selecting function
					input_speeds_subtask = subtask.solution[ solution_index_counter + 1 : solution_index_counter + length(variables_boundaries)]
					solution_index_counter += length(variables_boundaries)
					# input_speeds_subtask = subtask.solution[(variables_index-1)*scaling_coef_subtask_input_speeds+1:variables_index*scaling_coef_subtask_input_speeds]
					# println(input_speeds_subtask)
					new_subtask = Subtask(
						new_subtask_boundaries,
						[],
						SubtaskProblem(
							energy_in_system[new_subtask_boundaries.from],
							energy_in_system[new_subtask_boundaries.to],
							input_speeds_subtask,
							times[new_subtask_boundaries.from]
						),
						[]
					)
					push!(next_iteration.subtasks, new_subtask)
				end
			end

			push!(iterations, next_iteration)
		end
	end # 8. 	go to 2:
	# 9. final calculations?
	println("Calc done")
	return iterations

end

# ╔═╡ 4e3c8ce5-cbd2-4d4c-8bb9-67f4f8745068
@time res = iterative_optimization_novelty(
    track_full, segments_full,
    5,
    5,
    5100.,
    start_datetime
);

# ╔═╡ 50e63a4a-119a-4c17-b469-166e17acd2e3
plot_resulted_plan(
	track_full,
                segments_full,
                res[1].solution.speeds,
                res[1].solution.energies,
                res[1].solution.times,
                "iteration $(res[1].number)"
)

# ╔═╡ c8c23f2b-7690-4de3-bfe5-b98f323afcdc
simulate_run_finish_time(
	res[1].solution.speeds,
	track_full,
	segments_full,
	max_energy,
	start_datetime
)

# ╔═╡ 3aab8070-9b76-4e88-ba27-0555714f045d
plot_resulted_plan(
	track_full,
                segments_full,
                res[2].solution.speeds,
                res[2].solution.energies,
                res[2].solution.times,
                "iteration $(res[2].number)"
)

# ╔═╡ 54539a77-1e06-43d5-bc54-897377feb390
plot_resulted_plan(
	track_full,
                segments_full,
                res[3].solution.speeds,
                res[3].solution.energies,
                res[3].solution.times,
                "iteration $(res[3].number)"
)

# ╔═╡ bdfadbc5-b240-4616-b52b-36e7d477fa81
simulate_run_finish_time(
	res[3].solution.speeds,
	track_full,
	segments_full,
	max_energy,
	start_datetime
)

# ╔═╡ 409c92f6-cd3b-41c5-a69e-a2b3e16a3960
# plot_resulted_plan(
# 	track_full,
#                 segments_full,
#                 res[6].solution.speeds,
#                 res[6].solution.energies,
#                 res[6].solution.times,
#                 "iteration $(res[6].number)"
# )

simulate_run_finish_time(
	res[6].solution.speeds,
	track_full,
	segments_full,
	max_energy,
	start_datetime
)

# ╔═╡ 3bd9ea26-32ae-4c73-8606-546ebcbe8e99
md"## Высокая трасса с алгоритмом" 

# ╔═╡ b004e777-ba5c-4e64-8c2f-1c026fb09f0b
begin
	track_high = copy(track_full)
	track_high.altitude = track_high.altitude .* 10;
	segments_high = get_segments_for_track(track_high);
	segments_high.weather_coeff = weather_full_coeff;
end

# ╔═╡ 75d28697-8f3f-4c00-80a0-998b62694ab1
res_high = iterative_optimization_novelty(
    track_high, segments_high,
    5,
    5,
    5100.,
    start_datetime
);

# ╔═╡ a4c9be7c-fe85-4f78-9cfa-a805e3b29478
simulate_run_finish_time(
	res_high[1].solution.speeds,
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ ea926683-39a6-4bb3-a74f-21322a8a7812
simulate_run_finish_time(
	res_high[2].solution.speeds,
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ 6b103fc3-d4be-4834-9c10-d20c1210728d
simulate_run_finish_time(
	res_high[3].solution.speeds,
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ c856f853-f4b1-4f84-ac14-dc16ee86fdcd
simulate_run_finish_time(
	res_high[4].solution.speeds,
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ dead6b98-85cb-4ce6-85a4-a155bb3c0dae
res_high

# ╔═╡ 1291c857-d356-438d-a83a-b5ba5a257c79
simulate_run_finish_time(
	res_high[6].solution.speeds,
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ 08c98d56-de5a-449a-b5ab-0f330bc7034a
md"## Высокая трасса одна скорость" 

# ╔═╡ 863ed1cd-0f9b-42cd-b179-32b552fe477d
function minimize_single_speed_novelty(
	track, segments, start_energy, start_datetime, init_speed
)
	boundaries = calculate_boundaries(1, size(track, 1), 1)

	function f_single_speed(input_speed)
		# return solar_partial_trip_wrapper_iter(
		return solar_partial_trip_wrapper_iter_with_low_energy(
			input_speed, segments, boundaries,
			start_energy, 0.,
			start_datetime
		)
	end
	speeds = [init_speed]
	println("Calculating best single speed")
	result = optimize(
		f_single_speed,
		speeds 
		# .+ random_term .- 0.5
		,
		NelderMead()
	)
	
	minimized_speeds = Optim.minimizer(result)
	println("Got $(minimized_speeds) km/h")
	return minimized_speeds
end

# ╔═╡ b3b88186-beee-4e2b-bf16-c5149f481760
res_single_speed = minimize_single_speed_novelty(
	track_high,
	segments_high,
	5100.,
	start_datetime,
	40.
)

# ╔═╡ 65b59e55-9ace-47a5-9398-997e68da4120
simulate_run_finish_time(
	fill(res_single_speed[1], size(segments_high, 1)),
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ 4f630f4a-88f9-4cb8-8ba7-98e53ea9ac24
simulate_run(
	fill(res_single_speed[1], size(segments_high, 1)),
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ 6f1eccb3-9333-4884-a144-01a21f439301
md"Время в пути $(415020.821 / 3600.) часов"

# ╔═╡ 6fad67c5-ff1b-4280-bfcc-32b0e9c46f73
md"## Высокая трасса несколько скоростей"

# ╔═╡ 59adbef1-7e7d-4f11-844e-106b79c54358
function minimize_n_speeds_novelty(
	track :: DataFrame,
	segments :: DataFrame,
	n_variables :: Int64,
	start_energy :: Float64,
	start_datetime :: DateTime,
	init_speed :: Float64
)

	boundaries = calculate_boundaries(1, size(track, 1), n_variables)

	function f_speeds(input_speeds :: Vector{<: Real}) :: Real
		# return solar_partial_trip_wrapper_iter(
		return solar_partial_trip_wrapper_iter_with_low_energy_typed(
			input_speeds, segments, boundaries,
			start_energy, 0.,
			start_datetime
		)
	end
	speeds = fill(init_speed, n_variables)
	println("Calculating best speeds")
	td_0 = TwiceDifferentiable(f_speeds, speeds; autodiff = :forward)
	lower_bound_0 = fill(10.0, n_variables)
	upper_bound_0 = fill(100.0, n_variables)
	tdc_0 = TwiceDifferentiableConstraints(lower_bound_0, upper_bound_0)
	# line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, 1),
		#Newton(; linesearch = line_search),
	# result = optimize(td_0, tdc_0, speeds 
	# # .+ rand(1) .- 0.5
	# 	,
	# 	IPNewton(),
	# 	Optim.Options(
	# 		x_tol = 1e-12,
	# 		f_tol = 1e-12,
	# 		g_tol = 1e-12
	# 	)
	# )

	result = optimize(
		f_speeds,
		speeds 
		# .+ random_term .- 0.5
		,
		NelderMead()
	)
	
	
	minimized_speeds = Optim.minimizer(result)
	println("Got $(minimized_speeds) km/h")
	return minimized_speeds
end

# ╔═╡ 439d1038-d1e6-4160-9c71-aea97df8370c
n_variables = 50

# ╔═╡ e0416978-9f80-472c-a710-55baec62da36
res_n_speeds = minimize_n_speeds_novelty(
	track_high,
	segments_high,
	n_variables,
	5100.,
	start_datetime,
	40.
)

# ╔═╡ a8656ba2-d8f1-48c6-92a5-39f664f96f40
boundaries_n_speeds = calculate_boundaries(1, size(track_high, 1), n_variables)

# ╔═╡ 444da343-5611-4350-8741-39e1fe52041c
speed_vector_n_speeds = set_speeds_boundaries(res_n_speeds, boundaries_n_speeds)

# ╔═╡ cffdd67e-9d46-4ea9-a6ed-12d21ae9d8c3
simulate_run_finish_time(
	speed_vector_n_speeds,
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ 69a5937b-6c33-4661-adf2-619381b4b109
simulate_run(
	speed_vector_n_speeds,
	track_high,
	segments_high,
	max_energy,
	start_datetime
)

# ╔═╡ b7492047-8dbf-4e68-b7eb-1e18c614d356
md"Время в пути $(351778.687 / 3600.) часов"

# ╔═╡ a937b4ed-d5d2-4c44-87fc-8a05e1c34e02
md"# 4 - Всё разом "

# ╔═╡ 3af88372-a833-4b6c-b140-c31c2e2e6490
md"## С маршрутом Peaks"

# ╔═╡ 2a834193-ca43-459b-822d-c85c881ce872
md"Конструируем высокую трассу с peaks"

# ╔═╡ 8911b6b2-ebe6-41fa-ac5f-947f1bbd2097
begin
	track_peaks_high, segments_peaks_high = keep_extremum_only_peaks_segments(track_high)
	segments_peaks_high.weather_coeff = weather_coeff;
end

# ╔═╡ a8ce7c9f-1953-4d80-9336-fda268ee2bc7
md"Запускаем задачу оптимизации"

# ╔═╡ d300b1de-d91b-4c9c-ab7d-83f3c4e5e893
res_peaks_high = iterative_optimization(
    track_peaks_high, segments_peaks_high,
    5,
    5,
    max_energy,
    start_datetime
);

# ╔═╡ 4000c7e9-6ae3-4a2c-a414-92216b036e3c
md"35.7 секунд по времени " 

# ╔═╡ 749a3887-0237-47cc-90b0-170266a07aaf
simulate_run_finish_time(
	res_peaks_high[1].solution.speeds,
	track_peaks_high,
	segments_peaks_high,
	max_energy,
	start_datetime
)

# ╔═╡ 7f5c1bbd-292a-4635-9e0b-ca8414d19315
simulate_run_finish_time(
	res_peaks_high[2].solution.speeds,
	track_peaks_high,
	segments_peaks_high,
	max_energy,
	start_datetime
)

# ╔═╡ e38159a0-788f-443b-8030-0c6c85217a66
simulate_run_finish_time(
	res_peaks_high[3].solution.speeds,
	track_peaks_high,
	segments_peaks_high,
	max_energy,
	start_datetime
)

# ╔═╡ 0d65ff1d-be79-435e-8781-133dc5584680
simulate_run_finish_time(
	res_peaks_high[4].solution.speeds,
	track_peaks_high,
	segments_peaks_high,
	max_energy,
	start_datetime
)

# ╔═╡ 176813e3-aced-4c07-9bd7-f32f1cebf3dd
simulate_run_finish_time(
	res_peaks_high[5].solution.speeds,
	track_peaks_high,
	segments_peaks_high,
	max_energy,
	start_datetime
)

# ╔═╡ 7f24f371-6bae-4649-ba8f-9c3a1b5e52ed
simulate_run_finish_time(
	res_peaks_high[6].solution.speeds,
	track_peaks_high,
	segments_peaks_high,
	max_energy,
	start_datetime
)

# ╔═╡ 16eecde5-e534-45e2-b6c6-206d936f0cef
md"Красота!" 

# ╔═╡ a61db582-c7ad-49e5-8d9c-9363cf62cdc1
simulate_run(
	res_peaks_high[6].solution.speeds,
	track_peaks_high,
	segments_peaks_high,
	max_energy,
	start_datetime
)

# ╔═╡ 719ac3a8-400d-490c-b576-147eed156782
last(res_peaks_high[6].solution.seconds)

# ╔═╡ a5592cb8-51dc-4648-a722-d2dbfb4743f3
md"Время в пути $(last(res_peaks_high[6].solution.seconds) / 3600.) часов"

# ╔═╡ 3fe01e09-b832-4fc7-a805-9a643f3cd352


# ╔═╡ 1c8f6b62-b72d-48db-ac72-2f16e61554a0
md"## С маршрутом k"

# ╔═╡ 6a023ad1-fe87-4db3-a27f-41d77bcf6d93
md"Конструируем высокую трассу с k"

# ╔═╡ e22b4c20-8bd5-43ae-a98b-9ec550ad54d2
begin
	track_k_high, points_k_high = parametrized_track_simplification(track_high, k);
	segments_k_high = get_segments_for_track(track_k_high);
	weather_k_coeff = calculate_weather_weights_for_segments(
	    w,
	    elat,
	    elon,
	    segments_k_high
	);
	segments_k_high.weather_coeff = weather_k_coeff;
end

# ╔═╡ ef829540-cd8d-4556-b1e0-46782f8e782f
md"Запускаем залдачу оптимизации"

# ╔═╡ 9357318b-fa4f-49d5-b08c-006c682d6a90
size(track_k_high,1)

# ╔═╡ 5562b93e-0ce8-48cf-aab4-9f81113d7186
res_k_high = iterative_optimization(
    track_k_high, segments_k_high,
    5,
    5,
    max_energy,
    start_datetime
);

# ╔═╡ f5f3a114-1edb-40a1-b571-9a38711f56cf
simulate_run_finish_time(
	res_k_high[1].solution.speeds,
	track_k_high,
	segments_k_high,
	max_energy,
	start_datetime
)

# ╔═╡ 2482f572-e095-4f2e-9fa0-273b0160acd6
simulate_run_finish_time(
	res_k_high[2].solution.speeds,
	track_k_high,
	segments_k_high,
	max_energy,
	start_datetime
)

# ╔═╡ 8849e213-e3c5-4f9f-9762-dd8519779018
simulate_run_finish_time(
	res_k_high[3].solution.speeds,
	track_k_high,
	segments_k_high,
	max_energy,
	start_datetime
)

# ╔═╡ ac4db70d-ee48-453b-82ec-98527914195d
simulate_run_finish_time(
	res_k_high[4].solution.speeds,
	track_k_high,
	segments_k_high,
	max_energy,
	start_datetime
)

# ╔═╡ 925ad924-a051-4000-a434-a7d3fb7b2381
simulate_run_finish_time(
	res_k_high[5].solution.speeds,
	track_k_high,
	segments_k_high,
	max_energy,
	start_datetime
)

# ╔═╡ ec754126-47e8-4b50-9680-0b127b7ca643
simulate_run_finish_time(
	res_k_high[6].solution.speeds,
	track_k_high,
	segments_k_high,
	max_energy,
	start_datetime
)

# ╔═╡ 84cfa1c8-f950-4fe5-9dd7-20c1cdcc40cd
md"Время в пути $(last(res_k_high[6].solution.seconds) / 3600.) часов"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
NLSolversBase = "d41bc354-129a-5804-8e4c-c37616107c6c"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
TimeZones = "f269a46b-ccf7-5d73-abea-4c690281aa53"
WebIO = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"

[compat]
CSV = "~0.10.11"
DataFrames = "~1.6.1"
Distributions = "~0.25.100"
LineSearches = "~7.2.0"
NLSolversBase = "~7.8.3"
Optim = "~1.7.7"
Peaks = "~0.4.4"
PlotlyJS = "~0.18.10"
Plots = "~1.39.0"
PlutoUI = "~0.7.52"
ProgressMeter = "~1.9.0"
StatsBase = "~0.34.0"
TimeZones = "~1.13.0"
WebIO = "~0.8.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "d2d85f5744e9caca19a66a6ce96f35abec79a75a"

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
git-tree-sha1 = "b1c61fd7e757c7e5ca6521ef41df8d929f41e3af"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.8"

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

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

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
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"

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

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "938fe2981db009f531b6332e31c58e9584a2f9bd"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.100"

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

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

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
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "a20eaa3ad64254c61eeb5f230d9306e937405434"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.6.1"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

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
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

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
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

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
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

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
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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
git-tree-sha1 = "0d097476b6c381ab7906460ef1ef1638fbce1d91"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.2"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

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
git-tree-sha1 = "4cc0c5a83933648b615c36c2b956d94fda70641e"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "821e918c170ead5298ff84bffee41dd28929a681"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.17"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "3129380a93388e5062e946974246fe3f2e7c73e2"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.18"

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

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "7452869933cd5af22f59557390674e8679ab2338"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.10"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

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
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "ee094908d720185ddbdc58dbe0c1cbe35453ec7a"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "7c29f0e8c575428bd84dc3c72ece5178caa67336"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.2+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "eeab25344bf9901146c0200a7ca64ea479f8bf5c"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.0"

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
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

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
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

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

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "d39314cdbaf5b90a047db33858626f8d1cc973e1"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.0.0+2023c"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a1f34829d5ac0ef499f6d84428bd6b4c71f02ead"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.0"

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
deps = ["Artifacts", "Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "89e64d61ef3cd9e80f7fc12b7d13db2d75a23c03"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.13.0"
weakdeps = ["RecipesBase"]

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

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
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

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

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

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
git-tree-sha1 = "0eef0765186f7452e52236fa42ca8c9b3c11c6e3"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.21"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "04a51d15436a572301b5abbb9d099713327e9fc4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.4+0"

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

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

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

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

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

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

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

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

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

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

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
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╠═cef231a9-6584-440f-926c-6d33d7df1250
# ╠═82888790-7eff-11ee-1ea8-69148277dca2
# ╠═621782b4-0289-41b3-8702-02b96ab5a37d
# ╠═651b3d97-232a-4d9a-8534-f35828db5380
# ╠═8e4c579d-d9cf-47be-8987-6e35f7aff9ff
# ╠═95d48aa1-d500-4ed0-85da-65184ccf3364
# ╠═309ee8a3-f136-46ab-bb05-9fae649ce940
# ╠═a2d44ee3-a967-4c3e-8c28-782c58384870
# ╠═3fb0db03-1cbe-4711-9ff7-2ccf8b0ee08a
# ╠═a48f9dc8-654b-4072-a8ea-27b0e812324e
# ╠═e5dfb032-07c2-4d2e-bdcb-370c6d9f7c82
# ╠═5a7b2e00-3819-4d48-8165-af5f0ab69011
# ╠═146ede7c-64dd-443e-922e-e5fddc788d55
# ╠═e39e5375-4d8f-4e24-b0d9-97d743428f0f
# ╠═bba824d1-3ee6-4ccb-8a0b-3ebf73a08099
# ╠═6e727497-bec2-4b77-93ff-2344512795aa
# ╠═c5ef2c9b-fc8f-450b-9edf-e543c85723c2
# ╠═8f0207f2-f7b9-4e88-acfa-bf142e50e6a9
# ╠═6a3c26ea-5de0-49e4-9067-af429173aed4
# ╠═94f81fb6-92ad-4b92-99cc-8453fe5b1589
# ╠═aef683c8-fd5b-479d-83b2-e047ba0f8713
# ╠═6a9ac47d-fee4-4609-b3ed-fc5598b88762
# ╠═35950f4c-a4f3-4e94-9ca9-b3e3c1550d09
# ╠═7a30a132-9b75-4aef-9752-c5fcefbca922
# ╠═6f54da25-f910-441e-8851-b4bf3e847a4e
# ╠═e2737024-4cb2-4eec-90a1-d1137bd5e487
# ╠═a82feccb-57c8-4ab8-a9de-ec6117f3a38c
# ╠═8608e4d1-8ff6-4c57-80eb-45e8b721c4b3
# ╠═9df91593-d0dc-42ed-a88a-b22c92f91486
# ╠═dfdf787e-5fd4-430f-875a-b538182b085e
# ╠═fbd63615-69b1-4f1b-bfae-3a295c46c406
# ╠═ea338a3a-620a-4971-80bf-cc47b576ca63
# ╠═50bdbde7-ed0f-437f-a563-d20eb4314f60
# ╠═068f2adf-63da-4cd9-9c17-aefc92e5f3b0
# ╠═84eaf793-abb4-4dd9-8708-2ccbf55806dc
# ╠═494a4eec-af31-47a3-b552-7faf371d81d2
# ╠═18cf4375-dff5-4912-a34b-d2b929475443
# ╠═08313d1e-7822-4139-9495-1641836716fe
# ╠═1b6df0b2-572b-4255-8dd5-ee42d1a9a08d
# ╠═c58758ea-4fd4-42ce-8a82-0d1a09f4f07b
# ╠═4ab6e5b4-d16c-4016-a4d8-d2b7c8db5b9b
# ╠═49bf8d7f-605d-4dbd-8077-8238573aea92
# ╠═4ccac8b1-070c-43cd-aeaf-69955d3ce2c3
# ╠═54ff073c-c751-4a13-b6be-05d744537291
# ╠═f405aaa5-e520-4d0c-bd63-ea4276a1be93
# ╠═978810c4-8f70-47c3-b1e3-0608ca855b1a
# ╠═2039e2ff-a56d-4a44-89ec-9073ae98480b
# ╠═f04149e4-a1f3-4a1d-bc3c-7c4e520ba976
# ╠═2106d6cb-9a70-4ccb-907d-3f157cfe1c46
# ╠═82deaf73-4205-4b43-8e7e-4f6500012834
# ╠═8bab99e6-b6de-4fa1-b5ae-0a31af1cf59a
# ╠═cedabc1e-0cb6-4920-9bed-200f85b861c7
# ╠═e0a006aa-cb6b-474d-91af-be792a03deaa
# ╠═a77ac935-1c8f-43af-b0ff-62bb39204562
# ╠═ebb6de9a-eedb-4352-af99-ba2f95c057d3
# ╠═cefe8ebf-cbec-40b1-a906-b4ca714538ae
# ╠═18b859fa-b004-4406-bbf3-81249e85be40
# ╠═aad06b97-aade-47ae-8196-58bb8a0f9c75
# ╠═70c353c1-c103-45ee-9ab0-b355326f96a6
# ╠═18d5923d-eb7f-48d1-80ea-5d147d99d8a5
# ╠═190d79b4-fcbc-4473-8e30-8316f2cded96
# ╠═3b499e5e-5235-4133-b35d-123f2ca18ada
# ╠═08493508-869b-46a2-aeb3-098b1ca99380
# ╠═674597af-873d-44f6-bab7-8960bfa19ee3
# ╠═4f298aa8-de16-4414-a794-e5a532c733b5
# ╠═a5688643-48a7-459d-8ee4-f32e3e86211c
# ╠═ab2f7e0b-8aa9-4ddd-aec5-16c3e5de70b3
# ╠═8295f55e-69ec-4a17-8cb6-5e5441a632ce
# ╠═c0b97d5e-4a54-41ad-83d3-0d8923290d70
# ╠═a477f042-55ac-4a8c-9c7d-fb3f6a85bc1a
# ╠═b0685b3a-42eb-4312-86a2-c9ab9432b565
# ╠═3d0e6774-90ad-4026-abf1-6feaa7579f30
# ╠═da104d8b-67f8-406f-b6b7-095a40596df8
# ╠═95d3b5a0-a4e6-4cc1-8b23-48e1d18961cf
# ╠═9d5fc300-afdc-4389-be3e-2f548fcafa76
# ╠═e2409a9e-580b-4fc4-bb1a-3d07e53dceb4
# ╠═530ec0eb-87d9-408b-9ca0-4f76ba29307d
# ╠═7c7f9b0b-3648-4ae7-8686-0546294337fe
# ╠═5b45416b-ab8e-4127-a637-99dbd76cef98
# ╠═439882df-7137-4157-89cf-d8ac92a6a8e6
# ╠═ad5d9933-b68d-4eb4-816f-b083306189ef
# ╠═bf28df12-d9bf-4310-951a-bd22eab4b12b
# ╠═d50301b9-de09-4b91-9104-af1a149ca6c0
# ╠═690226f5-a3df-4aad-8cc7-c1f16c0ab033
# ╠═c6798503-84ac-4d29-b1f2-f3a2037ced11
# ╠═89203abf-3941-41f1-8361-aae96f867c58
# ╠═3680d006-9322-41b9-93c6-51498d56889f
# ╠═fde27fab-aaa2-45bf-8c71-e92bb53dc25c
# ╠═6fcc5460-c724-4f6e-bb04-8785cddb7f79
# ╠═23906645-c779-4eb1-9201-6302fec045ea
# ╠═9ffbda77-5314-493e-89f6-bcaaff48a737
# ╠═62499f26-9863-47f4-9aab-37d2f6f1b255
# ╠═35b780d1-bf92-42bc-ac72-0a95e2585f56
# ╠═955df201-d9d4-4da2-b604-16619c92dd85
# ╠═127ca813-a652-4cc7-9995-588a0ab863cb
# ╠═b53eec1b-c644-4ef9-8086-6d90eb89935b
# ╠═f497ffca-5630-4d30-811f-bd1e15e8d34b
# ╠═1246d429-eb90-41df-bb14-df94f67a152f
# ╠═83934f3d-bf9c-4cd2-976c-0d50a91cbbb8
# ╠═d7ab2d22-c546-4b6c-b2d8-7aad0af65261
# ╠═e269d75f-bc01-4621-b193-555f54746ead
# ╠═93b53e5d-69e1-47d5-a6bc-4d1756467b2d
# ╠═990ac9a0-fc6f-4afb-aed2-e679ddcadfbf
# ╠═429fcd9c-300f-437b-ae15-16ff5040f3d3
# ╠═b11ef547-15b2-4425-bc83-fa6b8071e0a9
# ╠═bf08d4cd-d34b-4a42-82df-edbbb8d02437
# ╠═0c17cb82-91e8-423a-8764-87b768e71a46
# ╠═251b48a8-3a82-4859-a1b7-95306e484375
# ╠═4807e032-e508-4795-9ec8-e830d9abfcc5
# ╠═cf34f9ff-744a-498b-9835-b478a2c2d6c5
# ╠═9d1b383c-5dee-45b7-bfce-f1f4bd9ef1e1
# ╠═468f734c-5f23-420a-8b87-e068e44cc977
# ╠═db538803-f8c5-48f3-9051-adcfa1c4160d
# ╠═cc8085a3-8ce6-45b0-a353-2bc3eddc4381
# ╠═9482aafe-37b6-482d-ae11-a9a64b183bfd
# ╠═63989b8f-fa2f-4151-b2cf-394e79b80474
# ╠═362f6a0b-b2e2-4bf0-a74d-59c1ec38578a
# ╠═09032708-dbde-43b1-9d18-70a881d65b7d
# ╠═215e45a7-e1f3-45d5-af34-ec7e89fb96fb
# ╠═88694580-4e7f-4e35-8c33-e90bf5679c7e
# ╠═480e5187-19d7-486e-ad08-3a340bd7a760
# ╠═e3f80864-8fd2-48bd-acb8-8c91aba28e76
# ╠═985db9d8-b795-45c7-991f-8b68d81feee9
# ╠═b32881c8-59c1-4806-91e9-f6447f2aa22b
# ╠═d1146d90-664a-4235-95b0-8570e70713fc
# ╠═f8780df8-5868-49da-a5a8-a556afdb3f18
# ╠═75a0a6e3-4462-4f1d-baa6-3364d211a4b5
# ╠═a0775a71-0cca-43a6-835a-d66850605a45
# ╠═fd6a8bb0-7691-4a03-a928-999ddf6e3076
# ╠═f9d7c72f-c666-4a4a-b31a-c92b37507da9
# ╠═62043829-4ca6-4dd7-bf4f-85c7343ee177
# ╠═46c90566-c373-40c1-b4ca-1ae1c9872a87
# ╠═dcdf64c2-01df-4d73-a365-1b6381c9c3e1
# ╠═3a6eb336-4ec3-400a-b58b-418c7429ca2d
# ╠═2f0d3995-cbfd-420f-8936-380fc5d5da88
# ╠═a694736a-a239-4b3c-8670-9982aaf0204e
# ╠═714c393b-bb42-4622-848f-1dca1868169d
# ╠═a13bb6c3-0bed-4c2f-868f-5acded90aada
# ╠═75983ca9-e54a-419d-895c-892916135484
# ╠═75ec8028-db00-4e56-936c-dddf5c286735
# ╠═d105fc37-e079-4c5d-aa6b-288e4828c993
# ╠═2e35d429-cbb1-4e63-a03b-eaa9662229ab
# ╠═2d8a9498-aa45-46d8-947e-94da8a07c27f
# ╠═31ac9c67-77f5-422c-b57a-6cf0c1021f9b
# ╠═caa3f733-8766-4674-8d1e-f3b67fdb84a8
# ╠═b3437d4d-eb4b-418b-a139-1bc5bf4e45a7
# ╠═5cf54f86-1377-4342-a7b8-700af58bbe0a
# ╠═4f73ca8a-97d7-47af-bda8-003496ce38a4
# ╠═800979ad-68c4-48ea-8bac-7bb0e6cb10e5
# ╠═c7c4ad0e-600b-47f3-ad26-9f76cecf2931
# ╠═1eecf941-1a74-45f8-a104-1db2d6c8e377
# ╠═1fb8c4d4-1cdb-4ee8-860c-dd82cc0e7043
# ╠═b1effadf-79ca-422a-91e3-08f72b0239ad
# ╠═4e3c8ce5-cbd2-4d4c-8bb9-67f4f8745068
# ╠═50e63a4a-119a-4c17-b469-166e17acd2e3
# ╠═c8c23f2b-7690-4de3-bfe5-b98f323afcdc
# ╠═3aab8070-9b76-4e88-ba27-0555714f045d
# ╠═54539a77-1e06-43d5-bc54-897377feb390
# ╠═bdfadbc5-b240-4616-b52b-36e7d477fa81
# ╠═409c92f6-cd3b-41c5-a69e-a2b3e16a3960
# ╠═3bd9ea26-32ae-4c73-8606-546ebcbe8e99
# ╠═b004e777-ba5c-4e64-8c2f-1c026fb09f0b
# ╠═75d28697-8f3f-4c00-80a0-998b62694ab1
# ╠═a4c9be7c-fe85-4f78-9cfa-a805e3b29478
# ╠═ea926683-39a6-4bb3-a74f-21322a8a7812
# ╠═6b103fc3-d4be-4834-9c10-d20c1210728d
# ╠═c856f853-f4b1-4f84-ac14-dc16ee86fdcd
# ╠═dead6b98-85cb-4ce6-85a4-a155bb3c0dae
# ╠═1291c857-d356-438d-a83a-b5ba5a257c79
# ╠═08c98d56-de5a-449a-b5ab-0f330bc7034a
# ╠═863ed1cd-0f9b-42cd-b179-32b552fe477d
# ╠═b3b88186-beee-4e2b-bf16-c5149f481760
# ╠═65b59e55-9ace-47a5-9398-997e68da4120
# ╠═4f630f4a-88f9-4cb8-8ba7-98e53ea9ac24
# ╠═6f1eccb3-9333-4884-a144-01a21f439301
# ╠═6fad67c5-ff1b-4280-bfcc-32b0e9c46f73
# ╠═59adbef1-7e7d-4f11-844e-106b79c54358
# ╠═439d1038-d1e6-4160-9c71-aea97df8370c
# ╠═e0416978-9f80-472c-a710-55baec62da36
# ╠═a8656ba2-d8f1-48c6-92a5-39f664f96f40
# ╠═444da343-5611-4350-8741-39e1fe52041c
# ╠═cffdd67e-9d46-4ea9-a6ed-12d21ae9d8c3
# ╠═69a5937b-6c33-4661-adf2-619381b4b109
# ╠═b7492047-8dbf-4e68-b7eb-1e18c614d356
# ╠═a937b4ed-d5d2-4c44-87fc-8a05e1c34e02
# ╠═3af88372-a833-4b6c-b140-c31c2e2e6490
# ╠═2a834193-ca43-459b-822d-c85c881ce872
# ╠═8911b6b2-ebe6-41fa-ac5f-947f1bbd2097
# ╠═a8ce7c9f-1953-4d80-9336-fda268ee2bc7
# ╠═d300b1de-d91b-4c9c-ab7d-83f3c4e5e893
# ╠═4000c7e9-6ae3-4a2c-a414-92216b036e3c
# ╠═749a3887-0237-47cc-90b0-170266a07aaf
# ╠═7f5c1bbd-292a-4635-9e0b-ca8414d19315
# ╠═e38159a0-788f-443b-8030-0c6c85217a66
# ╠═0d65ff1d-be79-435e-8781-133dc5584680
# ╠═176813e3-aced-4c07-9bd7-f32f1cebf3dd
# ╠═7f24f371-6bae-4649-ba8f-9c3a1b5e52ed
# ╠═16eecde5-e534-45e2-b6c6-206d936f0cef
# ╠═a61db582-c7ad-49e5-8d9c-9363cf62cdc1
# ╠═719ac3a8-400d-490c-b576-147eed156782
# ╠═a5592cb8-51dc-4648-a722-d2dbfb4743f3
# ╠═3fe01e09-b832-4fc7-a805-9a643f3cd352
# ╠═1c8f6b62-b72d-48db-ac72-2f16e61554a0
# ╠═6a023ad1-fe87-4db3-a27f-41d77bcf6d93
# ╠═e22b4c20-8bd5-43ae-a98b-9ec550ad54d2
# ╠═ef829540-cd8d-4556-b1e0-46782f8e782f
# ╠═9357318b-fa4f-49d5-b08c-006c682d6a90
# ╠═5562b93e-0ce8-48cf-aab4-9f81113d7186
# ╠═f5f3a114-1edb-40a1-b571-9a38711f56cf
# ╠═2482f572-e095-4f2e-9fa0-273b0160acd6
# ╠═8849e213-e3c5-4f9f-9762-dd8519779018
# ╠═ac4db70d-ee48-453b-82ec-98527914195d
# ╠═925ad924-a051-4000-a434-a7d3fb7b2381
# ╠═ec754126-47e8-4b50-9680-0b127b7ca643
# ╠═84cfa1c8-f950-4fe5-9dd7-20c1cdcc40cd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
