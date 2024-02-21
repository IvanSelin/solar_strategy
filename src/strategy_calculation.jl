using DataFrames

include("utils.jl")


mutable struct Boundaries
	from::Int64
	to::Int64
	size::Int64
	# points::DataFrame # maybe better another sub-struct or arrays?
	# segments::DataFrame
	# there is no sense in placing points(track) and segments here, 
	# since they will need another indexing
	# it is better to just use indexes

	# overloaded constructor will only work without @proto macro
	Boundaries(from, to) = new(from, to, to - from)
	Boundaries(from, to, size) = new(from, to, size)
end

mutable struct SubtaskProblem
	start_energy::Float64
	finish_energy::Float64
	initial_speeds::Vector{Float64}
	start_datetime::DateTime
end

mutable struct SubtaskVariable
	boundaries::Boundaries
	speed::Float64
end

struct IterationSolution
	speeds::Vector{Float64} # n-1 speeds (for segments)
	energies::Vector{Float64} # n energies (for points)
	seconds::Vector{Float64} # n-1 times (for segments)
	times::Vector{DateTime} # n times (for points)
end

mutable struct Subtask
	subtask_boundaries::Boundaries
	variables_boundaries::Vector{Boundaries}
	problem::SubtaskProblem
	solution::Vector{Float64} # n-1 speeds (for segments)
end

mutable struct Iteration
	subtasks::Vector{Subtask}
	number::Int64
	solution::IterationSolution
end

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

function calculate_split_indexes(start_index, end_index, chunks_amount)
	return start_index .-1 .+ calculate_split_indexes(end_index - start_index + 1, chunks_amount)
end

function calculate_boundaries(start_point, end_point, segments_amount)::Vector{Boundaries}
	# участки - segments, sections 
	length_in_segments = end_point - start_point;
	length_in_points = end_point - start_point + 1;

	# calculates split points (without start)
	split_points = calculate_split_indexes(length_in_segments, segments_amount) .+ start_point
	# so, we are adding start
	pushfirst!(split_points, start_point)

	boundaries_array::Vector{Boundaries} = []
	for i=2:length(split_points)
		segment = Boundaries(split_points[i-1], split_points[i], split_points[i] - split_points[i-1]);
		push!(boundaries_array, segment)
	end
	# return array of tuples?
	return boundaries_array
end

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
    power_use_accumulated_wt_h = mechanical_power + electrical_power
	cumsum!(power_use_accumulated_wt_h, power_use_accumulated_wt_h)
	power_use_accumulated_wt_h = power_use_accumulated_wt_h / 3600.

    # get solar energy income
	# @debug "track size is $(size(track.latitude, 1))"
    solar_power = solar_power_income_alloc.(
		track.latitude,
		track.longitude, 
		track.altitude, 
		time_df.utc_time,
		track.diff_distance,
		input_speed
	)
    solar_power_accumulated = calculate_power_income_accumulated(solar_power)
    energy_in_system = start_energy .+ solar_power_accumulated .- power_use_accumulated_wt_h

    # TODO: calculate night charging - do it later since it is not critical as of right now
    time_seconds = calculate_travel_time_seconds(input_speed, track)
    return power_use_accumulated_wt_h, solar_power_accumulated, energy_in_system, time_df, time_seconds
end

"Returns arrays of segments size (so, for task from 13 to 16th point, it will yield a size of 3)"
function solar_trip_boundaries_typed(
		input_speed :: Vector{<: Real},
		segments :: DataFrame,
		start_datetime :: DateTime,
		solar_car :: SolarCar,
		env :: Environment
	)
    # input speed in m/s

    #### calculcations
	mechanical_power = mechanical_power_calculation_alloc_typed.(
			input_speed, segments.slope, segments.diff_distance,
			solar_car, env
		)

    # electical losses
	electrical_power = electrical_power_calculation_typed(segments.diff_distance, input_speed, solar_car)
    # converting mechanical work to electrical power and then power use
    power_use_accumulated_wt_h = mechanical_power + electrical_power
	cumsum!(power_use_accumulated_wt_h, power_use_accumulated_wt_h)
	power_use_accumulated_wt_h = power_use_accumulated_wt_h ./ 3600.

	# calculating travel time
	intervals_seconds = segments.diff_distance ./ input_speed
	time_seconds = cumsum(intervals_seconds)
	mean_seconds = get_mean_data_typed(time_seconds)
	milliseconds = round.(mean_seconds .* 1000)
	mean_segment_utc = Vector{DateTime}(undef, size(input_speed, 1))
	mean_segment_utc .= start_datetime .+ Dates.Millisecond.(milliseconds)

    # get solar energy income
	solar_power = solar_power_income_alloc_typed_vector(
		segments.latitude,
		segments.altitude, 
		mean_segment_utc,
		intervals_seconds,
		solar_car
	)
	solar_power_adjusted = solar_power .* segments.weather_coeff
	calculate_power_income_accumulated!(solar_power_adjusted)

    # TODO: calculate night charging - do it later since it is not critical as of right now
	return power_use_accumulated_wt_h, solar_power_adjusted, time_seconds
end

function set_speeds(speeds, track, divide_at)
	println("set speed line 507")
	output_speeds = fill(last(speeds), size(track.distance, 1))
	for i=1:size(divide_at,1)-1
		if i==1
			output_speeds[1:divide_at[1]] .= speeds[1]
		else
			output_speeds[divide_at[i-1]:divide_at[i]] .= speeds[i]
		end
	end
	return output_speeds
end

function set_speeds_boundaries(speeds, variables_boundaries::Vector{Boundaries})
	if size(speeds,1) != size(variables_boundaries,1)
		throw(BoundsError("Input speeds and segments info mismatch!"))
	end
	output_speeds = []

	for i in eachindex(variables_boundaries)
		append!(
			output_speeds,
			fill(
				speeds[i],
				variables_boundaries[i].size
				)
			)
	end
	return output_speeds
end

function set_speeds_boundaries_typed(
	speeds :: Vector{T},
	variables_boundaries::Vector{Boundaries},
	output_length :: Integer
	) :: Vector{T} where {T <: Real} 
	# still having GC issues!
	if size(speeds,1) != size(variables_boundaries,1)
		throw(BoundsError("Input speeds and segments info mismatch!"))
	end
	output_speeds = Vector{T}(undef, output_length)
	counter_index::Int64 = 1
	counter_size::Int64 = 0

	for i in eachindex(variables_boundaries)
		# for j = 1:variables_boundaries[i].size
		# 	output_speeds[counter_index]=speeds[i]
		# 	counter_index += 1
		# end
		counter_size = variables_boundaries[i].size
		output_speeds[counter_index:(counter_index + counter_size - 1)] .= speeds[i]
		counter_index += counter_size
	end

	# for i in eachindex(variables_boundaries)
	# 	append!(
	# 		output_speeds,
	# 		fill(
	# 			speeds[i],
	# 			variables_boundaries[i].size
	# 			)
	# 		)
	# end
	return output_speeds
end

function solar_partial_trip_wrapper_alloc(speeds, track, indexes, start_energy, finish_energy, start_datetime)
	speeds_ms = convert_kmh_to_ms(speeds)
	speed_vector = set_speeds(speeds_ms, track, indexes)
	power_use, solar_power, energy_in_system, time, time_s = solar_trip_calculation_bounds_alloc(speed_vector, track, start_datetime, start_energy)
	cost = last(time_s) + 10 * abs(finish_energy - last(energy_in_system))
	return cost
	# return solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
end

function solar_partial_trip_wrapper_iter_typed(
		speeds :: Vector{<: Real},
		segments :: DataFrame,
		variables_boundaries :: Vector{Boundaries},
		start_energy :: Real,
		finish_energy :: Real,
		start_datetime :: DateTime,
		solar_car :: SolarCar,
		env :: Environment
	) :: Real
	speeds_ms :: Vector{<: Real} = convert_kmh_to_ms_typed(speeds)
	speed_vector :: Vector{<: Real} = set_speeds_boundaries_typed(speeds_ms, variables_boundaries, size(segments, 1))
	power_use, solar_power, time_s = solar_trip_boundaries_typed(
		speed_vector, segments, start_datetime,
		solar_car, env
	)
	# track points, not segments, that's why it is size is +1 
	# energy_in_system = start_energy .+ solar_power .- power_use
	last_energy_in_system = start_energy + last(solar_power) - last(power_use)

	# energy_capacity = 5100.

	# cost = sum(segments.diff_distance ./ speed_vector) + 100 * (finish_energy - last_energy_in_system)^2;
	cost = last(time_s) + 10000 * abs(finish_energy - last_energy_in_system);

	# cost = last(time_s) + (
	# 	10000 * (finish_energy - last(energy_in_system))^2 +
	# 	100 * max(0, maximum(energy_in_system) - energy_capacity)
	# )
	return cost
	# return solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
end

function solar_partial_trip_wrapper_iter_with_low_energy_typed(
		speeds :: Vector{<: Real},
		segments :: DataFrame,
		variables_boundaries :: Vector{Boundaries},
		start_energy :: Real,
		finish_energy :: Real,
		start_datetime :: DateTime,
		solar_car :: SolarCar,
		env :: Environment
	) :: Real
	speeds_ms :: Vector{<: Real} = convert_kmh_to_ms_typed(speeds)
	speed_vector :: Vector{<: Real} = set_speeds_boundaries_typed(speeds_ms, variables_boundaries, size(segments, 1))
	power_use, solar_power, time_s = solar_trip_boundaries_typed(
		speed_vector, segments, start_datetime, solar_car, env
	)
	# track points, not segments, that's why it is size is +1 
	# energy_in_system = start_energy .+ solar_power .- power_use
	last_energy = last(solar_power) - last(power_use) + start_energy
	solar_power -= power_use
	min_penalty = abs(minimum(solar_power) + start_energy)
	# not adding first energy element, since it is not needed for cost function
	# only last and minimum is needed
	# pushfirst!(energy_in_system, start_energy)

	# final_size = size(segments, 1) + 1
	# energy_in_system = Vector{}(undef, final_size)
	# energy_in_system[1] = start_energy
	# for i=2:final_size
	# 	energy_in_system[i] = start_energy + solar_power[i-1] - power_use[i-1]
	# end

	# energy_capacity = 5100.

	# cost = sum(segments.diff_distance ./ speed_vector) + 10000 * abs(minimum(energy_in_system)) + 100 * (finish_energy - last(energy_in_system));
	# cost = sum(segments.diff_distance ./ speed_vector) + 150000 * abs(minimum(energy_in_system))^2 + 10000 * (finish_energy - last(energy_in_system))^2;
	cost = sum(segments.diff_distance ./ speed_vector) + 150000 * min_penalty^2 + 10000 * (finish_energy - last_energy)^2;
	# println("f cost min energy $((minimum(energy_in_system)))")
	# cost = last(time_s) + (
	# 	10000 * (finish_energy - last(energy_in_system))^2 +
	# 	100 * max(0, maximum(energy_in_system) - energy_capacity)
	# )
	return cost
	# return solar_partial_trip_cost(speed_vector, track, start_energy, finish_energy, start_datetime)
end

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
	# line_search = LineSearches.BackTracking();
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
	split_energies = energy_in_system[split_indexes]
	pushfirst!(split_energies, start_energy)
	split_times = time[split_indexes, :utc_time]
	pushfirst!(split_times, start_datetime)
	
	# 5 - go though each track piece and enter function again
	# @debug "split_energies size is $(size(split_energies, 1)), chunks_amount is $(chunks_amount)"
	result_speeds = []
	for i=1:chunks_amount
		result_speeds_chunk = hierarchical_optimization_alloc(minimized_speeds[i], tracks[i], chunks_amount, split_energies[i], split_energies[i+1], split_times[i], iteration + 1 , start_index + split_indexes[i] - 1)
		append!(result_speeds, result_speeds_chunk)
	end

	return result_speeds
end

function fill_array(source_array, desired_size)
	source_size = size(source_array, 1)
	elements_left_to_repeat = mod(desired_size, source_size)
	# is_spreadable_evenly = elements_left_to_repeat == 0
	full_repeats = div(desired_size, source_size)
	output_array = []
	for elem in source_array
		append!(output_array, fill(elem, full_repeats))
		if elements_left_to_repeat > 0
			push!(output_array, elem)
			elements_left_to_repeat -= 1
		end
	end
	return Float64.(output_array)
end

function iterative_optimization(
		track :: DataFrame,
		segments :: DataFrame,
		scaling_coef_subtasks :: Integer,
		scaling_coef_subtask_input_speeds :: Integer,
		start_energy :: Real,
		solar_car::SolarCar,
		env::Environment,
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
		31.,
		solar_car,
		env
	)

	start_n_speeds = minimize_n_speeds(
		track,
		segments,
		2,
		start_energy,
		start_datetime,
		first(start_speeds),
		solar_car,
		env
	)

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

		# @floop for subtask in iteration.subtasks
		start_time_next_subtask = start_datetime
		# @floop for subtask in iteration.subtasks
		@showprogress for subtask in iteration.subtasks
			subtask.problem.start_datetime = start_time_next_subtask
			is_divisible = process_subtask!(subtask, scaling_coef_variables, segments, solar_car, env)
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
		power_use, solar_power, time_seconds = solar_trip_boundaries_typed(
			convert_kmh_to_ms(iteration_speeds),
			segments,
			start_datetime,
			solar_car,
			env
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

function process_subtask!(subtask::Subtask, scaling_coef_variables:: Real, segments::DataFrame,
		solar_car::SolarCar, env::Environment
	)
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
			subtask.problem.start_datetime,
			solar_car, env
		)
	end

	function f_iter_low_energy(input_speeds :: Vector{<: Real})
		# return solar_partial_trip_wrapper_iter(
		return solar_partial_trip_wrapper_iter_with_low_energy_typed(
			input_speeds, subtask_segments, subtask.variables_boundaries,
			subtask.problem.start_energy, subtask.problem.finish_energy,
			subtask.problem.start_datetime,
			solar_car, env
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

	lower_bound = fill(0.0, vars_amount)
	upper_bound = fill(100.0, vars_amount)
	tdc = TwiceDifferentiableConstraints(lower_bound, upper_bound)
	result = optimize(td, tdc, prev_iter_speeds 
	.+ random_term .- 0.5
		,
		IPNewton(),
		Optim.Options(
			x_tol = 1e-12,
			f_tol = 1e-12,
			g_tol = 1e-12
		)
	)

	# result = optimize(
	# 	optim_func,
	# 	prev_iter_speeds 
	# 	.+ random_term .- 0.5
	# 	,
	# 	# ConjugateGradient(),
	# 	# LBFGS(),
	# 	# autodiff = :forward
	# 	# NelderMead()
	# 	SimulatedAnnealing()
	# )

	minimized_speeds = Optim.minimizer(result)
	# println(Optim.minimizer(result))
	# println(Optim.minimum(result))
	subtask.solution = minimized_speeds
	return is_track_divisible_further
end

function minimize_single_speed(
		track::DataFrame,
		segments::DataFrame,
		start_energy::Real,
		start_datetime::DateTime,
		init_speed::Real,
		solar_car::SolarCar,
		env::Environment
	)

	boundaries = calculate_boundaries(1, size(track, 1), 1)

	function f_single_speed(input_speed :: Vector{<: Real})
		# return solar_partial_trip_wrapper_iter(
		return solar_partial_trip_wrapper_iter_with_low_energy_typed(
			input_speed,
			segments,
			boundaries,
			start_energy,
			0.,
			start_datetime,
			solar_car,
			env
		)
	end
	speeds = [init_speed]
	println("Calculating best single speed")
	td_0 = TwiceDifferentiable(f_single_speed, speeds; autodiff = :forward)
	lower_bound_0 = fill(5.0, 1)
	upper_bound_0 = fill(100.0, 1)
	# upper_bound_0 = fill(150.0, 1)
	tdc_0 = TwiceDifferentiableConstraints(lower_bound_0, upper_bound_0)
	# line_search = LineSearches.BackTracking();
	# result = optimize(td, fill(speed, 1),
		#Newton(; linesearch = line_search),
	result = optimize(td_0, tdc_0, speeds 
	# .+ rand(1) .- 0.5
		,
		IPNewton(),
		Optim.Options(
			x_tol = 1e-12,
			f_tol = 1e-12,
			g_tol = 1e-12
		)
	)
	minimized_speeds = Optim.minimizer(result)
	println("Got $(minimized_speeds) km/h")
	return minimized_speeds
end

function minimize_n_speeds(
	track :: DataFrame,
	segments :: DataFrame,
	n_variables :: Int64,
	start_energy :: Float64,
	start_datetime :: DateTime,
	init_speed :: Float64,
	solar_car::SolarCar,
	env::Environment
	)

	boundaries = calculate_boundaries(1, size(track, 1), n_variables)

	function f_speeds(input_speeds :: Vector{<: Real}) :: Real
		# return solar_partial_trip_wrapper_iter(
		return solar_partial_trip_wrapper_iter_with_low_energy_typed(
			input_speeds, segments, boundaries,
			start_energy, 0.,
			start_datetime,
			solar_car, env
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
	result = optimize(td_0, tdc_0, speeds 
	# .+ rand(1) .- 0.5
		,
		IPNewton(),
		Optim.Options(
			x_tol = 1e-12,
			f_tol = 1e-12,
			g_tol = 1e-12
		)
	)
	minimized_speeds = Optim.minimizer(result)
	println("Got $(minimized_speeds) km/h")
	return minimized_speeds
end

function calculate_boundaries_from_split_points(start_point, end_point, split_points)::Vector{Boundaries}
	splits_set = Set(split_points)
	push!(splits_set, start_point)
	push!(splits_set, end_point)

	splits_array = collect(splits_set)
	sort!(splits_array)

	boundaries_array::Vector{Boundaries} = []

	for i=2:length(splits_array)
		segment = Boundaries(splits_array[i-1], splits_array[i], splits_array[i] - splits_array[i-1]);
		push!(boundaries_array, segment)
	end
	# return array of tuples?
	return boundaries_array
end

function minimize_speeds_split_points_typed(
		track :: DataFrame,
		segments :: DataFrame,
		split_points :: Vector{Int64},
		start_energy :: Float64,
		start_datetime :: DateTime,
		init_speed :: Float64,
		solar_car :: SolarCar,
		env :: Environment
	) :: Vector{<: Real}
	boundaries = calculate_boundaries_from_split_points(1, size(track, 1), split_points)
	n_variables = size(boundaries, 1)
	
	function f_speeds(input_speeds :: Vector{<: Real}) :: Real
		# return solar_partial_trip_wrapper_iter(
		return solar_partial_trip_wrapper_iter_with_low_energy_typed(
			input_speeds, segments, boundaries,
			start_energy, 0.,
			start_datetime,
			solar_car, env
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
	result = optimize(td_0, tdc_0, speeds 
	# .+ rand(1) .- 0.5
		,
		IPNewton(),
		Optim.Options(
			x_tol = 1e-12,
			f_tol = 1e-12,
			g_tol = 1e-12
		)
	)
	minimized_speeds = Optim.minimizer(result)
	println("Got $(minimized_speeds) km/h")
	return minimized_speeds
end
