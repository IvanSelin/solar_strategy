function travel_time_to_datetime(time_s, start_datetime)
    start_date = Dates.Date(start_datetime)
	start_seconds = Dates.hour(start_datetime) * 3600 + Dates.minute(start_datetime) * 60 + Dates.second(start_datetime)
    time_expanded = []
    push!(time_expanded, start_seconds)
	time_s .+= start_seconds
    append!(time_expanded, time_s)

	daily_start_hour_time = 8
    daily_finish_hour_time = 16

    start_time_seconds = daily_start_hour_time * 60 * 60
    finish_time_seconds = daily_finish_hour_time * 60 *60
    seconds_in_a_day = 24 * 60 * 60

    # adjust travel time so it happend only between daily start hour time and daily finish hour time
    # TODO: think of in-place operations to reduce memory consumption
    day_length = finish_time_seconds - start_time_seconds
    days_amount = last(time_expanded) ÷ day_length
	if first(time_expanded) < start_time_seconds
		δ = start_time_seconds - first(time_expanded)
		time_expanded .+= δ
	end
    for day=1:days_amount
		# @debug day*finish_time_seconds + (day-1)*seconds_in_a_day
		time_expanded[time_expanded .> finish_time_seconds + (day-1)*seconds_in_a_day] .+= (seconds_in_a_day - finish_time_seconds) + start_time_seconds
	end

    # day = time_s .÷ day_length .+ 1
    # time_s .+= start_time_seconds .* day .+
    # ( seconds_in_a_day .- finish_time_seconds) .* (day .- 1)

    # create a DataFrame for time information
    # adds seconds for 
    return Dates.DateTime(start_date) .+ Dates.Millisecond.(round.(time_expanded .* 1000))
end

function calculate_travel_time_seconds(speed_vector, track_df)
    #### time manipulation
    # time needed to spend on each part of the track
    time_s_intervals = track_df.diff_distance ./ speed_vector
    # base time, seconds driven from start of the race
    return cumsum(time_s_intervals)
end

function calculate_travel_time_seconds_intervals_typed(
    speed_vector :: Vector{<: Real},
    diff_distance :: Vector{<: Real})
    #### time manipulation
    # time needed to spend on each part of the track
    return diff_distance ./ speed_vector
end

function calculate_travel_time_seconds_intervals_typed_night(
    speed_vector :: Vector{<: Real},
    diff_distance :: Vector{<: Real},
    start_datetime :: DateTime)
    #### time manipulation
    # time needed to spend on each part of the track
    # TODO: night charge and time diff
    time_s = diff_distance ./ speed_vector
    return travel_time_to_datetime(time_s, start_datetime)
    # return diff_distance ./ speed_vector
end

function calculate_travel_time_datetime(speed_vector, segments, start_datetime)
    time_s = calculate_travel_time_seconds(speed_vector, segments)
	time_utc = travel_time_to_datetime(time_s, start_datetime)
	return DataFrame(utc_time=time_utc[2:end], time_s=time_s)
end

function calculate_travel_time_datetime_and_seconds_no_gaps(speed_vector, segments, start_datetime)
    time_s_segments = calculate_travel_time_seconds(speed_vector, segments)
    time_utc = start_datetime .+ Dates.Millisecond.(round.(time_s_segments .* 1000))
    pushfirst!(time_utc, start_datetime)
    return time_utc, time_s_segments
end