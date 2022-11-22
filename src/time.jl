function travel_time_to_real_time(time_s)
    start_day = 140
    daily_start_time = 8
    daily_finish_time = 16

    daily_start_time_s = daily_start_time * 60 * 60
    daily_finish_time_s = daily_finish_time * 60 * 60

    day_seconds = 24 * 60 * 60

    # initializing the array for days
    day = fill(0, length(time_s))

    # set the day length in seconds (depends on start and finish time)
    day_length_s = daily_finish_time_s - daily_start_time_s

    # create a DataFrame for time information
    time_df = DataFrame(seconds=time_s)

    # calculate day number for each time frame
    time_df.day = Int64.(div.(time_df.seconds, day_length_s))
    time_df.year_day = time_df.day .+ start_day

    # adjust running time to time-in-day
    time_df.day_seconds = (time_df.seconds .+ daily_start_time_s .-
        day_length_s .* time_df.day
        )

    # day in float
    time_df.year_day_float = time_df.year_day .+ time_df.day_seconds ./ day_seconds

    # TODO: get rid of return
    return time_df
end

function travel_time_to_datetime(time_s)
    return travel_time_to_datetime(time_s, DateTime(2022,7,1,0,0,0))
    # TODO: make use of daily_start_hour_time and adjust millis for that
end

function travel_time_to_datetime(time_s, start_datetime)
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
    return start_datetime .+ Dates.Millisecond.(round.(time_s .* 1000))
end

function generate_year_time_dataframe(time_step_millis::Int64)

    start_datetime = DateTime(2022,1,1,0,0,0)
    end_datetime = DateTime(2022,12,31,23,59,59)
    # every 10 seconds (10000ms)
    interval_datetime = convert.(
        DateTime,
        Millisecond.(
            Dates.value(start_datetime):time_step_millis:Dates.value(end_datetime)
        )
    )
    time_df = DataFrame(utc_time=interval_datetime)
    return time_df
end

function calculate_travel_time_datetime(speed_vector, track_df)
    time_s = calculate_travel_time_seconds(speed_vector, track_df)
    # coverting the journey time to real time
    time_utc = travel_time_to_datetime(time_s)
    return DataFrame(utc_time=time_utc, time_s=time_s)
end

function calculate_travel_time_seconds(speed_vector, track_df)
    #### time manipulation
    # time needed to spend on each part of the track
    time_s_intervals = track_df.diff_distance ./ speed_vector
    # base time, seconds driven from start of the race
    return cumsum(time_s_intervals)
end

function calculate_travel_time_datetime(speed_vector, track, start_datetime)
    time_s = calculate_travel_time_seconds(speed_vector, track)
	time_utc = travel_time_to_datetime(time_s, start_datetime)
	return DataFrame(utc_time=time_utc, time_s=time_s)
end