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
