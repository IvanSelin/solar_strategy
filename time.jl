function travel_time_to_real_time(time_s)
    daily_start_time = 8
    daily_finish_time = 16

    daily_start_time_s = daily_start_time * 60 * 60
    daily_finish_time_s = daily_finish_time * 60 * 60

    daily_time_s = 24 * 60 * 60

    # initializing the array for days
    day = fill(0, length(time_s))

    # set the day length in seconds (depends on start and finish time)
    day_length_s = daily_finish_time_s - daily_start_time_s

    # create a DataFrame for time information
    time_df = DataFrame(seconds=time_s)

    # calculate day number for each time frame
    time_df.day = div.(time_df.seconds, day_length_s)

    # adjust running time to time-in-day
    time_df.day_seconds = time_df.seconds .+ daily_start_time_s .-
        .- day_length_s .* time_df.day

    return time_df
end
