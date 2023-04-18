using Peaks

function get_track_data(path_to_data)
    track_csv = CSV.read(path_to_data, DataFrame)
    
    ####### preprocessing of DataFrame
    # meters from km
    track_csv.distance = track_csv.distance * 1000
    # create new dataframe with size of n-1 (and without spoiled sin column)
    # sin column is ignored, since it is incomplete
    # btw, it's better to store the slope angle
    # track is for storing diff data, there are n-1 elements
    track = select(track_csv[2:size(track_csv,1),:], Not(:sin))
    track.diff_distance = diff(track_csv.distance)
    track.diff_altitude = diff(track_csv.altitude)
    track.slope = atand.(track.diff_altitude ./ track.diff_distance)
    return track
end

function get_track_and_segments(path_to_data)
    track_csv = CSV.read(path_to_data, DataFrame)
    
    ####### preprocessing of DataFrame
    # meters from km
    track_csv.distance = track_csv.distance * 1000

    # select everything except sin column
    points_df = select(track_csv, Not(:sin))

    return points_df, get_segments_for_track(points_df)
end

function keep_extremum_only(track)
    track_copy = copy(track)

    extremum_df = DataFrame()
    previous_row = first(track_copy)
    previous_diff_altitude = previous_row[:diff_altitude]
    # previous_index = 1
    # for i in 1:size(track, 1)
    for row in eachrow(track_copy)
        # if track[i,:diff_altitude] * previous_diff_altitude < 0
        if row.diff_altitude * previous_row.diff_altitude < 0
            # merge track pieces since last into one
            constructed_row = row
            # distance = track[previous_index:i, :distance] - track[previous_index, :distance]
            # constructed_row[:distance] = distance
            constructed_row.diff_distance = row.distance - previous_row.distance
            constructed_row.diff_altitude = row.altitude - previous_row.altitude
            constructed_row.slope = atand(constructed_row.diff_distance / constructed_row.diff_altitude)
            push!(extremum_df, constructed_row)
            previous_diff_altitude = previous_row.diff_altitude
            # previous_index = i
            previous_row = row
        end
    end
    # extremum_df.diff_distance = diff(extremum_df.distance)
    # extremum_df.diff_altitude = diff(extremum_df.altitude)
    # extremum_df.slope = atand.(extremum_df.diff_altitude ./ extremum_df.diff_distance)
    return extremum_df
end

function keep_extremum_only_peaks(track)
    # from Peaks
    # pks, vals = findmaxima(array)
    # pks, vals = findminima(array)
    # also https://juliapackages.com/p/findpeaks

    # track_copy = copy(track)
    max_altitude_peaks_indexes, max_altitude_peaks = findmaxima(track.altitude)
    min_altitude_peaks_indexes, min_altitude_peaks = findminima(track.altitude)
    peaks_indexes = []
    append!(peaks_indexes, max_altitude_peaks_indexes)
    append!(peaks_indexes, min_altitude_peaks_indexes)
    sort!(peaks_indexes)
    slice_df = track[peaks_indexes,:]


    temp_df = DataFrame(track[1,:])
    append!(temp_df, copy(slice_df))
    temp_df[1, :distance] = 0
    temp_df[1, :altitude] = 0

    slice_df.diff_distance = diff(temp_df.distance)
    slice_df.diff_altitude = diff(temp_df.altitude)
    slice_df.slope = atand.(slice_df.diff_altitude ./ slice_df.diff_distance)

    return slice_df
end

function keep_extremum_only_peaks_segments(track)
    # from Peaks
    # pks, vals = findmaxima(array)
    # pks, vals = findminima(array)
    # also https://juliapackages.com/p/findpeaks

    # track_copy = copy(track)
    max_altitude_peaks_indexes, max_altitude_peaks = findmaxima(track.altitude)
    min_altitude_peaks_indexes, min_altitude_peaks = findminima(track.altitude)
    peaks_indexes = []
    append!(peaks_indexes, max_altitude_peaks_indexes)
    append!(peaks_indexes, min_altitude_peaks_indexes)
    sort!(peaks_indexes)

    # add 1st and last point
    if first(peaks_indexes) != 1
        pushfirst!(peaks_indexes, 1)
    end
    if last(peaks_indexes) != size(track.distance, 1)
        push!(peaks_indexes, size(track.distance, 1))
    end
    track_peaks = track[peaks_indexes,:]

    # and make new segments
    # maybe in separate function, since this logic is used in get track data

    return track_peaks, get_segments_for_track(track_peaks)
end

function get_segments_for_track(track)
    segments_df = DataFrame(
        from = 1:size(track.distance,1) - 1,
        to = 2:size(track.distance,1),
        diff_distance = diff(track.distance),
        diff_altitude = diff(track.altitude)
    )

    segments_df.slope = atand.(segments_df.diff_altitude ./ segments_df.diff_distance)

    return segments_df
end

function get_track_interval(track, from_point, to_point)
    return track[from_point : to_point, :]
end

function get_segments_interval(segments, from_point, to_point)
    return segments[from_point : to_point - 1, :]
end
