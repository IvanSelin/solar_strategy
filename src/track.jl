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

function keep_extremum_only_peaks_segments_with_points(track)
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

    return track_peaks, get_segments_for_track(track_peaks), peaks_indexes
end

function parametrized_track_simplification(track, threshold)
    segments = get_segments_for_track(track);
    points = segments.from;
    points_to_delete = []
    # prev_row = segments[1,:];
    for i=1:size(segments,1) - 1
        if abs(segments.slope[i+1] - segments.slope[i]) < threshold
            push!(points_to_delete, segments.to[i]);
        end
    end
    new_points = setdiff(points, points_to_delete);
    return track[new_points,:], new_points
end

function get_mean_data(series)
    shifted = circshift(series,1)
    mean = (series .+ shifted) ./ 2
    return mean[2:end]
end

function get_mean_data_typed(series :: Vector{T}) where {T <: Real} 
    new_series_size = size(series,1)
    new_series = Vector{T}(undef, new_series_size)
    new_series[1] = 0.
    for i=2:new_series_size
        new_series[i] = (series[i-1] + series[i] ) / 2.
    end
    return new_series
end

function get_segments_for_track(track)
    segments_df = DataFrame(
        from = 1:size(track.distance,1) - 1,
        to = 2:size(track.distance,1),
        diff_distance = diff(track.distance),
        diff_altitude = diff(track.altitude)
    )

    segments_df.slope = atand.(segments_df.diff_altitude ./ segments_df.diff_distance)

    segments_df.latitude = get_mean_data(track.latitude)
    segments_df.longitude = get_mean_data(track.longitude)
    segments_df.altitude = get_mean_data(track.altitude)
    segments_df.weather_coeff .= 1.0 

    return segments_df
end

function get_peak_points(track)
    # находит max и min позиции, но не ловит плато.
    # то есть находит только первый элемент плоского плато, но не последний

    # переписал, стало чуть лучше, но всё равно находит не все плато
    max_altitude_peaks_indexes, max_altitude_peaks = findmaxima(track.altitude)
    min_altitude_peaks_indexes, min_altitude_peaks = findminima(track.altitude)
    reverse_alt = track.altitude[end:-1:1];
    reverse_max = argmaxima(reverse_alt);
    reverse_min = argminima(reverse_alt);
    reverse_max_reversed = length(track.altitude) .- reverse_max .+ 1;
    reverse_min_reversed = length(track.altitude) .- reverse_min .+ 1;
    peaks_indexes = Set();
    union!(peaks_indexes, Set(max_altitude_peaks_indexes))
    # append!(peaks_indexes, max_altitude_peaks_indexes)
    union!(peaks_indexes, Set(min_altitude_peaks_indexes))
    # append!(peaks_indexes, min_altitude_peaks_indexes)
    union!(peaks_indexes, Set(reverse_max_reversed))
    # append!(peaks_indexes, reverse_max_reversed)
    union!(peaks_indexes, Set(reverse_min_reversed))
    # append!(peaks_indexes, reverse_min_reversed)
    push!(peaks_indexes, 1)
    push!(peaks_indexes, length(track.altitude))
    peaks_array = sort(collect(peaks_indexes));
    # sort!(peaks_indexes)

    # # add 1st and last point
    # if first(peaks_indexes) != 1
    #     pushfirst!(peaks_indexes, 1)
    # end
    # if last(peaks_indexes) != size(track.distance, 1)
    #     push!(peaks_indexes, size(track.distance, 1))
    # end

    return peaks_array
end

function get_peak_points_plateau(altitudes)
    # не ловит пограничные случаи и плато
    points = Set();
    push!(points, 1);
    push!(points, length(altitudes));
    # current_peak = altitudes[1];
    # peak_positions = [];
    # max_positions = [];
    # push!(peak_positions, 1);
    for i = 2 : length(altitudes) - 1
        if altitudes[i] == altitudes[i-1] && altitudes[i] == altitudes[i+1]
            continue
        end
        if altitudes[i] >= altitudes[i-1] && altitudes[i] >= altitudes[i+1]
            # push!(max_positions, i)
            push!(points, i)
        end
        if altitudes[i] <= altitudes[i-1] && altitudes[i] <= altitudes[i+1]
            # push!(max_positions, i)
            push!(points, i)
        end
    end
    return sort(collect(points))
end


function get_track_and_segments_for_selected_points(track, points)
    temp_track = copy(track)
    temp_track.index = 1:size(temp_track, 1)
    temp_track = temp_track[points,:]



    segments_df = DataFrame(
        from = temp_track.index[1:size(temp_track, 1) - 1],
        to = temp_track.index[2:size(temp_track, 1)],
        diff_distance = diff(temp_track.distance),
        diff_altitude = diff(temp_track.altitude)
    )

    segments_df.slope = atand.(segments_df.diff_altitude ./ segments_df.diff_distance)

    segments_df.latitude = get_mean_data(temp_track.latitude)
    segments_df.longitude = get_mean_data(temp_track.longitude)
    segments_df.altitude = get_mean_data(temp_track.altitude)
    segments_df.weather_coeff .= 1.0 

    return temp_track, segments_df
end

function get_track_interval(track, from_point, to_point)
    return @view track[from_point : to_point, :]
end

function get_segments_interval(segments, from_point, to_point)
    return @view segments[from_point : to_point - 1, :]
end

function get_segments_interval_typed(segments::DataFrame, from_point::Integer, to_point::Integer)::DataFrame
    return segments[from_point : to_point - 1, :]
end

function get_average_on_segment(data_x, data_y, from, to)
	integrated_value = 0
	for i=from:to-1
		integrated_value += (data_y[i] + data_y[i+1]) / 2 * (data_x[i+1] - data_x[i])
	end
	result = integrated_value / (data_x[to] - data_x[from])
	return result
end