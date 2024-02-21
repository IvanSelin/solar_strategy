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
    # points_df = select(track_csv, Not(:sin))
    points_df = select(track_csv, names(track_csv, x -> x!="sin"))
    

    return points_df, get_segments_for_track(points_df)
end

function keep_extremum_only_peaks_segments(track)
    # from Peaks
    # pks, vals = findmaxima(array)
    # pks, vals = findminima(array)
    # also https://juliapackages.com/p/findpeaks

    # track_copy = copy(track)
    max_altitude_peaks_indexes = argmaxima(track.altitude)
    min_altitude_peaks_indexes = argminima(track.altitude)
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
    max_altitude_peaks_indexes = argmaxima(track.altitude)
    min_altitude_peaks_indexes = argminima(track.altitude)
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
    if new_points[end] != size(track, 1)
        push!(new_points, size(track, 1))
    end
    return track[new_points,:], new_points
end

function get_average_on_segment(data_x, data_y, from, to)
	integrated_value::Float64 = 0
	# for i=1:length(data_y)-1
	for i=from:to-1
		integrated_value += (data_y[i] + data_y[i+1]) / 2 * (data_x[i+1] - data_x[i])
	end
	result = integrated_value / (data_x[to] - data_x[from])
	return result
end

function get_track_and_segments_for_selected_points_modified(track, points)
	# constructing proper lat, lon and alt values with numerical integration
	new_altitude = [];
	new_longitude = [];
	new_latitude = [];
	for i=1:length(points)-1
		push!(new_altitude, get_average_on_segment(track.distance, track.altitude, points[i], points[i+1]))
		push!(new_longitude,  get_average_on_segment(track.distance, track.longitude, points[i], points[i+1]))
		push!(new_latitude,  get_average_on_segment(track.distance, track.latitude, points[i], points[i+1]))
	end

	new_track = copy(track)
	new_track.index = 1:size(new_track, 1)
	new_track = new_track[points,:]
	
	new_segments = DataFrame(
        from = new_track.index[1:size(new_track, 1) - 1],
        to = new_track.index[2:size(new_track, 1)],
        diff_distance = diff(new_track.distance),
        diff_altitude = diff(new_track.altitude)
    )

	new_segments.slope = atand.(new_segments.diff_altitude ./ new_segments.diff_distance)

	# new_segments.latitude = get_mean_data(new_track.latitude)
    # new_segments.longitude = get_mean_data(new_track.longitude)
	new_segments.latitude = Float64.(new_latitude)
	new_segments.longitude = Float64.(new_longitude)
	# new_segments.altitude = get_mean_data(new_track.altitude)
	new_segments.altitude = Float64.(new_altitude)
	new_segments.weather_coeff .= 1.0
	
	return new_track, new_segments
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