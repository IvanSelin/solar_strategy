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