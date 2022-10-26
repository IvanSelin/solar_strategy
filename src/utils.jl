function convert_single_kmh_speed_to_ms_vector(speed_kmh, len)
    return fill( speed_kmh / 3.6, len);
end

function convert_kmh_to_ms(speed_kmh)
    return speed_kmh / 3.6;
end

function propagate_speeds(speed_ms, track)
    if length(speed_ms) == size(track, 1)
        return speed_ms
    end
    # speed_ms size is a lot less than track size
    track_len = length(track.distance)
    speed_len = length(speed_ms)
    speeds = fill(last(speed_ms), track_len)
    sector_size = div(track_len, speed_len)
    for i = 1:speed_len-1
        speeds[(i-1)*sector_size + 1:i*sector_size] .= speed_ms[i]
    end
    return speeds
end

