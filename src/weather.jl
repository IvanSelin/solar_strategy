function generate_clouds(
	lat_from,
	lon_from,
	lat_to,
	lon_to,
	lat_peak,
	lon_peak,
	lat_std,
	lon_std,
	ndims,
	coef
)
	lat_distr = rand(Normal(lat_peak, 2), 10000 * ndims)
	lon_distr = rand(Normal(lon_peak, 2), 10000 * ndims)

	edges_lat = range(min(lat_from, lat_to), max(lat_from, lat_to), ndims + 1) #lat_from:step:lat_to
	edges_lon = range(min(lon_from, lon_to), max(lon_from, lon_to), ndims + 1)#lon_from:step:lon_to

	hist = fit(
		Histogram,
		(lat_distr, lon_distr),
		(edges_lat, edges_lon)
	)

	# return hist
	normed_weights = hist.weights / maximum(hist.weights) * coef#, edges_lat, edges_lon

	edges_lat_collected = collect(edges_lat)
	edges_lon_collected = collect(edges_lon)
	
	if lat_from > lat_to
		reverse!(normed_weights, dims=1)
		reverse!(edges_lat_collected)
	end

	if lon_from > lon_to
		reverse!(normed_weights, dims=2)
		reverse!(edges_lon_collected)
	end

	return normed_weights, edges_lat_collected, edges_lon_collected
	
	#hist.weights, hist.edges[1], hist.egdes[2]
end

function write_weather_json(weather, edges_lat, edges_lon, filename)
	weather_dict = Dict("weather"=>weather,"edges_lat"=>edges_lat,"edges_lon"=>edges_lon)
	open(joinpath("data",filename*".json"), "w") do f
		write(f, JSON.json(weather_dict))
	end
end

function read_weather_json(filename)
	# result_dict = Dict()
	# open(joinpath("data",filename*".json"), "r") do f
	# 	text_content = readall(f)  # file information to string
	# 	result_dict=JSON.parse(text_content)  # parse and transform data
	# end
	# return result_dict
	read_dict = JSON.parsefile(joinpath("data",filename*".json"))
	edges_lat = convert(Vector{Float64}, read_dict["edges_lat"])
	edges_lon = convert(Vector{Float64}, read_dict["edges_lon"])
	# weather_coeff = convert(Array{Float64,2}, read_dict["weather"])
	weather_coeff = reduce(hcat, convert(Vector{Vector{Float64}}, read_dict["weather"]))
	return weather_coeff, edges_lat, edges_lon
end


function calculate_weather_weights_for_segments(
    weather_weights, weather_edges_lat, weather_edges_lon, segments
)
    weather_coef = zeros(size(segments.latitude,1))
    for i in eachindex(weather_coef)
        _, lat_index = findmin(abs.(weather_edges_lat .- segments.latitude[i]))
        _, lon_index = findmin(abs.(weather_edges_lon .- segments.longitude[i]))
        lat_index = min(lat_index, size(w,1))
        lon_index = min(lon_index, size(w,1))
        weather_coef[i] = 1 .- weather_weights[lat_index, lon_index]
        # println("lat $(segments.latitude[i]), lon $(segments.longitude[i]), w is $(weather_coef[i])")
    end

    return weather_coef
end