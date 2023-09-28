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

function read_weather_json(path_to_file)
	# result_dict = Dict()
	# open(joinpath("data",filename*".json"), "r") do f
	# 	text_content = readall(f)  # file information to string
	# 	result_dict=JSON.parse(text_content)  # parse and transform data
	# end
	# return result_dict
	read_dict = JSON.parsefile(path_to_file)
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
        lat_index = min(lat_index, size(weather_weights,1))
        lon_index = min(lon_index, size(weather_weights,1))
        weather_coef[i] = 1 .- weather_weights[lat_index, lon_index]
        # println("lat $(segments.latitude[i]), lon $(segments.longitude[i]), w is $(weather_coef[i])")
    end

    return weather_coef
end

function generate_density_mapbox(w, edges_lat, edges_lon)
	ndims = size(w,1)
	w_arr = collect(Iterators.flatten(w))
	edges_lat_rep = repeat(get_mean_data(edges_lat), outer=ndims)
	edges_lon_rep = repeat(get_mean_data(edges_lon), inner=ndims)

	df = DataFrame(lat=edges_lat_rep, lon=edges_lon_rep, z=w_arr)
	traces_vector::AbstractVector{plotjs.AbstractTrace} = [];

	push!(
		traces_vector,
		plotjs.densitymapbox(
			lat=df.lat,
			lon=df.lon,
			z=df.z,
			opacity=0.5
		)
	)
	
	push!(
		traces_vector,
		plotjs.scattermapbox(
			lat=track_aus.latitude,
			lon=track_aus.longitude,
			marker_color="red",
			marker_size=1,
			mode="lines"
		)
	)
	
	plotjs.Plot(
		traces_vector,
		plotjs.Layout(
			width=650,
			height=600,
			geo_fitbounds="locations",
			autosize=true,
			# mapbox_style="stamen-terrain"
			mapbox_style=mapbox_style,
			mapbox_center_lat=-25.0,
			mapbox_center_lon=132.0,
			mapbox_zoom=3
		)
	)
end

function generate_density_data(w, edges_lat, edges_lon)
	ndims = size(w,1)
	w_arr = collect(Iterators.flatten(w))
	edges_lat_rep = repeat(get_mean_data(edges_lat), outer=ndims)
	edges_lon_rep = repeat(get_mean_data(edges_lon), inner=ndims)

	df = DataFrame(lat=edges_lat_rep, lon=edges_lon_rep, z=w_arr)
	return df
end

function generate_heatmap_traces(w, edges_lat, edges_lon)
	traces_vector::AbstractVector{plotjs.AbstractTrace} = [];
	push!(
		traces_vector,
		plotjs.scattermapbox(
			lat=track_aus.latitude,
			lon=track_aus.longitude,
			marker_color="red",
			marker_size=1,
			mode="lines"
		)
	)
	for i=1:length(edges_lat)-1
		for j=1:length(edges_lon)-1
			# println("lat: $(edges_lat[i]), lon: $(edges_lon[j]), w: $(w[i,j])")
			# println("lat+1: $(edges_lat[i+1]), lon+1: $(edges_lon[j+1]), w+1: $(w[i,j])")
			trace = plotjs.scattermapbox(
				fill="toself",
				lat = [
					edges_lat[i],
					edges_lat[i],
					edges_lat[i+1],
					edges_lat[i+1],
					edges_lat[i]
				],
				lon = [
					edges_lon[j],
					edges_lon[j+1],
					edges_lon[j+1],
					edges_lon[j],
					edges_lon[j]
				],
				# lat=[-11.,-11.,-12.,-12.5,-12., -12.],
				# lon=[130.,131.,131.,130.5,130., 130.],
				marker_size=1,
				# marker_color="orange",
				opacity=0.5,
				showlegend=false,
				# marker_colorscale=w[i,j]
				# marker_colorscale="Viridis",
				# marker_color=w[i,j]
				marker_color="rgb($(w[i,j]*255),$((1-w[i,j])*255),0)",
				name="$(w[i,j])"
			);
			push!(traces_vector, trace)
		end
		# println()
	end
	plotjs.Plot(
		traces_vector,
		plotjs.Layout(
			width=700,
			height=600,
			geo_fitbounds="locations",
			mapbox_style=mapbox_style,
			autosize=true,
			mapbox_center_lat=-25.0,
			mapbox_center_lon=132.0,
			mapbox_zoom=3
		)
	)
end