function solar_radiation_alloc_typed_vector(
    latitude :: Vector{<: Real},
    altitude :: Vector{<: Real},
    utc_time :: Vector{DateTime})
    # starts here: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-outside-the-earths-atmosphere

        # solar constant
    H_constant = 1353 # W/m^2
    # AM0 = 1366 # air mass zero, W/m^2

    # equation of time
    # B = 360 / 365 * (Dates.dayofyear(utc_time) - 81) # in degrees
    # equation_of_time = 9.87 * sind(2B) - 7.53 * cosd(B) - 1.5 * sind(B) # minutes

    # minutes_from_start_of_the_day = Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time);
    hour_angle = 15 .* ((Dates.hour.(utc_time) .* 60 .+ Dates.minute.(utc_time)) ./60 .- 12)

    sun_declination_angle = -23.45 .* cosd.(360 ./ 365 .* (Dates.dayofyear.(utc_time) .+ 10))

    elevation = asin.(
        sind.(sun_declination_angle) .* sind.(latitude) .+
        cosd.(sun_declination_angle) .* cosd.(latitude) .* cosd.(hour_angle)
    )
    # when elevation angle is negative it means that sun is below horizon 
    # so, we should nullify these elements
    # elevation_filtered = copy(elevation)
	# if elevation <= 0
	# 	elevation = 0
	# end

    for i in eachindex(elevation)
        if elevation[i] <= 0
            elevation[i] = 0
        end
    end

    # calculating zenith angle, needed for optical air mass
    zenith_angle = pi./2. .- elevation
    # plot(data_df.utc_time, zenith_angle, title="Zenith angle of the sun in rad")
    cos_zenith_angle = cos.(zenith_angle)

    # solar radiation consists of direct and diffuse radiation
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass 

    # calculate air mass - how much air light must pass through
    air_mass_full = 1 ./ (cos_zenith_angle .+ 0.50572 .* (96.07995 .- deg2rad.(zenith_angle)).^(-1.6364) )

    # with altitude (from sea level)
    intensity_direct_altitude = H_constant .* ( (1 .- 0.14 .* altitude ./ 1000) .* 0.7 .^ (air_mass_full .^ 0.678) .+ 0.14 .* altitude ./ 1000) # W/m^2
    # intensity_direct_altitude = H_constant .* ( (1 .- 0.14 .* data_df.altitude ./ 1000) .* 0.7 .^ (air_mass_full .^ 0.678) .+ 0.14 .* data_df.altitude ./ 1000) # W/m^2
    intensity_global = intensity_direct_altitude .* 1.1
    s_module_azimuth = intensity_global .* sin.(elevation)
    return s_module_azimuth
end

function solar_power_income_alloc_typed_vector(
    latitude :: Vector{<: Real},
    altitude :: Vector{<: Real},
    utc_time :: Vector{DateTime},
    time_s :: Vector{<: Real},
    solar_car :: SolarCar
    )
    # electrics_efficiency = 0.86
    # solar_panels_efficiency = 0.228
    # panels_area = 4 # m^2
	return solar_radiation_alloc_typed_vector(
        latitude, altitude, utc_time
        ) .* solar_car.electrics_efficiency .* solar_car.solar_panels_efficiency .* solar_car.panels_area .* time_s ./ 3600.0
end

function calculate_power_income_accumulated(power_income)
    return cumsum!(power_income, power_income)
end

function calculate_power_income_accumulated!(power_income)
    return cumsum!(power_income, power_income)
end