function calculate_noon_time()
    days = [1:365;]
    b = 360 * (days .- 81) / 365
    t = 7.53 * cosd.(b) + 1.5 * sind.(b) - 9.87 * sind.(2b)
    return t
end

function solar_power_old(latitude, time_df)
    transmittance=0.75 # transmittance (unitless)
    solar_constant=1367 # solar constant (w/m^2)
    p=101325 # normal atmospheric pressure
    latitude = -12.438056

    sun_maximum_declination = 23.45 # degrees

    # # some obsolete lines
    # hours = collect(0:23) # or hours = [0:23;]
    # angle_by_hour = (12 .- hours) * 15 * pi / 180.0
    #
    # decline_angle = 23.45 * sin(2 * pi * (284 + day) / 365) * pi / 180

    # Global Horizontal Irradiance (GHI) = Direct Normal Irradiace (DNI, beam) +
    # + Diffuse Horizontal Irradiance (DHI) * cos(z)
    # OR
    # GHI = beam + diffuse + albedo (rad. refl. by the ground)
    # Plane of Array = beam rad * cos(z) + diffuse + reflected
    # should we use PoA if the panel array tilt angle is 0?

    # https://brstu.ru/static/unit/journal_smt/docs/number-36/115-121.pdf
    # 360/365 - rotation of sun per day in degrees
    # solar radiation intensity
    s_0 = solar_constant * (1 .+ 0.033 * cosd.(360 / 365.25 * time_df.year_day))
    plot(s_0, title = "Solar radiation intensity")
    # fixed formula? more sun during summer, day shift by 10
    # s_0 = solar_constant * (1 - 0.033 * cosd(360 / 365 * (day + 10))
    # or should there be cos instead of cosd?

    # sun declination  angle
    # why 285 and not 284 as it is stated in fourmula? days in year from zero?
    # https://susdesign.com/popups/sunangle/declination.php
    sun_declination_angle = (sun_maximum_declination *
        sind.(360 * ( (285 .+ time_df.year_day) / 365.25) )
        )
    plot(sun_declination_angle, title = "Sun declination angle")

    # equation of time (уравнение времени)
    # a sum of 2 sinusoids with periods of 1 year and 6 months
    # B = 360 * (time_df.year_day_float.-81)/365
    B = 360 * (time_df.year_day.-81)/365
    # equation of time itself, E(t) ~ E(day)
    E = 7.53 * cosd.(B) + 1.5 * sind.(B) - 9.87 * sind.(2B)
    plot(E, title = "Equation of time, minutes (day of the year)")

    # sun hour angle (ω)
    # sun_hour_angle = ( 15 * (hour - hour_zenith) + E(t) + (phi - phi_zone)
    sun_hour_angle = 15 * (time_df.day_seconds/60/60 .- 12) + E # E in minutes,
    # while hour angle is calculated for seconds?
    plot(sun_hour_angle, title = "Sun hour angle")

    # sunrise/sunset angle
    sunrise_angle = acosd.(- tand.(latitude) .* tand.(sun_declination_angle) )
    plot(sunrise_angle, title = "Sunrise angle for given latitude")

    # sun hours in a day
    sun_hours = 2 * sunrise_angle / 15
    plot(sun_hours, title = "Sun hours in a day for given latitude")

    # sun altitude angle (α)
    # need hour angle!
    # maybe hour angle is this?
    hours=[7,8,9,10,11,12,13,14,15,16,17]
    hour_angle=(12.0.-hours)*15.0
    # nah, too simple. calculate through equation of time?
    altitude_angle = ( ( cosd.(latitude) .* cosd.(sun_declination_angle) .*
        cosd.(hour_angle) ) .+ sind.(latitude) .* sind.(sun_declination_angle)
        )





    w_s = acosd(-tand(59.93)*tand(sun_decl_spb))
    sun_hours_test = (
        2 * w_s / 15
    )

    sin_alpha = ( sind.(latitude) .* sind.(sun_declination_angle) .+
        cosd.(latitude) .* cosd.(sun_declination_angle) .*
        cosd.()
        )

    # sun angle (6)
    # plot(sind.(360*(284 .+[1:365;])/365))

end


function solar_radiation_matlab()
    transmittance=0.75 # transmittance (unitless)
    solar_constant=1367 # solar constant (w/m^2)
    p=101325 # normal atmospheric pressure
    latitude = -12.438056
    day = 200

    hours=[7,8,9,10,11,12,13,14,15,16,17]
    hangle=(12.0.-hours)*15.0*pi/180.0 # radians
    plot(hangle, title="Hour angle")

    declangle=23.45*sin(2.0*pi*(284.0+day)/365.0)*pi/180.0 # in radians
    cosz=sind(latitude)*sin(declangle).+cosd(latitude)*cos(declangle)*cos.(hangle)
    plot(cosz, title="Cosz")

    m=p./(101.3*cosz) # optical airmass
    plot(m, title="Optical airmass")

    Sb=cosz*solar_constant.*(transmittance.^m)
    plot(Sb, title="Beam radiation")
    # Sb too small, e-150

    Sd=0.3*(1.0.-transmittance.^m)*solar_constant.*cosz
    plot(Sd, title="Diffuse radiation")

    St=Sb+Sd
    plot(St, title="Total radiation")

    x = zeros(24)
    x[7:17] = St/1000.0
    plot(x, title="Total radiation")

    return x

    # problems:
    # 1: amount of solar hours is fixed, which is wrong. it should be calculated
    # 2: beam radiation is too small
    # 3: no ground-reflected radiation
end


function solar_radiation_matlab_download()
    transmittance=0.75 # transmittance (unitless)
    solar_constant=1367 # solar constant (w/m^2)
    p=101325 # normal atmospheric pressure
    sun_maximum_declination = 23.45 # degrees
    latitude = -12.438056
    # day = 200
    days = 1:365

    s_0 = solar_constant * (1 .+ 0.033 * cosd.(360 / 365.25 * days))
    plot(s_0, title = "Solar radiation intensity")

    sun_declination_angle = (sun_maximum_declination *
        sind.(360 * ( (285 .+ days) / 365.25) )
        )
    plot(sun_declination_angle, title = "Sun declination angle degrees")
    plot(sun_declination_angle*pi/180.0, title = "Sun declination angle radians")


    cos_sunrise_angle = - tand.(latitude) .* tand.(sun_declination_angle)
    # fixes for polar circles. acos can't be calculated from values more than 1/-1
    cos_sunrise_angle[cos_sunrise_angle .> 1] .= 1
    cos_sunrise_angle[cos_sunrise_angle .< -1] .= -1
    # sunrise_angle = acosd.(- tand.(latitude) .* tand.(sun_declination_angle) )
    sunrise_angle = acosd.(cos_sunrise_angle)
    plot(sunrise_angle, title="Sunrise angle")
    plot(sunrise_angle*pi/180.0, title="Sunrise angle radians")

    days_length = 12*(1 .+ sunrise_angle / 180.0) .- 12*(1 .- sunrise_angle / 180.0)
    plot(days_length, title = "Days length (hours)")


    for day in days
        day = 200 # for debug purposes only
        # day_length = days_length[day]
        # sunrise_angle_day = sunrise_angle[day]
        # sun_declination_angle_day = sun_declination_angle[day]

        s_0 = solar_constant * (1 .+ 0.033 * cosd.(360 / 365.25 * day))
        sun_declination_angle = (sun_maximum_declination *
            sind(360 * ( (285 + day) / 365.25) )
            )
        cos_sunrise_angle = - tand(latitude) * tand(sun_declination_angle)
        # fixes for polar circles. acos can't be calculated from values more than 1/-1
        if cos_sunrise_angle > 1
            cos_sunrise_angle = 1
        end
        if cos_sunrise_angle < -1
            cos_sunrise_angle = -1
        end
        # sunrise_angle = acosd.(- tand.(latitude) .* tand.(sun_declination_angle) )
        sunrise_angle = acosd(cos_sunrise_angle)
        day_length = 12*(1 + sunrise_angle / 180.0) - 12*(1 - sunrise_angle / 180.0)

        step = 0.05
        # loop over sunshine hours
        for hour=0:step:day_length
            hour_angle = sunrise_angle*pi/180.0 - (pi * hour/day_length)
            # values are slightly off, check for +-1 errors
            # solar angle and azimuth
            sin_alpha = sind(latitude) * sind(sun_declination_angle) +
            cosd(latitude) * cosd(sun_declination_angle) * cos(hour_angle)
            # optical air mass - or air mass ratio
            M = sqrt(1229 + (614*sin_alpha)^2)-614*sin_alpha
            # atmoshperic transmitivity
            tau_b = 0.56 * (exp(-0.65 * M) + exp(-0.095 * M))
            # radiation diffusion coefficient for diffuse insolation
            tau_d = 0.271 - 0.294 * tau_b
            # reflectance transmitivity
            tau_r = 0.271 + 0.706 * tau_b

            cos_i = sind(sun_declination_angle) + cosd(sun_declination_angle) *
            cos(hour_angle) + cosd(sun_declination_angle) * sin(hour_angle)
        end
    end
end

function solar_radiation_pveducation()
    # starts here: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-outside-the-earths-atmosphere

    # solar constant
    H_constant = 1353 # W/m^2
    AM0 = 1366 # air mass zero, W/m^2
    # days of the year
    days = 1:365
    day = 200
    latitude = -12.438056 # Darwin
    longitude = 130.841111
    latitude = 59.9375 # Saint Petersburg
    longitude = 30.308611

    # radiant power density otside the Erath's athmosphere in W/m^2
    H = H_constant * (1 .+ 0.033*cosd.( (360 * (days .- 2)) / 365 ) )
    plot(days, H, title = "Radial power density W/m^2")

    # calculate air mass - how much air must pass through
    # phi not defined
    air_mass = 1 / (cos(phi) + 0.50572 * (96.07995 - theta)^-1.6364 )

    # direct intensity
    intesnity_direct = 1.353 * 0.7 ^ air_mass ^ 0.678 # kW/m^2
    # H_constant * 70% of radiation ^ air_mass ^ coefficient to fit the experimental data

    # with sea level data included
    indensity_direct_sea_level = 1.353 * ( (1 - 0.14 * height_km) * 0.7 ^ air_mass ^ 0.678 +
        0.14 * height_km) # kW/m^2

    # diffuse radiation
    intensity_diffuse = intensity_direct * 0.1

    # TIME
    # local solar time - when sun is in zenith
    # local time - time according to time zone
    # local standard time meridian - offset to utc in hours (from local time)
    # equation of time - offset for local solar time because of elliptical orbit
    # so: all times should be in UTC
    # then corrected  to local time, and then to local summer time


    local_time = DateTime(2022,04,10,13,58)
    UTC_time = DateTime(2022,04,10,10,58)
    # time shift for current latitude, hours
    local_standard_time_meridian = 15 * Dates.Hour(local_time - UTC_time).value
    # local_standard_time_meridian = 15 * longitude

    # equation of time
    B = 360 / 365 * (days .- 81) # in degrees
    equation_of_time = 9.87 * sind.(2B) - 7.53 * cosd.(B) - 1.5 * sind.(B) # minutes
    plot(
        equation_of_time,
        title = "Equation of time",
        label = "Equation of time",
        xlabel = "Day of the Year",
        ylabel = "Minutes"
    )

    # because local standard time meridian is not real sun time , minutes
    time_correction_factor = 4 * (longitude - local_standard_time_meridian) + equation_of_time[day]
    time_correction_factor_days = 4 * (longitude - local_standard_time_meridian) .+ equation_of_time

    local_solar_time = local_time + Dates.Second(round(time_correction_factor * 60.0))
    # LST = 4 longitude + 60 UTC_minutes +EoT
    # so, we only need to know UTC time and longitude
    local_solar_time = UTC_time +
        Dates.Second( round( (4*longitude + equation_of_time[day]) * 60) )

    # hour angle, 0 at solar noon
    hour_angle = 15 * (local_solar_time - 12)
    minutes_day = 0:24*60
    hour_angle_test = 15 * (minutes_day/60.0 .- 12)
    plot(hour_angle_test, title = "Hour angle for local summer time")
    # TODO: hor angle calculation for UTC time

    sun_declination_angle = -23.45 * cosd.(360 / 365 * (days .+ 10))
    plot(sun_declination_angle, title = "Sun declination angle")
    # also consider using 23.45 * sind.(360/365 * (days .+ 284))

    # elevation (altitude) angle
    elevation_angle = 90 + latitude - sun_declination_angle

    elevation = asin.(
        sind.(sun_declination_angle) .* sind(latitude) +
        cosd.(sun_declination_angle) .* cosd(latitude) .* cosd.(hour_angle_test)
    )
    # TODO: hourly data inputs

    sunrise_hour = 12 .- 1 / 15.0 .*
        acosd.(
            (-sind(latitude) * sind.(sun_declination_angle) ) ./
            (cosd(latitude) * cosd.(sun_declination_angle) )
            )
        .- time_correction_factor_days / 60
    plot(sunrise_hour, title = "sunrise hour")

    sunset_hour = 12 .+ 1 / 15.0 .*
        acosd.(
            (-sind(latitude) * sind.(sun_declination_angle) ) ./
            (cosd(latitude) * cosd.(sun_declination_angle) )
            )
        .- time_correction_factor_days / 60
    plot(sunset_hour, title = "sunset hour")

    # azimuth?
    # azimuth = acosd.(
    #     () / cosd()
    # )

    # TODO: data inputs as lat, long, DateTime, which is later turned to day_hous and day_minutes
    # simulate the whole year
end


function solar_radiation_pvedication_time(data_df, track_dataframe)
    # starts here: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-outside-the-earths-atmosphere

    # copy generate_year_time_dataframe(100000)
    # data_df = copy(time_dataframe)
    data_df.latitude = track_dataframe.latitude
    data_df.longitude = track_dataframe.longitude
    data_df.altitude = track_dataframe.altitude

    # for debug purposes
    # data_df = copy(time_df)
    # data_df.latitude = track.latitude
    # data_df.longitude = track.longitude
    # data_df.altitude = track.altitude

    # setting up needed values for calculcations
    # data_df.day = Dates.day.(data_df.datetime)

    # solar constant
    H_constant = 1353 # W/m^2
    # AM0 = 1366 # air mass zero, W/m^2

    # radiant power density otside the Erath's athmosphere in W/m^2
    # H = H_constant * (1 .+ 0.033*cosd.( (360 * (Dates.dayofyear.(data_df.utc_time) .- 2)) / 365 ) )
    # plot(Dates.dayofyear.(data_df.utc_time),H, title = "Radial power density W/m^2")

    # TODO: air mass calculation with phi? and theta?

    # time shift for current latitude, hours
    # local_standard_time_meridian = 15 * Dates.Hour(local_time - UTC_time).value
    # very rough estimation
    # proper solution like https://stackoverflow.com/questions/5584602/determine-timezone-from-latitude-longitude-without-using-web-services-like-geona
    data_df.lstm = 15 * ceil.(data_df.longitude * 24 / 360)

    # equation of time
    B = 360 / 365 * (Dates.dayofyear.(data_df.utc_time) .- 81) # in degrees
    equation_of_time = 9.87 * sind.(2B) - 7.53 * cosd.(B) - 1.5 * sind.(B) # minutes
    # plot(
    #     equation_of_time,
    #     title = "Equation of time",
    #     label = "Equation of time",
    #     xlabel = "Day of the Year",
    #     ylabel = "Minutes"
    # )

    # TODO: resume from here
    # because local standard time meridian is not real sun time , minutes
    # time_correction_factor = 4 * (data_df.longitude - data_df.lstm) + equation_of_time
    # plot(time_correction_factor, title="Time correction factor")

    # local_solar_time = local_time + Dates.Second(round(time_correction_factor * 60.0))
    # # LST = 4 longitude + 60 UTC_minutes +EoT
    # # so, we only need to know UTC time and longitude
    # local_solar_time = UTC_time +
    #     Dates.Second( round( (4*longitude + equation_of_time[day]) * 60) )

    data_df.local_solar_time = data_df.utc_time +
        Dates.Second.( round.(4 * data_df.longitude + equation_of_time) * 60)
    
    # TODO: should hour angle be only in integer values? 
    # minutes_from_start_of_the_day = Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time);
    hour_angle = 15 * ((Dates.hour.(data_df.utc_time) * 60 .+ Dates.minute.(data_df.utc_time)) ./60 .- 12)
    # TODO: hour angle calculation for UTC time
    # plot(data_df.utc_time, hour_angle, title="Hour angle vs UTC time")
    # TODO: hour angle calculation for local solar time ??? check if needed

    # TODO: continue to sun declination angle
    # TODO: continous sun declination angle calculation? (as for now quantified by days)
    sun_declination_angle = -23.45 * cosd.(360 / 365 * (Dates.dayofyear.(data_df.utc_time) .+ 10))
    # plot(sun_declination_angle, title = "Sun declination angle")

    # elevation (altitude) angle
    # elevation_angle = 90 .+ data_df.latitude .- sun_declination_angle
    # plot(data_df.utc_time, elevation_angle, title="Elevation angle (deg)")

    # elevation?
    elevation = asin.(
        sind.(sun_declination_angle) .* sind.(data_df.latitude) +
        cosd.(sun_declination_angle) .* cosd.(data_df.latitude) .* cosd.(hour_angle)
    )
    # plot(elevation, title="Elevation of the sun in rad")
    # when elevation angle is negative it means that sun is below horizon 
    # so, we should nullify these elements
    elevation_filtered = copy(elevation)
    elevation_filtered[elevation_filtered .<=0] .= 0
    # plot(data_df.utc_time, elevation_filtered, title="elevation angle when sun above horizon")

    # TODO: sunrise and sunset hours

    # sunrise_hour = 12 .- 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60
    # # plot(data_df.utc_time, sunrise_hour, title = "sunrise hour")

    # sunset_hour = 12 .+ 1 / 15.0 .*
    #     acosd.(
    #         (-sind.(data_df.latitude) .* sind.(sun_declination_angle) ) ./
    #         (cosd.(data_df.latitude) .* cosd.(sun_declination_angle) )
    #         )
    #     .- time_correction_factor / 60

    # plot(data_df.utc_time, sunset_hour, title = "sunset hour")

    # calculating zenith angle, needed for optical air mass
    zenith_angle = pi/2 .- elevation_filtered
    # plot(data_df.utc_time, zenith_angle, title="Zenith angle of the sun in rad")
    cos_zenith_angle = cos.(zenith_angle)

    # solar radiation consists of direct and diffuse radiation
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass 

    # calculate air mass - how much air light must pass through
    air_mass_full = 1 ./ (cos_zenith_angle .+ 0.50572 .* (96.07995 .- deg2rad.(zenith_angle)).^-1.6364 )

    # direct intensity, not actually used
    # intensity_direct = H_constant .* 0.7 .^ (air_mass_full .^ 0.678) # W/m^2
    # H_constant * 70% of radiation ^ air_mass_full ^ coefficient to fit the experimental data
    # plot(data_df.utc_time, intensity_direct, title="Direct intensity at Earth ground, W/m^2")

    # with altitude (from sea level)
    intensity_direct_altitude = H_constant .* ( (1 .- 0.14 .* data_df.altitude ./ 1000) .* 0.7 .^ (air_mass_full .^ 0.678) .+ 0.14 .* data_df.altitude ./ 1000) # W/m^2
    # plot(data_df.utc_time, intensity_direct_altitude, title="Direct intensity at altitude, W/m^2")
    # diffuse radiation
    intensity_diffuse = intensity_direct_altitude * 0.1
    # global irradiance
    intensity_global = intensity_direct_altitude + intensity_diffuse
    # plot(data_df.utc_time, intensity_global, title = "Global irradiance at altitude")
    # irradiance * cos(zenith_angle) is incident radiation ?

    
    s_incident = intensity_global
    # TODO: calculate radiation on a tilted surface
    # can be done through perpendecular (incident) or horizontal
    # s_horizontal = s_incident .* sin.(elevation)
    # module_angle_rad - at what angle to surface panels are located
    data_df.module_angle_rad .= 0.0 # just on the ground
    data_df.azimuth_angle .= 0.0 # just something, since it is on the ground, will not be used
    # TODO: calculate sun azimuth angle according to https://www.pveducation.org/pvcdrom/properties-of-sunlight/azimuth-angle 
    azimuth_cos = ( sind.(sun_declination_angle) .* cos.(data_df.latitude) .- 
    cosd.(sun_declination_angle) .* sin.(data_df.latitude) .* cos.(hour_angle) ) ./ cos.(elevation)
    azimuth_cos[azimuth_cos .> 1] .= 1
    azimuth_cos[azimuth_cos .< -1] .= -1
    data_df.azimuth_angle .= acos.( azimuth_cos )
    data_df.lst .= Dates.hour.(data_df.local_solar_time)
    data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ] .= 2*pi .- data_df.azimuth_angle[Dates.hour.(data_df.local_solar_time) .> 12 ]
    # module azimuth angle
    module_azimuth_angle = data_df.azimuth_angle
    # s_module stand for solar module, tilted surface
    # s_module = s_incident .* sin.(elevation .+ data_df.module_angle_rad)
    # plot(data_df.utc_time, s_module, title="Solar intensity without module azimuth")
    # OR https://www.pveducation.org/pvcdrom/properties-of-sunlight/arbitrary-orientation-and-tilt 
    s_module_azimuth = s_incident .* (
        cos.(elevation) .* sin.(data_df.module_angle_rad) .* cos.(data_df.azimuth_angle .- module_azimuth_angle) .+
        sin.(elevation) .* cos.(data_df.module_angle_rad)
        )
    # plot(data_df.utc_time, s_module_azimuth, title="Solar intensity with module azimuth")

    # next - simulate the race
    return s_module_azimuth
end

function solar_power_income(time_df, track_df, speed_vector)
    electrics_efficiency = 0.86
    solar_panels_efficiency = 0.228
    panels_area = 4 # m^2
    solar_intensity = solar_radiation_pvedication_time(time_df, track_df)
    power_income = electrics_efficiency .* solar_panels_efficiency .* solar_intensity .* panels_area .* track_df.diff_distance ./ speed_vector # W*s
    power_income_wt_h = power_income ./ 3600
    # power_income(i) = electrics_efficiency * panels_efficiency * solar_rad_h * panels_area * ...
    # ( dist_diff(i) * 3.6 / ( speed * 3600 ) ); % kWt*h
    return power_income_wt_h 
end

function calculate_power_income_accumulated(power_income)
    return cumsum!(power_income, power_income)
end