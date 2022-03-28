function calculate_noon_time()
    days = [1:365;]
    b = 360 * (days .- 81) / 365
    t = 7.53 * cosd.(b) + 1.5 * sind.(b) - 9.87 * sind.(2b)
    return t
end

function solar_power(latitude, time_df)
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
