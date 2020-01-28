

function solar_power(latitude, time)
    transmittance=0.75 # transmittance (unitless)
    solar_constant=1367 # solar constant (w/m^2)
    p=101325 # normal atmospheric pressure

    sun_maximum_inclination = 23.45

    # TODO: calculate day from time
    # time as an input parameter, not day



    hours = collect(0:23) # or hours = [0:23;]
    angle_by_hour = (12 .- hours) * 15 * pi / 180.0

    decline_angle = 23.45 * sin(2 * pi * (284 + day) / 365) * pi / 180

    # Global Horizontal Irradiance (GHI) = Direct Normal Irradiace (DNI, beam) +
    # + Diffuse Horizontal Irradiance (DHI) * cos(z)
    # OR
    # GHI = beam + diffuse + albedo (rad. refl. by the ground)
    # Plane of Array = beam rad * cos(z) + diffuse + reflected
    # should we use PoA if the panel array tilt angle is 0?

    # https://brstu.ru/static/unit/journal_smt/docs/number-36/115-121.pdf
    # 360/365 - rotation of sun per day in degrees
    s_0 = solar_constant * (1 + 0.033 * cosd(360 / 365 * day))
    # fixed formula? more sun during summer, day shift by 10
    # s_0 = solar_constant * (1 - 0.033 * cosd(360 / 365 * (day + 10))
    # or should there be cos instead of cosd?

    # sun inclination angle
    sun_inclination_angle = sun_maximum_inclination *
        * sind.(360 * ( (284 .+ day) / 365) )



    sin_a = sind.(latitude) * sind.()

    # sun angle (6)
    # plot(sind.(360*(284 .+[1:365;])/365))

end
