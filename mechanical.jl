
function mechanical_work(speed_ms, slope, diff_distance)
    drag = 0.18
    frontal_area = 1 # m^2
    ro = 1.18 # air density

    mass = 390 # kg
    g = 9.8019 # at start: [41.2646201567207,-95.9244249307473,301.540649414063];
    # 9.80147 at finish: [43.9660024736000,-121.345052439700,1229.07763671875]
    friction_1 = 0.0023;
    friction_2 = 0.000041; # total friction = friction_1 + friction_2*speed

    engine_efficiency = 0.87

    # mechanical force = drag force + friction force + gravitational force
    # newtons
    mechanical_force = (
        drag * frontal_area * speed_ms^2 * ro / 2 .+
        .+ mass * g * (friction_1 + friction_2 * 4 * speed_ms) * cosd.(slope) .+
        .+ mass * g * sind.(slope)
        )

    # mechanical power = mechanical force * distance delta / engine efficiency
    # watts * s
    mechanical_power = (
        mechanical_force .* diff_distance / (engine_efficiency)
        )

    return mechanical_power
end
