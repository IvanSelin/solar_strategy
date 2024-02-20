struct SolarCar
    mass::Float64 # kg
    drag::Float64 # dimensionless
    frontal_area::Float64 # m²
    tire_friction::Float64 # dimensionless
    tire_friction_speed::Float64 # dimensionless
    engine_efficiency::Float64 # dimensionless
    electrics_power::Float64 # Wt
    electrics_efficiency::Float64 # dimensionless
    solar_panels_efficiency::Float64 # dimensionless
    panels_area::Float64 # m²
end

struct Environment
    g::Float64
    ro::Float64
end

# enabling broadcasting for custom types

Broadcast.broadcastable(sc::SolarCar) = Ref(sc)
Broadcast.broadcastable(e::Environment) = Ref(e)

function mechanical_power_calculation(speed_ms, slope, diff_distance)
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
        drag .* frontal_area .* speed_ms .^ 2 .* ro ./ 2 .+
        mass .* g .* (friction_1 .+ friction_2 * 4 .* speed_ms) .* cosd.(slope) .+
        mass .* g .* sind.(slope)
        )

    # mechanical power = mechanical force * distance delta / engine efficiency
    # watts * s
    mechanical_power = (
        mechanical_force .* diff_distance / (engine_efficiency)
        )

    # TODO: get rid of return, or at least make it type-stable
    # see https://docs.julialang.org/en/v1/manual/faq/#Types,-type-declarations,-and-constructors-1
    return mechanical_power
end

function mechanical_power_calculation_alloc(speed_ms, slope, diff_distance)
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
        drag * frontal_area * speed_ms ^ 2 * ro / 2. +
        mass * g * (friction_1 + friction_2 * 4 * speed_ms) * cosd(slope) +
        mass * g * sind(slope)
        )

    # mechanical power = mechanical force * distance delta / engine efficiency
    # watts * s
    mechanical_power = (
        mechanical_force * diff_distance / engine_efficiency
        )

    # TODO: get rid of return, or at least make it type-stable
    # see https://docs.julialang.org/en/v1/manual/faq/#Types,-type-declarations,-and-constructors-1
    return mechanical_power
end

function mechanical_power_calculation_alloc_typed(
        speed_ms :: Real, slope :: Real, diff_distance :: Real,
        solar_car :: SolarCar,
        env :: Environment
    ) :: Real
    # drag = 0.18
    # frontal_area = 1. # m^2
    # ro = 1.18 # air density

    # mass = 390. # kg
    # g = 9.8019 # at start: [41.2646201567207,-95.9244249307473,301.540649414063];
    # # 9.80147 at finish: [43.9660024736000,-121.345052439700,1229.07763671875]
    # friction_1 = 0.0023;
    # friction_2 = 0.000041; # total friction = friction_1 + friction_2*speed

    # engine_efficiency = 0.87

    # mechanical force = drag force + friction force + gravitational force
    # newtons
    mechanical_force = (
        solar_car.drag * solar_car.frontal_area * speed_ms ^ 2 * env.ro / 2. +
        solar_car.mass * env.g * (solar_car.tire_friction + solar_car.tire_friction_speed * 4 * speed_ms) * cosd(slope) +
        solar_car.mass * env.g * sind(slope)
        )

    # mechanical power = mechanical force * distance delta / engine efficiency
    # watts * s
    mechanical_power = (
        mechanical_force * diff_distance / solar_car.engine_efficiency
        )

    # TODO: get rid of return, or at least make it type-stable
    # see https://docs.julialang.org/en/v1/manual/faq/#Types,-type-declarations,-and-constructors-1
    return mechanical_power
end

function mechanical_power_calculation_typed(
    speed_ms :: Vector{<: Real},
    slope :: Vector{<: Real},
    diff_distance :: Vector{<: Real}
    )
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
        drag .* frontal_area .* speed_ms .^ 2 .* ro ./ 2 .+
        mass .* g .* (friction_1 .+ friction_2 * 4 .* speed_ms) .* cosd.(slope) .+
        mass .* g .* sind.(slope)
        )

    # mechanical power = mechanical force * distance delta / engine efficiency
    # watts * s
    mechanical_power = (
        mechanical_force .* diff_distance / (engine_efficiency)
        )

    # TODO: get rid of return, or at least make it type-stable
    # see https://docs.julialang.org/en/v1/manual/faq/#Types,-type-declarations,-and-constructors-1
    return mechanical_power
end

function electrical_power_calculation(diff_distance, speed_ms)
    power_onboard = 40; # Wt, 0.04kWt
    return power_onboard .* diff_distance ./ speed_ms;
end

function electrical_power_calculation_typed(
        diff_distance :: Vector{<: Real},
        speed_ms :: Vector{<: Real},
        solar_car :: SolarCar
    )
    return solar_car.electrics_power .* diff_distance ./ speed_ms;
end

function calculate_power_use_accumulated(mechanical, electrical)
    power_use = mechanical .+ electrical;
    power_use_accumulated = cumsum(power_use);
    return  power_use_accumulated / 3600;
end

function calculate_power_use(mechanical, electrical)
    power_use = mechanical .+ electrical;
    return  power_use / 3600;
end