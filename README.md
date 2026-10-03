# solar_strategy
Energy management for solar car racing written in Julia

Example usage can be seen in a src/strategy_example.jl file

## TODOs

List of stuff to do:
1. Think of a better way to switch between km/h and m/s. Either make a proper guidelines (but will require a lot of function variations for both m/s and km/h, maybe consider a macro/decorator?) or switch to proper units module with type checks and dynamic switching (introduce proper units typing with function accepting several unit and converting them on-the-fly)
2. 