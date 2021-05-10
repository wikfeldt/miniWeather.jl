using miniWeather

using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    choices_weather_types = ["collision", "thermal", "mountain_waves", 
        "turbulence", "density_current", "injection"]
    @add_arg_table s begin
        "--nx_glob"
        help = "Total number of cells in the x-direction"
        arg_type = Int
        default = 100
        "--nz_glob"
        help = "Total number of cells in the z-direction"
        arg_type = Int
        default = 50
        "--sim_time"
        help = "How many seconds to run the simulation"
        arg_type = Int
        default = 50
        "--output_freq"
        help = "How frequently to output data to file (in seconds)"
        arg_type = Int
        default = 10
        "--weather_type"
        help = "Choose weather scenario and how to initialize the data. " * 
                "Must be one of " * join(choices_weather_types, ", ", " or ")
        range_tester = (x->x ∈ choices_weather_types)
        default = "thermal"
    end

    return parse_args(s)
end




function main()

    config = parse_commandline()

    println("Parsed args:")
    for (arg, val) in config
        println("  $arg  =>  $val")
    end

    model, grid = init(config["nx_glob"], config["nz_glob"], config["sim_time"], config["weather_type"])

    mass0, te0 = reductions(model, grid)
    println("mass = $mass0, te = $te0")
    
    run!(model, grid, config["output_freq"], verbose=true)

    mass, te = reductions(model, grid)
    println("d_mass: ", (mass - mass0) / mass0)
    println("d_te:   ", (te - te0) / te0)
    
end



main()
