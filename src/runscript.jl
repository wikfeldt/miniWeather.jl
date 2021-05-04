using Pkg
Pkg.activate("miniWeather")

using miniWeather

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

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
        "--data_spec_int"
        help = "How to initialize the data"
        default = Int(DATA_SPEC_THERMAL)
    end

    return parse_args(s)
end




function main()

    if !isinteractive()
        config = parse_commandline()
    else
        config = Dict(
            "nx_glob" => 100,
            "nz_glob" => 50,
            "sim_time" => 50,
            "output_freq" => 10,
            "data_spec_int" => Int(DATA_SPEC_THERMAL),
        )
    end

    println("Parsed args:")
    for (arg, val) in config
        println("  $arg  =>  $val")
    end

    model, grid = init(config["nx_glob"], config["nz_glob"], config["sim_time"], config["data_spec_int"])

    mass0, te0 = reductions(model, grid)
    println("mass = $mass0, te = $te0")

    #output(model, etime)

    #=
    anim_dens = Animation()
    anim_uwnd = Animation()
    anim_wwnd = Animation()
    anim_theta = Animation()
    anim = [anim_dens, anim_uwnd, anim_wwnd, anim_theta]
    #output_gif!(model, etime, grid, anim);
    =#
    
    run(model, grid, config["output_grid"], verbose=true)

    mass, te = reductions(model, grid)
    println("d_mass: ", (mass - mass0) / mass0)
    println("d_te:   ", (te - te0) / te0)
    
end



main()
