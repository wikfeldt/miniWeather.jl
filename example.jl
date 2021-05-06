using miniWeather

nx_glob = 100
nz_glob = 50
sim_time = 500.0
output_freq = 10.0
data_spec_int = Int(DATA_SPEC_THERMAL)

model, grid = init(nx_glob, nz_glob, sim_time, data_spec_int);

mass0, te0 = reductions(model, grid)

simulate(model, grid, output_freq)

mass, te = reductions(model, grid)

println("d_mass: ", (mass - mass0) / mass0)
println("d_te:   ", (te - te0) / te0)



#miniWeather.run(model, grid)

