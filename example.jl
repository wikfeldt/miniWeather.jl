using miniWeather

nx_glob = 100
nz_glob = 50
sim_time = 1000.0
output_freq = 10.0
# possible options: collision, thermal, mountain_waves, turbulence, density_current, injection
#weather_type = "thermal"
weather_type = "density_current"


model, grid = init(nx_glob, nz_glob, sim_time, weather_type);

mass0, te0 = reductions(model, grid)

run!(model, grid, output_freq, "output1.nc")

mass, te = reductions(model, grid)
println("d_mass: ", (mass - mass0) / mass0)
println("d_te:   ", (te - te0) / te0)

#run!(model, grid, output_freq, "output2.nc")
#mass2, te2 = reductions(model, grid)
#println("d_mass: ", (mass2 - mass) / mass)
#println("d_te:   ", (te2 - te) / te)






