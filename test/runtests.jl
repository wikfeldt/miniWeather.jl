using miniWeather
using Test

@testset "test mass and energy conservation" begin
    nx = 100
    nz = 50
    sim_time = 400
    out_freq = 400
    weather_type = "thermal"

    model, grid = init(nx, nz, sim_time, weather_type);
    mass0, te0 = reductions(model, grid)
    run!(model, grid, out_freq, "output.nc")
    mass, te = reductions(model, grid)

    @test (mass - mass0) / mass0 < 1.0e-9
    @test (te - te0) / te0 < 4.5e-5

end
