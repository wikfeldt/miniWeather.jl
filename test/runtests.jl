using miniWeather
using Test

@testset "test mass and energy conservation" begin
    nx = 100
    nz = 50
    sim_time = 400
    out_freq = 400
    data_spec_int = Int(DATA_SPEC_THERMAL)

    model, grid = init(nx_glob, nz_glob, sim_time, data_spec_int);
    mass0, te0 = reductions(model, grid)
    simulate(model, grid, output_freq)
    mass, te = reductions(model, grid)

    @test (mass - mass0) / mass0 < 1.0e-9
    @test (te - te0) / te0 < 4.5e-5

end
