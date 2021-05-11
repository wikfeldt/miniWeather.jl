# miniWeather

This package simulates weather-like flows and can be used for training in parallelizing 
and porting numerical Julia code to accelerated HPC architectures. It is directly derived 
from and inspired by [this miniWeather Mini App](https://github.com/mrnorman/miniWeather) 
developed by [Matt Norman](https://mrnorman.github.io), which provides example serial code 
in C, C++ and Fortran and solutions for how to parallelize it with MPI, OpenACC-Offloading, OpenMP-Threading and OpenMP-Offloading.

## Installation

The miniWeather package can be installed with the following command from the Julia REPL:
```julia
julia> using Pkg; Pkg.clone("https://github.com/wikfeldt/miniWeather.jl")
```

then it can be imported with the command:
```julia
julia> using miniWeather
```

## Usage

An example for how to use the package is found in this [example](./example.jl).

First load the package and set simulation parameters:
```julia
using miniWeather

# number of cells in x-direction
nx_glob = 100
# number of cells in x-direction
nz_glob = 50
# total simulation time in seconds
sim_time = 1000.0
# output frequency in seconds
output_freq = 10.0
# what type of phenomenon to simulate 
# (options: collision, thermal, mountain_waves, turbulence, density_current, injection)
weather_type = "thermal"
```

Then initialize the model and calculate an initial reduction (summation) of total 
mass and kinetic energy (which should be conserved):
```julia
model, grid = init(nx_glob, nz_glob, sim_time, weather_type);
mass0, te0 = reductions(model, grid)
```

Run the simulation with
```julia
# simulation results will be saved in netCDF file which can be visualized
run!(model, grid, output_freq, "output.nc")
```

Finally, compute another reduction and make sure that mass and kinetic energy 
are relatively well conserved (see below for rigorous integration testing):
```julia
mass, te = reductions(model, grid)
println("Δ_mass: ", (mass - mass0) / mass0)
println("Δ_te:   ", (te - te0) / te0)
```

## Testing

WRITEME

## Parallelization 

WRITEME

## License and credit

WRITEME