module miniWeather

using ArgParse
using NetCDF
using ProgressMeter
import Base.Threads.@spawn
#using LoopVectorization: @avx

include("const.jl")
include("utils.jl")
include("Initialize.jl")
include("Timestep.jl")

export Grid, Model, init, run!, reductions

end
