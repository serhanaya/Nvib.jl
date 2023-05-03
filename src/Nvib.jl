__precompile__()

module Nvib

using ProgressMeter

export IPBeam, OPBeam, freqTab

include("analysis_type.jl")
include("frequency_table.jl")
include("coords.jl")
include("fundamental_matrix.jl")
include("detrmt.jl")

if !isfile(joinpath(homedir(), ".julia/environments/v1.7/Project.toml"))
    using Pkg
    Pkg.activate(joinpath(homedir(), ".julia/environments/v1.7"))
    Pkg.add("ProgressMeter")
end
using ProgressMeter

end # module
