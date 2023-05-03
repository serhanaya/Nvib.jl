__precompile__()

module Nvib

include("analysis_type.jl")
include("frequency_table.jl")
include("coords.jl")
include("fundamental_matrix.jl")
include("detrmt.jl")







using Pkg
Pkg.activate(joinpath(homedir(), ".julia/environments/v1.7"))
Pkg.add("ProgressMeter")
end
using ProgressMeter
export IPBeam, OPBeam, freqTab

end 