module doubleLogistics
import Libdl
# using Printf
using Parameters


include("curveFitting/FitDL.jl")
include("optim/optim_pheno.jl")
include("weights/wFUN.jl")

end # module
