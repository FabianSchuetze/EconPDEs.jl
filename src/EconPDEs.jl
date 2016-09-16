module EconPDEs


import NLsolve: nlsolve
import Distributions: Normal, Gamma, Categorical
using Interpolations
using Compat
import Compat.view
##############################################################################
##
## Load files
##
##############################################################################
include("Ψtc.jl")
include("models/utils.jl")
include("models/solve.jl")
include("models/CampbellCochrane.jl")
include("models/BansalYaron.jl")
include("models/GarleanuPanageas.jl")
include("models/DiTella.jl")
include("models/WangWangYang.jl")

##############################################################################
##
## Exported methods and types 
##
##############################################################################
export Ψtc,
hjb!,
StateGrid,
ReflectingArray,
EconPDEModel,
initialize, 
solve,
fullsolve,
GarleanuPanageasModel,
CampbellCochraneModel,
BansalYaronModel,
DiTellaModel,
WangWangYangModel
end