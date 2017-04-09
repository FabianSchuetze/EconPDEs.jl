module EconPDEs


import NLsolve: nlsolve
import Distributions: Normal, Gamma, Categorical
using Interpolations
using Compat
import Compat.view
import Combinatorics: with_replacement_combinations
import NamedTuples: @NT, NamedTuple
##############################################################################
##
## Load files
##
##############################################################################
include("nl_solve.jl")
include("pde_solve.jl")
include("utils.jl")
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
export nl_solve,
pde_solve,
EconPDEModel,
state_grid,
initialize, 
pde,
@NT,
simulate,
GarleanuPanageasModel,
CampbellCochraneModel,
BansalYaronModel,
DiTellaModel,
WangWangYangModel
end