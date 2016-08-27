##############################################################################
##
## A EconPDE model must define:
## 1. A function derive
## 2. A function pde
##
## See examples in the folder /models
##
##############################################################################
abstract EconPDEModel

function hjb!(apm::EconPDEModel, grid::StateGrid, y::Array, ydot::Array)
    n_functions = div(length(y), prod(size(grid)))
    y = ReflectingArray(y)
    for i in eachindex(grid)
        functionsi = derive(apm, grid, y, i)
        outi, drifti, othersi = pde(apm, grid[i], functionsi)
        functionsi = derive(apm, grid, y, i, drifti)
        outi, drifti, othersi = pde(apm, grid[i], functionsi)
        for k in 1:n_functions
            ydot[i, k] = outi[k]
        end
    end
    return ydot
end

function get_names(apm, grid, y)
    y = ReflectingArray(y)
    i = start(eachindex(grid))
    functionsi = derive(apm, grid, y, i)
    outi, drifti, othersi = pde(apm, grid[i], functionsi)
    collect(map(first, othersi))
end

function compute_arrays(apm, grid, y)
    n_functions = div(length(y), prod(size(grid)))
    y = ReflectingArray(y)
    names = get_names(apm, grid, y)
    len = length(names)
    A = Dict(name => zeros(size(grid)) for name in names)
    for i in eachindex(grid)
        functionsi = derive(apm, grid, y, i)
        outi, drifti, othersi = pde(apm, grid[i], functionsi)
        functionsi = derive(apm, grid, y, i, drifti)
        outi, drifti, othersi = pde(apm, grid[i], functionsi)
        for j in 1:len
            A[first(othersi[j])][i] = last(othersi[j])
        end
    end
    return A
end


function solve(apm::EconPDEModel, grid::StateGrid, y0; kwargs...)
    Ψtc((y, ydot) -> hjb!(apm, grid, y, ydot), y0; kwargs...)
end

function fullsolve(apm::EconPDEModel, grid::StateGrid, y0; kwargs...)
    y, distance = solve(apm, grid, y0; kwargs...)
    a = compute_arrays(apm, grid, y)
    return a, distance
end
