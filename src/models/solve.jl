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
    fy = ReflectingArray(y)
    for i in eachindex(grid)
        functionsi = derive(apm, grid, fy, i)
        outi, drifti, othersi = pde(apm, grid[i], functionsi)
        functionsi = derive(apm, grid, fy, i, drifti)
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
    A = Dict([Pair(name => zeros(size(grid))) for name in names])
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

#========================================================================================

Minimal function to solve for stationary distribution and simulate in the case of one state variable
TO DO : Version with several state variables

========================================================================================#
function stationary_distribution(grid, a)
    if length(grid.name) > 2
        throw("simulate does not work with multiple state variables")
    end
    xn, = size(grid)
    Δx, = grid.Δx
    Δxp, = grid.Δxp
    Δxm, = grid.Δxm
    A = ReflectingArray(zeros(xn, xn))
    μname = Symbol(:μ, grid.name[1])
    σname = Symbol(:σ, grid.name[1])
    for i in 1:xn
        if a[μname][i] >= 0
            A[i + 1, i] += a[μname][i] / Δxp[i]
            A[i, i] -= a[μname][i] / Δxp[i]
        else
            A[i, i] += a[μname][i] / Δxm[i] 
            A[i - 1, i] -= a[μname][i] / Δxm[i] 
        end
        A[i - 1, i] += 0.5 * a[σname][i]^2 * Δxp[i] / (Δx[i] * Δxm[i] * Δxp[i])
        A[i, i] -= 0.5 * a[σname][i]^2 * 2 * Δx[i] / (Δx[i] * Δxm[i] * Δxp[i])
        A[i + 1, i] += 0.5 * a[σname][i]^2 * Δxm[i] / (Δx[i] * Δxm[i] * Δxp[i])
    end
    for j in 1:xn
        A[1, j] = 1.0
    end
    b = vcat(1.0, zeros(xn - 1))
    density = A.A \ b
    @assert all(density .> -1e-5)
    density = abs(density) ./ sumabs(density)  
    return density 
end

function simulate(grid, a, shocks; dt = 1.0, x0 = sum(grid.x[1] .* stationary_distribution(grid, a)))
    if length(grid.name) > 2
        throw("simulate does not work with multiple state variables")
    end
    # interpolate all functions
    @compat ai = Dict([Pair(k => interpolate(grid.x, a[k], Gridded(Linear()))) for k in keys(a)])
    y = zeros(shocks)
    @compat aT = Dict([Pair(k => zeros(shocks)) for k in keys(a)])

    μname = Symbol(:μ, grid.name[1])
    σname = Symbol(:σ, grid.name[1])
    name = grid.name[1]
    sqrtdt = sqrt(dt)
    xt = x0
    for id in 1:size(shocks, 2)
        for t in 1:size(shocks, 1)
            for k in keys(a)
                aT[k][t, id] = ai[k][xt]
            end
            xt = xt + ai[μname][xt] * dt + ai[σname][xt] * shocks[t] * sqrtdt
        end
    end
    return aT
end
