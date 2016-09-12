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

function hjb!{N1}(apm::EconPDEModel, grid::StateGrid{N1}, y::Array, ydot::Array)
    N = div(length(y), prod(size(grid)))
    colontuple = ((Colon() for i in 1:N1)...)
    fy = ((ReflectingArray(view(y, colontuple..., i)) for i in 1:N)...)
    hjb!(apm, grid, fy, ydot)
end

@generated function hjb!{N1, N, T}(apm::EconPDEModel, grid::StateGrid{N1}, fy::NTuple{N, T}, ydot::Array)
    if N == 1
        exp1 = :(derive(apm, grid, fy[1], i))
        exp2 = :(derive(apm, grid, fy[1], i, drifti))
    else
        exp1 = Expr(:call, :tuple, (:(derive(apm, grid, fy[$k], i)) for k in 1:N)...)
        exp2 = Expr(:call, :tuple, (:(derive(apm, grid, fy[$k], i, drifti)) for k in 1:N)...)
    end
    exp3 = Expr(:block, (:(setindex!(ydot, outi[$k], i, $k)) for k in 1:N)...)
    quote
        for i in eachindex(grid)
            outi, drifti, othersi = pde(apm, grid[i], simplify($exp1))
            outi, drifti, othersi = pde(apm, grid[i], simplify($exp2))
            $exp3
        end
        return ydot
    end
end

simplify{T}(x::NTuple{1, T}) = x[1]
simplify{N, T}(x::NTuple{N, T}) = x


function get_names{N, T}(apm, grid, fy::NTuple{N, T})
    i = start(eachindex(grid))
    functionsi = simplify(((derive(apm, grid, fy[k], i) for k in 1:N)...))
    outi, drifti, othersi = pde(apm, grid[i], functionsi)
    collect(map(first, othersi))
end

function compute_arrays{N1}(apm, grid::StateGrid{N1}, y)
    N = div(length(y), prod(size(grid)))
    colontuple = ((Colon() for i in 1:N1)...)
    fy = ((ReflectingArray(view(y, colontuple..., i)) for i in 1:N)...)
    names = get_names(apm, grid, fy)
    len = length(names)
    A = Dict([Pair(name => zeros(size(grid))) for name in names])
    for i in eachindex(grid)
        functionsi = simplify(((derive(apm, grid, fy[k], i) for k in 1:N)...))
        outi, drifti, othersi = pde(apm, grid[i], functionsi)
        functionsi = simplify(((derive(apm, grid, fy[k], i) for k in 1:N)...))
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

function simulate(grid, a, shocks; dt = 1.0, x0 = grid.x[1][rand(Categorical(stationary_distribution(grid, a)), size(shocks, 2))])
    if length(grid.name) > 2
        throw("simulate does not work with multiple state variables")
    end
    # interpolate all functions
    ai = Dict([Pair(k => interpolate(grid.x, a[k], Gridded(Linear()))) for k in keys(a)])
    aT = Dict([Pair(k => zeros(shocks)) for k in keys(a)])
    aT[:id] = zeros(shocks)
    aT[:t] = zeros(shocks)
    aT[:shock] = zeros(shocks)
    μx = Symbol(:μ, grid.name[1])
    σx = Symbol(:σ, grid.name[1])
    sqrtdt = sqrt(dt)
    for id in 1:size(shocks, 2)
        xt = x0[id]
        for t in 1:size(shocks, 1)
            for k in keys(a)
                aT[k][t, id] = ai[k][xt]
            end
            aT[:id][t, id] = id
            aT[:t][t, id] = t
            aT[:shock][t, id] = shocks[t, id]
            xt = xt + ai[μx][xt] * dt + ai[σx][xt] * shocks[t, id] * sqrtdt
        end
    end
    return aT
end
