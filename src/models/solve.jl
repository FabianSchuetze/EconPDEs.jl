##############################################################################
##
## A EconPDE model must define:
## A function pde
##
## See examples in the folder /models
##
##############################################################################

abstract EconPDEModel

function hjb!{N1}(apm::EconPDEModel, grid::StateGrid{N1}, y::Array, ydot::Array)
    N = div(length(y), prod(size(grid)))
    colontuple = ([Colon() for i in 1:N1]...)
    fy = ([ReflectingArray(view(y, colontuple..., i)) for i in 1:N]...)
    hjb!(apm, grid, fy, ydot)
end

function hjb!{N1, N, T}(apm::EconPDEModel, grid::StateGrid{N1}, fy::NTuple{N, T}, ydot::Array)
    for i in eachindex(grid)
        outi, drifti, othersi = pde(apm, grid, fy, i)
        outi, drifti, othersi = pde(apm, grid, fy, i, drifti)
        _setindex!(ydot, outi, i)
    end
    return ydot
end

@generated function _setindex!{N, T}(ydot, outi::NTuple{N, T}, i)
    Expr(:block, [:(setindex!(ydot, outi[$k], i, $k)) for k in 1:N]...)
end
_setindex!(ydot, outi, i) = setindex!(ydot, outi, i)

function solve(apm::EconPDEModel, grid::StateGrid, y0; kwargs...)
    Ψtc((y, ydot) -> hjb!(apm, grid, y, ydot), y0; kwargs...)
end
##############################################################################
##
## Derive default
##
##############################################################################


function derive(grid::StateGrid, y::ReflectingArray, ituple::CartesianIndex{1}, drift = (0.0,))
  is = ituple[1]
  μs = drift[1]
  Δs, = grid.Δx
  Δsm, = grid.Δxm
  Δsp, = grid.Δxp
  p = y[is]
  if μs >= 0.0
      ps = (y[is + 1] - y[is]) / Δsp[is]
  else
      ps = (y[is] - y[is - 1]) / Δsm[is]
  end
  pss = (Δsm[is] * y[is + 1] + Δsp[is] * y[is - 1] - 2 * Δs[is] * y[is]) / (Δs[is] * Δsm[is] * Δsp[is])
  return p, ps, pss
end

function derive(grid::StateGrid, y::ReflectingArray, ituple::CartesianIndex{2}, drift = (0.0, 0.0))
    iμ, iσ = ituple[1], ituple[2]
    ix, iν = ituple[1], ituple[2]
    μX, μν = drift
    Δx, Δν = grid.Δx
    if μX <= 0.0
      indx1 = 0
      indx2 = -1
    else
     indx1 = 1
     indx2 = 0
    end
    if μν <= 0.0
      indν1 = 0
      indν2 = -1
    else
      indν1 = 1
      indν2 = 0
    end
    p = y[ix, iν]
    px = (y[ix + indx1, iν] - y[ix + indx2, iν]) / Δx[ix]
    pν = (y[ix, iν + indν1] - y[ix, iν + indν2]) / Δν[iν]
    pxx = (y[ix + 1, iν] + y[ix - 1, iν] - 2 * y[ix, iν]) / Δx[ix]^2
    pνν = (y[ix, iν + 1] + y[ix, iν - 1] - 2 * y[ix, iν]) / Δν[iν]^2
    pxν = (y[ix + indx1, iν + indν1] - y[ix + indx1, iν + indν2] - y[ix + indx2, iν + indν1] + y[ix + indx2, iν + indν2]) / (Δν[iν] * Δx[ix])
    return p, px, pν, pxx, pxν, pνν
end
##############################################################################
##
## Full solve (i.e. comptue all quantities)
## 
##############################################################################

function get_info{N, T}(apm, grid, fy::NTuple{N, T})
    i = start(eachindex(grid))
    outi, drifti, othersi = pde(apm, grid, fy, i)
    collect(map(first, othersi)), collect(map(typeof, map(last, othersi)))
end

function compute_arrays{N1}(apm, grid::StateGrid{N1}, y)
    N = div(length(y), prod(size(grid)))
    colontuple = ([Colon() for i in 1:N1]...)
    fy = ([ReflectingArray(view(y, colontuple..., i)) for i in 1:N]...)
    names, types = get_info(apm, grid, fy)
    len = length(names)
    A = Dict([Pair(names[i] => Array(types[i], size(grid))) for i in 1:length(names)])
    for i in eachindex(grid)
        outi, drifti, othersi = pde(apm, grid, fy, i)
        outi, drifti, othersi = pde(apm, grid, fy, i, drifti)
        for j in 1:len
            A[first(othersi[j])][i] = last(othersi[j])
        end
    end
    return A
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

function simulate(grid, a, shocks; dt = 1 / 12, x0 = grid.x[1][rand(Categorical(stationary_distribution(grid, a)), size(shocks, 2))])
    # interpolate all functions
    ai = Dict([Pair(k => interpolate(grid.x, a[k], Gridded(Linear()))) for k in keys(a)])
    aT = Dict([Pair(k => zeros(shocks)) for k in keys(a)])
    aT[:id] = zeros(shocks)
    aT[:t] = zeros(shocks)
    aT[:shock] = zeros(shocks)
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
            xt = xt + ai[:μx][xt] * dt + ai[:σx][xt] * shocks[t, id] * sqrtdt
        end
    end
    return aT
end
