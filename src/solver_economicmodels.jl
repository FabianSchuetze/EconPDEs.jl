
#========================================================================================

Type Reflecting Array (i.e. 0 is 1 and N+1 is N)

========================================================================================#
type ReflectingArray{P, T, N} <: AbstractArray{T, N}
    A::P
end
ReflectingArray{T, N}(A::AbstractArray{T, N}) = ReflectingArray{typeof(A), T, N}(A)
Base.size(y::ReflectingArray, args...) = size(y.A, args...)
Base.eltype{P, T, N}(y::ReflectingArray{P, T, N}) = T
Base.eachindex(y) = eachindex(y.A)
@generated function Base.getindex{P, T, N}(A::ReflectingArray{P, T, N}, args...)
    Expr(:call, getindex, :(A.A), [_helper(args[i], i) for i in 1:N]...)
end
@generated function Base.setindex!{P, T, N}(A::ReflectingArray{P, T, N}, value, args...)
    Expr(:call, :setindex!, :(A.A), :value, [_helper(args[i], i) for i in 1:N]...)
end
_helper(x::Type{Int}, i) = :(clamp(args[$i], 1, size(A.A, $i)))
_helper(x::Type{Colon}, i) = :(:)



#========================================================================================

Type State Grid

========================================================================================#

type StateGrid{N}
    x::NTuple{N, Vector{Float64}}
    Δx::NTuple{N, Vector{Float64}}
    Δxm::NTuple{N, Vector{Float64}}
    Δxp::NTuple{N, Vector{Float64}}
    name::NTuple{N, Symbol}
end

function make_Δ(x)
    n = length(x)
    Δxm = zeros(x)
    for i in 2:n
        Δxm[i] = x[i] - x[i-1]
    end
    Δxp = zeros(x)
    for i in 1:(n-1)
        Δxp[i] = x[i+1] - x[i]
    end
    Δx = zeros(x)
    Δx[1] = Δxp[1]
    Δxm[1] = Δxp[1]
    Δx[end] = Δxm[end]
    Δxp[end] = Δxm[end]
    for i in 2:(n-1)
        Δx[i] = 0.5 * (Δxm[i] + Δxp[i])
    end
    return Δx, Δxm, Δxp
end

function StateGrid(; kwargs...)
    x = Vector{Vector{Float64}}()
    Δx = Vector{Vector{Float64}}()
    Δxm = Vector{Vector{Float64}}()
    Δxp = Vector{Vector{Float64}}()
    name = Vector{Symbol}()
    for (k, v) in kwargs
        push!(x, v)
        push!(name, k)
        Δv, Δvm, Δvp = make_Δ(v)
        push!(Δx, Δv)
        push!(Δxm, Δvm)
        push!(Δxp, Δvp)
    end
    return StateGrid{length(x)}(tuple(x...), tuple(Δx...), tuple(Δxm...), tuple(Δxp...), tuple(name...))
end

Base.eachindex(grid::StateGrid) = CartesianRange(size(grid))

@generated function Base.size{N}(grid::StateGrid{N}, args...)
    Expr(:call, :tuple, [:(length(grid.x[$i])) for i in 1:N]...)
end

function Base.getindex{N}(grid::StateGrid{N}, args...)
    getindex(grid, CartesianIndex{N}(args...))
end

@generated function Base.getindex{N}(grid::StateGrid{N}, args::CartesianIndex)
    Expr(:call, :tuple, [:(getindex(grid.x[$i], args[$i])) for i in 1:N]...)
end

function Base.getindex{N}(grid::StateGrid{N}, x::Symbol)
    grid.x[find(collect(grid.name) .== x)[1]]
end

#========================================================================================

Derive

========================================================================================#

# Case with 1 state variable
function derive(grid::StateGrid{1}, y::ReflectingArray, ituple::NTuple{1, Int}, drift = (0.0,))
  i = ituple[1]
  μx = drift[1]
  Δx, = grid.Δx
  Δxm, = grid.Δxm
  Δxp, = grid.Δxp
  p = y[i]
  if μx >= 0.0
      px = (y[i + 1] - y[i]) / Δxp[i]
  else
      px = (y[i] - y[i - 1]) / Δxm[i]
  end
  pxx = (Δxm[i] * y[i + 1] + Δxp[i] * y[i - 1] - 2 * Δx[i] * y[i]) / (Δx[i] * Δxm[i] * Δxp[i])
  return p, px, pxx
end

# Case with 2 state variables
function derive(grid::StateGrid{2}, y::ReflectingArray, ituple::NTuple{2, Int}, drift = (0.0, 0.0))
    i1, i2 = ituple[1], ituple[2]
    μx1, μx2 = drift
    Δx1, Δx2 = grid.Δx
    if μx1 >= 0.0
        i1h = i1 + 1
        i1l = i1
    else
      i1h = i1
      i1l = i1 - 1
    end
    if μx2 >= 0.0
        i2h = i2 + 1
        i2l = i2
    else
      i2h = i2
      i2l = i2 - 1
    end
    p = y[i1, i2]
    px1 = (y[i1h, i2] - y[i1l, i2]) / Δx1[i1]
    px2 = (y[i1, i2h] - y[i1, i2l]) / Δx2[i2]
    px1x1 = (y[i1 + 1, i2] + y[i1 - 1, i2] - 2 * y[i1, i2]) / Δx1[i1]^2
    px2x2 = (y[i1, i2 + 1] + y[i1, i2 - 1] - 2 * y[i1, i2]) / Δx2[i2]^2
    px1x2 = (y[i1h, i2h] - y[i1h, i2l] - y[i1l, i2h] + y[i1l, i2l]) / (Δx1[i1] * Δx2[i2])
    return p, px1, px2, px1x1, px1x2, px2x2
end

#========================================================================================

Type EconPDEModel

========================================================================================#

abstract EconPDEModel

function hjb!{N1}(apm::EconPDEModel, grid::StateGrid{N1}, y::Array, ydot::Array)
    N = div(length(y), prod(size(grid)))
    colontuple = ([Colon() for i in 1:N1]...)
    fy = ([ReflectingArray(view(y, colontuple..., i)) for i in 1:N]...)
    hjb!(apm, grid, fy, ydot)
end

function hjb!{N1, N, T}(apm::EconPDEModel, grid::StateGrid{N1}, fy::NTuple{N, T}, ydot::Array)
    for i in eachindex(grid)
        outi, drifti, othersi = pde(apm, grid, fy, i.I)
        outi, drifti, othersi = pde(apm, grid, fy, i.I, drifti)
        _setindex!(ydot, outi, i)
    end
    return ydot
end

@generated function _setindex!{N, T}(ydot, outi::NTuple{N, T}, i)
    Expr(:block, [:(setindex!(ydot, outi[$k], i, $k)) for k in 1:N]...)
end

_setindex!(ydot, outi, i) = setindex!(ydot, outi, i)


#========================================================================================

Solve

========================================================================================#

function solve(apm::EconPDEModel, grid::StateGrid, y0; kwargs...)
    y, distance = Ψtc((y, ydot) -> hjb!(apm, grid, y, ydot), y0; kwargs...)
    a = create_dictionary(apm, grid, y)
    return a, distance
end

function create_dictionary{N1}(apm, grid::StateGrid{N1}, y)
    N = div(length(y), prod(size(grid)))
    colontuple = ([Colon() for i in 1:N1]...)
    fy = ([ReflectingArray(view(y, colontuple..., i)) for i in 1:N]...)
    names, types = get_namestypes(apm, grid, fy)
    len = length(names)
    A = Dict([Pair(names[i] => Array(types[i], size(grid))) for i in 1:length(names)])
    for i in eachindex(grid)
        outi, drifti, othersi = pde(apm, grid, fy, i.I)
        outi, drifti, othersi = pde(apm, grid, fy, i.I, drifti)
        for j in 1:len
            A[first(othersi[j])][i] = last(othersi[j])
        end
    end
    return A
end

function get_namestypes(apm, grid, fy::NTuple)
    i = start(eachindex(grid))
    outi, drifti, othersi = pde(apm, grid, fy, i.I)
    collect(map(first, othersi)), collect(map(typeof, map(last, othersi)))
end

#========================================================================================

Stationary Distribution

========================================================================================#

# Case with 1 state variable
function stationary_distribution(grid::StateGrid{1}, a)
    n, = size(grid)
    Δx, = grid.Δx
    Δxp, = grid.Δxp
    Δxm, = grid.Δxm
    A = ReflectingArray(zeros(n, n))
    for i in 1:n
        μ = a[Symbol(:μ, grid.name[1])][i]
        σ2 = a[Symbol(:σ, grid.name[1])][i]^2
        if μ >= 0
            A[i + 1, i] += μ / Δxp[i]
            A[i, i] -= μ / Δxp[i]
        else
            A[i, i] += μ / Δxm[i] 
            A[i - 1, i] -= μ / Δxm[i] 
        end
        A[i - 1, i] += 0.5 * σ2 * Δxp[i] / (Δx[i] * Δxm[i] * Δxp[i])
        A[i, i] -= 0.5 * σ2 * 2 * Δx[i] / (Δx[i] * Δxm[i] * Δxp[i])
        A[i + 1, i] += 0.5 * σ2 * Δxm[i] / (Δx[i] * Δxm[i] * Δxp[i])
    end
    for j in 1:size(A, 2)
        A[1, j] = 1.0
    end
    b = vcat(1.0, zeros(n - 1))
    density = A.A \ b
    @assert all(density .> -1e-5)
    density = abs(density) ./ sumabs(density)  
    return density 
end

# Case with 2 state variables
function stationary_distribution(grid::StateGrid{2}, a)
    n1, n2 = size(grid)
    A = ReflectingArray(zeros(n1, n2, n1, n2))
    for i2 in 1:n2
        for i1 in 1:n1
            μ = a[Symbol(:μ, grid.name[1])][i1, i2]
            σ2 = a[Symbol(:σ, grid.name[1], :2)][i1, i2]
            if μ >= 0
                i1h = i1 + 1
                i1l = i1
            else
               i1h = i1
               i1l = i1 - 1
            end
            A[i1h, i2, i1, i2] += μ /  grid.Δx[1][i1]
            A[i1l, i2, i1, i2] -= μ /  grid.Δx[1][i1]
            A[i1 - 1, i2, i1, i2] += 0.5 * σ2 /  grid.Δx[1][i1]^2
            A[i1, i2, i1, i2] -= 0.5 * σ2 * 2/  grid.Δx[1][i1]^2
            A[i1 + 1,i2, i1, i2] += 0.5 * σ2 /  grid.Δx[1][i1]^2

            μ = a[Symbol(:μ, grid.name[2])][i1, i2]
            σ2 = a[Symbol(:σ, grid.name[2], :2)][i1, i2]
            if μ >= 0
                i2h = i2 + 1
                i2l = i2
            else
               i2h = i2
               i2l = i2 - 1
            end
            A[i1, i2h, i1, i2] += μ /  grid.Δx[1][i1]
            A[i1, i2l, i1, i2] -= μ /  grid.Δx[1][i1]
            A[i1, i2 - 1, i1, i2] += 0.5 * σ2 /  grid.Δx[2][i2]^2
            A[i1, i2, i1, i2] -= 0.5 * σ2 * 2/  grid.Δx[2][i2]^2
            A[i1,i2 + 1, i1, i2] += 0.5 * σ2 /  grid.Δx[2][i2]^2

            σ12 = a[Symbol(:σ, grid.name[1], :σ, grid.name[2])][i1, i2]
            A[i1h, i2h, i1, i2] += σ12 /  (grid.Δx[1][i1] * grid.Δx[2][i2])
            A[i1l, i2h, i1, i2] -= σ12 /  (grid.Δx[1][i1] * grid.Δx[2][i2])
            A[i1h, i2l, i1, i2] -= σ12 /  (grid.Δx[1][i1] * grid.Δx[2][i2])
            A[i1l, i2l, i1, i2] += σ12 /  (grid.Δx[1][i1] * grid.Δx[2][i2])
        end
    end
    A = reshape(A.A, (n1 * n2, n1 * n2))
    for j in 1:size(A, 2)
        A[1, j] = 1.0
    end
    b = vcat(1.0, zeros(size(A, 2) - 1))
    density = A \ b
    density = abs(density) ./ sumabs(density)  
    return reshape(density, (n1, n2))
end

#========================================================================================

Simulate

========================================================================================#


function simulate{N}(grid::StateGrid{N}, a, shocks::Dict; dt = 1 / 12, x0 = nothing)
    T = size(shocks[first(keys(shocks))], 1)
    I = size(shocks[first(keys(shocks))], 2)
    if x0 == nothing
        i0 = rand(Categorical(vec(stationary_distribution(grid, a))), I)
        if N == 1
            x0 = Dict(grid.name[1] => grid.x[1][i0])
        elseif N == 2
            i10 = mod(i0 - 1, length(grid.x[1])) + 1
            i20 = div(i0 - 1, length(grid.x[1])) + 1
            x0 = Dict(grid.name[1] => grid.x[1][i10], grid.name[2] => grid.x[2][i20])
        end
    end
    # interpolate all functions
    ai = Dict([Pair(k => interpolate(grid.x, a[k], Gridded(Linear()))) for k in keys(a)])
    aT = Dict([Pair(k => zeros(T, I)) for k in keys(a)])
    aT[:id] = zeros(T, I)
    aT[:t] = zeros(T, I)
    for k in keys(shocks)
        aT[k] = zeros(T, I)
    end
    sqrtdt = sqrt(dt)
    for id in 1:I
        xt = tuple([x0[grid.name[i]][id] for i in 1:N]...)
        for t in 1:T
            for k in keys(a)
                aT[k][t, id] = ai[k][xt...]
            end
            aT[:id][t, id] = id
            aT[:t][t, id] = t
            for k in keys(shocks)
                aT[k][t, id] = shocks[k][t, id]
            end
            xt = tuple([xt[i] + _update_state(xt, grid.name[i], shocks, ai, t, id, dt, sqrtdt) for i in 1:N]...)
        end
    end
    return aT
end

function _update_state(xt, name, shocks, ai, t, id, dt, sqrtdt)
    out = ai[Symbol(:μ, name)][xt...] * dt
    if length(keys(shocks)) == 1
        out += ai[Symbol(:σ, name)][xt...] * shocks[first(keys(shocks))][t, id] * sqrtdt
    else
        for k in keys(shocks)
            out += ai[Symbol(:σ, name, :_, k)][xt...] * shocks[k][t, id] * sqrtdt
        end
    end
    return out
end

