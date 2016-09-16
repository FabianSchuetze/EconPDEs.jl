
##############################################################################
##
## Reflecting Array (i.e. 0 is 1 and N+1 is N)
##
##############################################################################
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

##############################################################################
##
## State Grid
##
##############################################################################

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

##############################################################################
##
## Derive default
##
##############################################################################

function derive(grid::StateGrid, y::ReflectingArray, ituple::CartesianIndex{1}, drift = (0.0,))
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

function derive(grid::StateGrid, y::ReflectingArray, ituple::CartesianIndex{2}, drift = (0.0, 0.0))
    i1, i2 = ituple[1], ituple[2]
    μx1, μx2 = drift
    Δx1, Δx2 = grid.Δx
    if μx1 <= 0.0
      i1h = i1
      i1l = i1 - 1
    else
     i1h = i1 + 1
     i1l = i1
    end
    if μx2 <= 0.0
      i2h = i2
      i2l = i2 - 1
    else
      i2h = i2 + 1
      i2l = i2
    end
    p = y[i1, i2]
    px1 = (y[i1h, i2] - y[i1l, i2]) / Δx1[i1]
    px2 = (y[i1, i2h] - y[i1, i2l]) / Δx2[i2]
    px1x1 = (y[i1 + 1, i2] + y[i1 - 1, i2] - 2 * y[i1, i2]) / Δx1[i1]^2
    px2x2 = (y[i1, i2 + 1] + y[i1, i2 - 1] - 2 * y[i1, i2]) / Δx2[i2]^2
    px1x2 = (y[i1h, i2h] - y[i1h, i2l] - y[i1l, i2h] + y[i1l, i2l]) / (Δx1[i1] * Δx2[i2])
    return p, px1, px2, px1x1, px1x2, px2x2
end