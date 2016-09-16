using Optim
function rosenbrock_res(x, r)
    r[1] = 10.0 * (x[2] - x[1]^2 )
    r[2] =  1.0 - x[1]
    return r
end

function rosenbrock_jac(x, j)
    j[1, 1] = -20.0 * x[1]
    j[1, 2] =  10.0
    j[2, 1] =  -1.0
    j[2, 2] =   0.0
    return j
end

r = zeros(2)
j = zeros(2,2)

frb(x) = rosenbrock_res(x, r)
grb(x) = rosenbrock_jac(x, j)


result = Optim.levenberg_marquardt(frb, grb, [150.0, 150.0]; lower = [10.0, 10.0], upper = [200.0, 200.0])

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

