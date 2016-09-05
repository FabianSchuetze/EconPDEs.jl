
##############################################################################
##
## Reflecting Array (i.e. 0 is 1 and N+1 is N)
##
##############################################################################

type ReflectingArray{T, N}
    A::Array{T, N}
end

Base.size(y::ReflectingArray, args...) = size(y.A, args...)
Base.eltype(y::ReflectingArray) = eltype(y.A)
Base.eachindex(y) = eachindex(y.A)
@generated function Base.getindex{T, N}(A::ReflectingArray{T, N}, args...)
    Expr(:call, :getindex, :(A.A), tuple([_helper(args[i], i) for i in 1:N]...))
end
@generated function Base.setindex!{T, N}(A::ReflectingArray{T, N}, value, args...)
    Expr(:call, :setindex!, :(A.A), :value, tuple([_helper(args[i], i) for i in 1:N]...))
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

