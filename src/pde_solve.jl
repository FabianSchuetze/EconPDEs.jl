##NamedTuple create type in the module
# at parsing type
# 1. creates a type by looking at expression, so something like T_NT_a_b{T1, T2}
# 2. replace @NT by T_NT_ab



#========================================================================================

Type State Grid

========================================================================================#
type StateGrid{N}
    x::NTuple{N, Vector{Float64}}
    invΔx::NTuple{N, Vector{Float64}}
    invΔxm::NTuple{N, Vector{Float64}}
    invΔxp::NTuple{N, Vector{Float64}}
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
    return 1./Δx, 1./Δxm, 1./Δxp
end

function StateGrid(x)
    names = tuple((k for k in keys(x))...)
    x = tuple(x...)
    StateGrid{length(x)}(
        x,
        map(x -> make_Δ(x)[1], x),
        map(x -> make_Δ(x)[2], x),
        map(x -> make_Δ(x)[3], x),
        names
    )
end

Base.eachindex(grid::StateGrid) = CartesianRange(size(grid))
Base.size{N}(grid::StateGrid{N}, args...) = map(length, grid.x)
Base.ndims(grid) = length(size(grid))
@generated function Base.getindex{T, N}(grid::StateGrid{N}, ::Type{T}, args::CartesianIndex)
    quote
        $(Expr(:meta, :inline))
        $(Expr(:call, T, [:(getindex(grid.x[$i], args[$i])) for i in 1:N]...))
    end
end
Base.getindex{N}(grid::StateGrid{N}, x::Symbol) = grid.x[find(collect(grid.name) .== x)[1]]


#========================================================================================

Type Reflecting Array (i.e. for a vector A[0] is A[1] and A[N+1] is A[N])

========================================================================================#
type ReflectingArray{P, T, N} <: AbstractArray{T, N}
    A::P
end
ReflectingArray{T, N}(A::AbstractArray{T, N}) = ReflectingArray{typeof(A), T, N}(A)
Base.size(y::ReflectingArray, args...) = size(y.A, args...)
Base.eltype{P, T, N}(y::ReflectingArray{P, T, N}) = T
Base.eachindex(y) = eachindex(y.A)
@generated function Base.getindex{P, T, N}(A::ReflectingArray{P, T, N}, args...)
    quote
        $(Expr(:meta, :inline))
        $(Expr(:call, getindex, :(A.A), [_helper(args[i], i) for i in 1:N]...))
    end
end
@generated function Base.setindex!{P, T, N}(A::ReflectingArray{P, T, N}, value, args...)
    Expr(:call, :setindex!, :(A.A), :value, [_helper(args[i], i) for i in 1:N]...)
end
_helper(x::Type{Int}, i) = :(clamp(args[$i], 1, size(A.A, $i)))
_helper(x::Type{Colon}, i) = :(:)


# Case with 1 state variable
@generated function derive{Tsolution}(::Type{Tsolution}, grid::StateGrid{1}, y::ReflectingArray, icar, drift = (0.0,))
    N = div(nfields(Tsolution), 3)
    expr = Expr[]
    for k in 1:N
        push!(expr, :(y[i, $k]))
        push!(expr, :(μx >= 0.0 ? (y[i + 1, $k] - y[i, $k]) * invΔxp[i] : (y[i, $k] - y[i - 1, $k]) * invΔxm[i]))
        push!(expr, :(y[i + 1, $k] * invΔxp[i] * invΔx[i] + y[i - 1, $k] * invΔxm[i] * invΔx[i] - 2 * y[i, $k] * invΔxp[i] * invΔxm[i]))
    end
    out = Expr(:call, Tsolution, expr...)
    quote
        $(Expr(:meta, :inline))
        i = icar[1]
        μx = drift[1]
        invΔx = grid.invΔx[1]
        invΔxm = grid.invΔxm[1]
        invΔxp = grid.invΔxp[1]
        $out
    end
end

# Case with 2 state variables
@generated function derive{Tsolution}(::Type{Tsolution}, grid::StateGrid{2}, y::ReflectingArray, icar, drift = (0.0, 0.0))
    N = div(nfields(Tsolution), 6)
    expr = Expr[]
    for k in 1:N
        push!(expr, :(y[i1, i2, $k]))
        push!(expr, :((y[i1h, i2, $k] - y[i1l, i2, $k]) * invΔx1[i1]))
        push!(expr, :((y[i1, i2h, $k] - y[i1, i2l, $k]) * invΔx2[i2]))
        push!(expr, :((y[i1 + 1, i2, $k] + y[i1 - 1, i2, $k] - 2 * y[i1, i2, $k]) * invΔx1[i1]^2))
        push!(expr, :((y[i1h, i2h, $k] - y[i1h, i2l, $k] - y[i1l, i2h, $k] + y[i1l, i2l, $k]) * invΔx1[i1] * invΔx2[i2]))
        push!(expr, :((y[i1, i2 + 1, $k] + y[i1, i2 - 1, $k] - 2 * y[i1, i2, $k]) * invΔx2[i2]^2))
    end
    out = Expr(:call, Tsolution, expr...)
    quote
        $(Expr(:meta, :inline))
        i1, i2 = icar[1], icar[2]
        μx1, μx2 = drift[1], drift[2]
        invΔx1, invΔx2 = grid.invΔx[1], grid.invΔx[2]
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
        $out
    end
end


#========================================================================================

Type EconPDEModel

========================================================================================#

abstract EconPDEModel


function hjb!{Ngrid, Tstate, Tsolution}(apm, grid::StateGrid{Ngrid}, ::Type{Tstate}, ::Type{Tsolution}, y, ydot)
    y = ReflectingArray(y)
    for i in eachindex(grid)
        state = getindex(grid, Tstate, i)
        solution = derive(Tsolution, grid, y, i)
        drifti = pde(apm, state, solution)[2]
        #upwind
        solution = derive(Tsolution, grid, y, i, drifti)
        outi = pde(apm, state, solution)[1]
        _setindex!(ydot, outi, i)
    end
    return ydot
end

@generated function _setindex!{N, T}(ydot, outi::NTuple{N, T}, i)
    quote
        $(Expr(:meta, :inline))
        $(Expr(:block, [:(setindex!(ydot, outi[$k], i, $k)) for k in 1:N]...))
    end
end
_setindex!(ydot, outi, i) = setindex!(ydot, outi, i)


function create_dictionary{Ngrid, Tstate, Tsolution}(apm, grid::StateGrid{Ngrid}, ::Type{Tstate}, ::Type{Tsolution}, y)
    y = ReflectingArray(y)
    i0 = start(eachindex(grid))
    state = getindex(grid, Tstate, i0)
    solution = derive(Tsolution, grid, y, i0)
    x = pde(apm, state, solution)
    A = @NT()
    if length(x) == 3
        names = keys(x[3])
        A = _NT(names, Array(Float64, size(grid)))
        for i in eachindex(grid)
            state = getindex(grid, Tstate, i)
            solution = derive(Tsolution, grid, y, i)
            drifti = pde(apm, state, solution)[2]
            # upwind
            solution = derive(Tsolution, grid, y, i, drifti)
            othersi = pde(apm, state, solution)[3]
            for j in 1:length(names)
                A[names[j]][i] = othersi[j]
            end
        end
    end
    return A
end


#========================================================================================

Solve

========================================================================================#

all_symbol(names) = vcat((map(x -> Symbol(x...), with_replacement_combinations(names, k)) for k in 0:2)...)
all_symbol(sols, states) = vec([Symbol(a, s) for s in all_symbol(states), a in sols])
_NT(names::Vector{Symbol}) =  eval(Expr(:macrocall, Symbol("@NT"), (x for x in names)...))
_NT(names::Vector{Symbol}, element) =  eval(Expr(:macrocall, Symbol("@NT"), (Expr(:kw, names[i], element) for i in 1:length(names))...))

function pde_solve(apm, grid::NamedTuple, y0::NamedTuple; is_algebraic = nothing, kwargs...)
    Tstate = _NT(keys(grid))
    Tsolution = _NT(all_symbol(keys(y0), keys(grid)))
    stategrid = StateGrid(grid)
    if is_algebraic == nothing
        is_algebraic = _NT(keys(y0))((fill(false, size(x)) for x in y0)...)
    end
    y, distance = nl_solve((y, ydot) -> hjb!(apm, stategrid, Tstate, Tsolution, y, ydot), _concatenate(y0); is_algebraic = _concatenate(is_algebraic), kwargs...)
    a = create_dictionary(apm, stategrid, Tstate, Tsolution, y)
    y = _deconcatenate(typeof(y0), y)
    return merge(y, a), distance
end

function _concatenate(y)
    if length(y) == 1
        y[1]
    else
        cat(ndims(y[1]) + 1, values(y)...)
    end
end
function _deconcatenate(T, y)
    if nfields(T) == 1
        T(y)
    else
        N = ndims(y) - 1
        T(map(i -> y[(Colon() for k in 1: N)..., i], 1:nfields(T))...)
    end
end
