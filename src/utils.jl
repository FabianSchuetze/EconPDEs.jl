
#========================================================================================

Stationary Distribution

========================================================================================#

# Case with 1 state variable
function stationary_distribution(grid::StateGrid{1}, a)
    n, = size(grid)
    invΔx, = grid.invΔx
    invΔxp, = grid.invΔxp
    invΔxm, = grid.invΔxm
    A = ReflectingArray(zeros(n, n))
    for i in 1:n
        μ = a[Symbol(:μ, grid.name[1])][i]
        σ2 = a[Symbol(:σ, grid.name[1])][i]^2
        if μ >= 0
            A[i + 1, i] += μ * invΔxp[i]
            A[i, i] -= μ * invΔxp[i]
        else
            A[i, i] += μ * invΔxm[i] 
            A[i - 1, i] -= μ * invΔxm[i] 
        end
        A[i - 1, i] += 0.5 * σ2 * invΔx[i] * invΔxm[i] 
        A[i, i] -= 0.5 * σ2 * 2 * invΔxm[i] * invΔxp[i]
        A[i + 1, i] += 0.5 * σ2 * invΔx[i] * invΔxp[i]
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

function stationary_distribution(grid::StateGrid{2}, a)
    n1, n2 = size(grid)
    A = ReflectingArray(zeros(n1, n2, n1, n2))
    invΔx1, invΔx2 = grid.invΔx
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
            A[i1h, i2, i1, i2] += μ * invΔx1[i1]
            A[i1l, i2, i1, i2] -= μ * invΔx1[i1]
            A[i1 - 1, i2, i1, i2] += 0.5 * σ2 * invΔx1[i1]^2
            A[i1, i2, i1, i2] -= 0.5 * σ2 * 2 * invΔx1[i1]^2
            A[i1 + 1,i2, i1, i2] += 0.5 * σ2 * invΔx1[i1]^2

            μ = a[Symbol(:μ, grid.name[2])][i1, i2]
            σ2 = a[Symbol(:σ, grid.name[2], :2)][i1, i2]
            if μ >= 0
                i2h = i2 + 1
                i2l = i2
            else
               i2h = i2
               i2l = i2 - 1
            end
            A[i1, i2h, i1, i2] += μ * invΔx1[i1]
            A[i1, i2l, i1, i2] -= μ * invΔx1[i1]
            A[i1, i2 - 1, i1, i2] += 0.5 * σ2 * invΔx2[i2]^2
            A[i1, i2, i1, i2] -= 0.5 * σ2 * 2 * invΔx2[i2]^2
            A[i1,i2 + 1, i1, i2] += 0.5 * σ2 * invΔx2[i2]^2

            σ12 = a[Symbol(:σ, grid.name[1], :σ, grid.name[2])][i1, i2]
            A[i1h, i2h, i1, i2] += σ12 * invΔx1[i1] * invΔx2[i2]
            A[i1l, i2h, i1, i2] -= σ12 * invΔx1[i1] * invΔx2[i2]
            A[i1h, i2l, i1, i2] -= σ12 * invΔx1[i1] * invΔx2[i2]
            A[i1l, i2l, i1, i2] += σ12 * invΔx1[i1] * invΔx2[i2]
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

function simulate(grid, a, shocks::Dict; dt = 1 / 12, x0 = nothing)
    grid = StateGrid(grid)
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
