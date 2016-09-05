##############################################################################
##
## Type
##
##############################################################################
type WangWangYangModel <: EconPDEModel
    μ::Float64 
    σ::Float64
    r::Float64
    ρ::Float64  
    γ::Float64 
    ψ::Float64
end

function WangWangYangModel(;μ = 0.01, σ = 0.1, r = 0.05, ρ = 0.055, γ = 4, ψ = 0.5)
    WangWangYangModel(μ, σ, r, ρ, γ, ψ)
end

function StateGrid(m::WangWangYangModel; wn = 100)
    StateGrid(w = collect(linspace(0.0, 30.0, wn)))
end

function initialize(m::WangWangYangModel, grid::StateGrid)
    grid.x[1]
end
	
function derive(m::WangWangYangModel, stategrid::StateGrid, y::ReflectingArray, ituple, drifti = 0.0)
    μ = m.μ ;  σ = m.σ ;  r = m.r ;  ρ = m.ρ ;  γ = m.γ ;  ψ = m.ψ 
    iw = ituple[1]
    μw = drifti
    Δw, = stategrid.Δx
    p = max(1e-10, y[iw])
    n = length(Δw)
    # At boundaries, take derivatives forward / backward
    if μw >= 0.0
        pw = (y[iw + 1] - y[iw]) / Δw[iw]
    else
        pw = (y[iw] - y[iw-1]) / Δw[iw]
    end
    pw = max(1e-10, pw)
    if iw == 1
        pww = (y[iw + 2] + y[iw] - 2 * y[iw + 1]) / Δw[iw + 1]^2
    elseif iw == n
        pww = (y[iw] + y[iw - 2] - 2 * y[iw - 1]) / Δw[iw - 1]^2
    else
        pww = (y[iw + 1] + y[iw - 1] - 2 * y[iw]) / Δw[iw]^2
    end
    # financial friction: check consumption < 1 when w = 0
    if iw == 1 
        m = r + ψ * (ρ - r)
        c = m * p * pw^(-ψ)
        if c >= 1.0
            pw = (m * p)^(1 / ψ)
        end
    end
    return p, pw, pww
end

function pde(m::WangWangYangModel, gridi, functionsi)
    μ = m.μ ;  σ = m.σ ;  r = m.r ;  ρ = m.ρ ;  γ = m.γ ;  ψ = m.ψ 
    w, = gridi
    p, pw, pww = functionsi
    m = r + ψ * (ρ - r)
    c = m * p * pw^(-ψ)
    out = ((m * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = (r - μ + σ^2) * w + 1 - c
    return out, μw, (:w => w, :p => p, :pw => pw, :pww => pww, :μw => μw, :c => c)
end
