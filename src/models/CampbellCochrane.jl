##############################################################################
##
## Type
##
##############################################################################

type CampbellCochraneModel  <: EconPDEModel
    # consumption process parameters
    μ::Float64 
    σ::Float64

    # utility
    γ::Float64
    ρ::Float64

    # habit
    κs::Float64
    b::Float64
end

function CampbellCochraneModel(;μ = 0.0189, σ = 0.015, γ = 2.0, ρ = 0.116, κs = 0.138, b = 0.0)
    # I choose persistence so that monthly simulation of the model matches processes in CC (1999)
    # ρ = 12 * (1 - 0.89^(1/12))
    # κs = 12 * (1 - 0.87^(1/12))
    CampbellCochraneModel(μ, σ, γ, ρ, κs, b)
end


function StateGrid(m::CampbellCochraneModel; smin = -300.0, n = 1000)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    smax =  sbar + 0.5 * (1 - Sbar^2)
    # corresponds to Grid 3 in Wachter (2005)
    shigh = log(linspace(0.0, exp(smax), div(n, 10)))
    slow = linspace(smin, shigh[2], n - div(n, 10))
    s = vcat(slow[1:(end-1)], shigh[2:end])
    StateGrid(s = s)
end

function initialize(m::CampbellCochraneModel, grid::StateGrid)
    fill(1.0, size(grid)...)
end
	
function derive(m::CampbellCochraneModel, stategrid::StateGrid, y::ReflectingArray, ituple, drifti = 0.0)
    is = ituple[1]
    μs = drifti
    Δs, = stategrid.Δx
    Δsm, = stategrid.Δxm
    Δsp, = stategrid.Δxp
    p = y[is]
    if μs >= 0.0
        ps = (y[is + 1] - y[is]) / Δsp[is]
    else
        ps = (y[is] - y[is - 1]) / Δsm[is]
    end
    pss = (Δsm[is] * y[is + 1] + Δsp[is] * y[is - 1] - 2 * Δs[is] * y[is]) / (Δs[is] * Δsm[is] * Δsp[is])
    return p, ps, pss
end

function pde(m::CampbellCochraneModel, gridi, functionsi)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    s, = gridi
    p, ps, pss = functionsi[1]
    # evolution state variable
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    μs = - κs * (s - sbar)
    λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
    σs = λ * σ
    # sdf
    r = ρ + γ * μ - (γ * κs - b) / 2 + b * (sbar - s)
    κ = γ * (σ + σs)
    # wealth / consumption
    σp = ps / p * σs
    μp = ps / p * μs + 0.5 * pss / p * σs^2
    out = p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
    return out, μs, (:p => p, :κ => κ, :λ => λ, :r => r, :σp => σp, :μs => μs, :σs => σs)
end

