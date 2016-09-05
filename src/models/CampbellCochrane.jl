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


function StateGrid(m::CampbellCochraneModel; smin = -100.0, sn = 100)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    smax =  sbar + 0.5 * (1 - Sbar^2)
    s1 = linspace(smin, log(Sbar / 100), 100)
    s2 = log(linspace(Sbar / 100, exp(smax), 100))
    s = vcat(s1, s2[2:end])
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
    pss = (y[is + 1] + y[is - 1] - 2 * y[is]) / Δs[is]^2
    return p, ps, pss
end

function pde(m::CampbellCochraneModel, gridi, functionsi)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    s, = gridi
    p, ps, pss = functionsi
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    μs = - κs * (s - sbar)
    λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
    σs = λ * σ
    r = ρ + γ * μ - (γ * κs - b) / 2 + b * (sbar - s)
    κ = γ * (σ + σs)
    σp = ps / p * σs
    μp = ps / p * μs + 0.5 * pss / p * σs^2
    out = p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
    return out, μs, (:p => p, :κ => κ, :λ => λ, :r => r, :σp => σp)
end

