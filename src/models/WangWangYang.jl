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

function state_grid(m::WangWangYangModel; n = 100)
    @NT(w = collect(linspace(0.0, 30.0, n)))
end

function initialize(m::WangWangYangModel, grid)
    @NT(p = grid.w)
end
	
function pde(m::WangWangYangModel, state, solution)
    μ = m.μ ;  σ = m.σ ;  r = m.r ;  ρ = m.ρ ;  γ = m.γ ;  ψ = m.ψ 
    w = state.w
    p, pw, pww = solution.p, solution.pw, solution.pww
    p = max(1e-10, p)
    pw = max(1e-10, pw)
    # financial friction: check consumption < 1 when w = 0
    #if ituple[1] == 1 
    #    m = r + ψ * (ρ - r)
    #    c = m * p * pw^(-ψ)
    #    if c >= 1.0
    #        pw = (m * p)^(1 / ψ)
    #    end
    #end
    m = r + ψ * (ρ - r)
    c = m * p * pw^(-ψ)
    out = ((m * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = (r - μ + σ^2) * w + 1 - c
    return out, (μw,), @NT(w = w, p = p, pw = pw, pww = pww, μw = μw, c = c)
end
