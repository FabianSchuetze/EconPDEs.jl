type GarleanuPanageasModel <: EconPDEModel # Where is EconPDEModel coming from? Looks like something like object inheritance

  # utility function
  γA::Float64
  ψA::Float64
  γB::Float64
  ψB::Float64
  ρ::Float64
  δ::Float64

  # proportion a
  νA::Float64

  # consumption
  μ::Float64
  σ::Float64

  # earning function
  B1::Float64
  δ1::Float64
  B2::Float64
  δ2::Float64
  ω::Float64
end

function GarleanuPanageasModel(;γA  = 1.5, ψA = 0.7, γB = 10.0, ψB = 0.05, ρ = 0.001, δ = 0.02, νA = 0.01, μ = 0.02, σ = 0.041, B1 = 30.72, δ1 = 0.0525, B2 = -30.29, δ2 = 0.0611, ω = 0.92)
  """
  Parameters:
  -----------

  \rho: scalar(float):
    The subjective discount factor (surprisingly small! How can it become so small? )

  \delta: scalar(float):
    The hazard-rate of death. How does this affect the consumption-wealth ratio?

  \nuA: scalar(floar):
    Changes the mean return of each group
  """
  scale = δ / (δ + δ1) * B1 + δ / (δ + δ2) * B2  # what is this scale variable?
  B1 = B1 / scale
  B2 = B2 / scale
  GarleanuPanageasModel(γA , ψA, γB, ψB, ρ, δ, νA, μ, σ, B1, δ1, B2, δ2, ω)
end

function state_grid(m::GarleanuPanageasModel; n = 200) # It looks like this function accecpts all objectes of the type GarleanuPanageasModel!
  @NT(x = linspace(0.0, 1.0, n)) # Why does it makes sense to use a tuple here?
end

function initialize(m::GarleanuPanageasModel, grid)
    x = fill(1.0, length(grid.x))
    @NT(pA = x, pB = x, ϕ1 = x, ϕ2 = x) # what is ϕ1 and ϕ2?
end

@inline function pde(m::GarleanuPanageasModel, state, solution) #what is the solution object?
  γA = m.γA ; ψA = m.ψA ; γB = m.γB ; ψB = m.ψB ; ρ = m.ρ ; δ = m.δ ; νA = m.νA ; μ = m.μ ; σ = m.σ; B1 = m.B1 ; δ1 = m.δ1 ; B2 = m.B2 ; δ2 = m.δ2 ; ω = m.ω
  x = state.x
  pA, pAx, pAxx, pB, pBx, pBxx, ϕ1, ϕ1x, ϕ1xx, ϕ2, ϕ2x, ϕ2xx = solution.pA, solution.pAx, solution.pAxx, solution.pB, solution.pBx, solution.pBxx, solution.ϕ1, solution.ϕ1x, solution.ϕ1xx, solution.ϕ2, solution.ϕ2x, solution.ϕ2xx

  # volatility of X, pA, pB, ϕ1, ϕ2, CA, CB and market price of risk κ
  Γ = 1 / (x / γA + (1 - x) / γB) ## aggregate risk aversion
  p = x * pA + (1 - x) * pB ## aggregate wealth-consumption ratio
  σx = σ * x * (Γ / γA - 1) / (1 + Γ * x * (1 - x) / (γA * γB) * ((1 - γB * ψB) / (ψB - 1) * (pBx / pB) - (1 - γA * ψA) / (ψA - 1) * (pAx / pA)))
  σpA = pAx / pA * σx
  σpB = pBx / pB * σx
  σϕ1 = ϕ1x / ϕ1 * σx
  σϕ2 = ϕ2x / ϕ2 * σx
  κ = Γ * (σ - x * (1 - γA * ψA) / (γA * (ψA - 1)) * σpA - (1 - x) * (1 - γB * ψB) / (γB * (ψB - 1)) * σpB)
  σCA = κ / γA + (1 - γA * ψA) / (γA * (ψA - 1)) * σpA
  σCB = κ / γB + (1 - γB * ψB) / (γB * (ψB - 1)) * σpB

  # drift of X, pA, pB, ϕ1, ϕ2, CA, CB and interest rate r
  mcA = κ^2 * (1 + ψA) / (2 * γA) + (1 - ψA * γA) / (γA * (ψA - 1)) * κ * σpA - (1 - γA * ψA) / (2 * γA * (ψA - 1)) * σpA^2
  mcB = κ^2 * (1 + ψB) / (2 * γB) + (1 - ψB * γB) / (γB * (ψB - 1)) * κ * σpB - (1 - γB * ψB) / (2 * γB * (ψB - 1)) * σpB^2
  r =  ρ + 1 / (ψA * x  + ψB * (1 - x))  * (μ - x * mcA - (1 - x) * mcB - δ * ((νA / pA + (1 - νA) / pB) * (ϕ1 + ϕ2) - 1))
  μCA = ψA * (r - ρ) + mcA
  μCB = ψB * (r - ρ) + mcB
  μx = x * (μCA - δ - μ) + δ * νA / pA * (ϕ1 + ϕ2) - σ * σx
  μpA = pAx / pA * μx + 0.5 * pAxx / pA * σx^2
  μpB = pBx / pB * μx + 0.5 * pBxx / pB * σx^2
  μϕ1 = ϕ1x / ϕ1 * μx + 0.5 * ϕ1xx / ϕ1 * σx^2
  μϕ2 = ϕ2x / ϕ2 * μx + 0.5 * ϕ2xx / ϕ2 * σx^2

  # PDE
  ## == is this A tuple corresponding to the value of the system of PDEs at this grid point ?? == ##
  out1 = pA * (1 / pA + (μCA - δ) + μpA + σCA * σpA - r - κ * (σpA + σCA)) # this looks like the difference. The aim should be out = 0?
  out2 = pB * (1 / pB + (μCB - δ) + μpB + σCB * σpB - r - κ * (σpB + σCB))
  out3 = ϕ1 * (B1 * ω / ϕ1 + (μ - δ - δ1) + μϕ1 + σ * σϕ1 - r - κ * (σϕ1 + σ))
  out4 = ϕ2 * (B2 * ω / ϕ2 + (μ - δ - δ2) + μϕ2 + σ * σϕ2 - r - κ * (σϕ2 + σ))

  return (out1, out2, out3, out4), μx, @NT(p = p, pA = pA, pB = pB, κ = κ, r = r, μx = μx, σx = σx)
end
