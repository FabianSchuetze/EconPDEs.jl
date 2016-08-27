type GarleanuPanageasModel <: EconPDEModel

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
  scale = δ / (δ + δ1) * B1 + δ / (δ + δ2) * B2
  B1 = B1 / scale
  B2 = B2 / scale
  GarleanuPanageasModel(γA , ψA, γB, ψB, ρ, δ, νA, μ, σ, B1, δ1, B2, δ2, ω)
end

function StateGrid(apm::GarleanuPanageasModel; n = 200)
  StateGrid(x = linspace(0.0, 1.0, n))
end

function initialize(apm::GarleanuPanageasModel, grid::StateGrid)
    fill(1.0, size(grid)..., 4)
end

function derive(apm::GarleanuPanageasModel, stategrid::StateGrid, y::ReflectingArray, ituple, drifti = (0.0,))
  ix = ituple[1]
  μX, = drifti
  Δx, = stategrid.Δx

  pA = y[ix, 1]
  pB = y[ix, 2]
  ϕ1 = y[ix, 3]
  ϕ2 = y[ix, 4]
  if μX <= 0.0
    pAx = (y[ix, 1] - y[ix - 1, 1]) / Δx[ix]
    pBx = (y[ix, 2] - y[ix - 1, 2]) / Δx[ix]
    ϕ1x = (y[ix, 3] - y[ix - 1, 3]) / Δx[ix]
    ϕ2x = (y[ix, 4] - y[ix - 1, 4]) / Δx[ix] 
  else
    pAx = (y[ix + 1, 1] - y[ix, 1]) / Δx[ix]
    pBx = (y[ix + 1, 2] - y[ix, 2]) / Δx[ix]
    ϕ1x = (y[ix + 1, 3] - y[ix, 3]) / Δx[ix]
    ϕ2x = (y[ix + 1, 4] - y[ix, 4]) / Δx[ix]
  end
  pAxx = (y[ix + 1, 1] + y[ix - 1, 1] - 2 * y[ix, 1]) / Δx[ix]^2
  pBxx = (y[ix + 1, 2] + y[ix - 1, 2] - 2 * y[ix, 2]) / Δx[ix]^2
  ϕ1xx = (y[ix + 1, 3] + y[ix - 1, 3] - 2 * y[ix, 3]) / Δx[ix]^2
  ϕ2xx = (y[ix + 1, 4] + y[ix - 1, 4] - 2 * y[ix, 4]) / Δx[ix]^2
  return pA, pAx, pAxx, pB, pBx, pBxx, ϕ1, ϕ1x, ϕ1xx, ϕ2, ϕ2x, ϕ2xx
end

function pde(apm::GarleanuPanageasModel, gridi, functionsi)
  x, = gridi
  pA, pAx, pAxx, pB, pBx, pBxx, ϕ1, ϕ1x, ϕ1xx, ϕ2, ϕ2x, ϕ2xx = functionsi
  γA = apm.γA ; ψA = apm.ψA ; γB = apm.γB ; ψB = apm.ψB ; ρ = apm.ρ ; δ = apm.δ ; νA = apm.νA ; μ = apm.μ ; σ = apm.σ; B1 = apm.B1 ; δ1 = apm.δ1 ; B2 = apm.B2 ; δ2 = apm.δ2 ; ω = apm.ω ; 
  Γ = 1 / (x / γA + (1 - x) / γB)
  σX = σ * x * (Γ / γA - 1) / (1 + Γ * x * (1 - x) / (γA * γB) * ((1 - γB * ψB) / (ψB - 1) * (pBx / pB) - (1 - γA * ψA) / (ψA - 1) * (pAx / pA)))
  σpA = pAx / pA * σX
  σpB = pBx / pB * σX 
  σϕ1 = ϕ1x / ϕ1 * σX
  σϕ2 = ϕ2x / ϕ2 * σX
  κ = Γ * (σ - x * (1 - γA * ψA) / (γA * (ψA - 1)) * σpA - (1 - x) * (1 - γB * ψB) / (γB * (ψB - 1)) * σpB)
  σCA = κ / γA + (1 - γA * ψA) / (γA * (ψA - 1)) * σpA
  σCB = κ / γB + (1 - γB * ψB) / (γB * (ψB - 1)) * σpB
  mcA = κ^2 * (1 + ψA) / (2 * γA) + (1 - ψA * γA) / (γA * (ψA - 1)) * κ * σpA - (1 - γA * ψA) / (2 * γA * (ψA - 1)) * σpA^2
  mcB = κ^2 * (1 + ψB) / (2 * γB) + (1 - ψB * γB) / (γB * (ψB - 1)) * κ * σpB - (1 - γB * ψB) / (2 * γB * (ψB - 1)) * σpB^2
  r =  ρ + 1 / (ψA * x  + ψB * (1 - x))  * (μ - x * mcA - (1 - x) * mcB - δ * ((νA / pA + (1 - νA) / pB) * (ϕ1 + ϕ2) - 1))
  μCA = ψA * (r - ρ) + mcA
  μCB = ψB * (r - ρ) + mcB
  μX = x * (μCA - δ - μ) + δ * νA / pA * (ϕ1 + ϕ2) - σ * σX  
  μpA = pAx / pA * μX + 0.5 * pAxx / pA * σX^2
  μpB = pBx / pB * μX + 0.5 * pBxx / pB * σX^2
  μϕ1 = ϕ1x / ϕ1 * μX + 0.5 * ϕ1xx / ϕ1 * σX^2
  μϕ2 = ϕ2x / ϕ2 * μX + 0.5 * ϕ2xx / ϕ2 * σX^2

  out1 = pA * (1 / pA + μCA + μpA + σCA * σpA - r - δ - κ * (σpA + σCA))
  out2 = pB * (1 / pB + μCB + μpB + σCB * σpB - r - δ - κ * (σpB + σCB))
  out3 = ϕ1 * (B1 * ω / ϕ1 + (μ - δ - δ1) + μϕ1 + σ * σϕ1 - r - κ * (σϕ1 + σ))
  out4 = ϕ2 * (B2 * ω / ϕ2 + (μ - δ - δ2) + μϕ2 + σ * σϕ2 - r - κ * (σϕ2 + σ))

  p = x * pA + (1 - x) * pB
  return (out1, out2, out3, out4), (μX,), (:p => p, :pA => pA, :pB => pB, :κ => κ, :r => r)
end


