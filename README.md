[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)

# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

This package proposes a new, fast, and robust algorithm to solve economic models in continuous time. This package was developed while writing my PhD thesis.



More precisely, the package includes 
1. a fast and robust function `Ψtc` to solve systems of PDEs + algebraic equations
2. a higher-level function `solve` to solve economics models that uses `Ψtc` in the background. The goal is to solve a variety of well-known economic models in continuous time through a common framework. The function directly solves
	- Asset pricing model with *time varying habit* (Campbell Cochrane (1999), Wachter (2005))
	- Asset pricing model with *long run risk* (Bansal Yaron (2004), Bansal, Kiku, Yaron (2009))
	- Asset pricing model with *heterogeneous agents* (Garleanu Panageas (2015), DiTella (2016))

You can also use `solve` to solve a new model. You just need to write a few functions that specify the model. Examples can be found in the folder `src/models`.


# `Ψtc` solves systems of PDES
The function `Ψtc` allows to solve systems of PDEs + eventual algebraic equations.

 The solver `Ψtc` has the following syntax. 
 - The first argument is the model, given as a function `F!(y, out)`.
 - The second argument is an initial value for `y`
 - The option `is_algebraic` (defaults to an array of `false`) is an array indicating the eventual algebraic equations (typically market clearing conditions).

 Some options control the algorithm:
 - The option `Δ` (default to 1.0) specifies the initial time step. 
 - The option `inner_iterations` (default to `10`) specifies the number of inner Newton-Raphson iterations. 
 - The option `autodiff` (default to `true`) specifies that the Jacobian is evaluated using automatic differentiation.

 The algorithm underlying this PDE solver is explained in details [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf)


# `solve` solves  economic models

The function `solve` is a higher-level function to solve economic models. In the background, the function still relies on the PDE solver `Ψtc`. It is a higher-level function in the sense that it reduces the boilerplate needed to solve models (in particular, the model automatically computes the finite difference derivatives through upwinding).

To `solve` a economic model, the user only needs to define a type and three functions.
- A type that stores the parameters of the models. For the case of Campbell Cochrane (1999),
```julia
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
# initialize
function CampbellCochraneModel(;μ = 0.0189, σ = 0.015, γ = 2.0, ρ = 0.116, κs = 0.138, b = 0.0)
    # I choose persistence so that monthly simulation of the model matches processes in CC (1999)
    # ρ = 12 * (1 - 0.89^(1/12))
    # κs = 12 * (1 - 0.87^(1/12))
    CampbellCochraneModel(μ, σ, γ, ρ, κs, b)
end
```
- a `Stategrid` function that creates the state space grid. For the case of Campbell Cochrane (1999),
```julia
function StateGrid(m::CampbellCochraneModel; smin = -300.0, n = 1000)
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    smax =  sbar + 0.5 * (1 - Sbar^2)
    s = linspace(smin, smax, n)
    StateGrid(s = s)
end
```
- an `initialize` function that returns an initial guess. I simply takes a vector of ones for the case of Campbell Cochrane (1999):

```julia
function initialize(m::CampbellCochraneModel, grid::StateGrid)
    fill(1.0, size(grid)...)
end
```
- a `pde` function that returns the system of PDEs. More precisely, the function takes as argument a current guess and a grid position. It returns  a tuple of three terms.
	1. A tuple corresponding to the value of the PDEs at this grid point.
	2. A tuple corresponding to the drift of state variables at this grid point (used for upwinding).
	3. A dictionary from symbols to values. This dictionary simply stores side functions computed while writing the PDE.

For the case of Campbell Cochrane (1999),
```julia
function pde(m::CampbellCochraneModel, grid, y, ituple, idrift = (0.0, 0.0))
    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs ; b = m.b
    s, = grid[ituple]
    p, ps, pss  = derive(grid, y[1], ituple, idrift)
    
    # drift and volatility of state variable s
    Sbar = σ * sqrt(γ / (κs - b / γ))
    sbar = log(Sbar)
    λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
    μs = - κs * (s - sbar)
    σs = λ * σ

    # market price of risk κ
    κ = γ * (σ + σs)

    # risk free rate  r
    r = ρ + γ * μ - (γ * κs - b) / 2 + b * (sbar - s)

    # drift and volatility of p
    σp = ps / p * σs
    μp = ps / p * μs + 0.5 * pss / p * σs^2

    # PDE
    out = p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
    return out, (μs,), (:p => p, :κ => κ, :λ => λ, :r => r, :σp => σp, :μs => μs, :σs => σs)
end
```

Each model is coded as a system of PDEs, where each PDE corresponds to the no-arbitrage condition for an asset. All the models can be found in `src/models`. 

```julia
using EconPDEs 

# Habit Models
## Campbell Cochrane (1999)
m = CampbellCochraneModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)
using Plots
plotly()
plot(exp(grid[:s]), result[:p])
## Wachter (2005) calibration:
m = CampbellCochraneModel(μ = 0.022, σ = 0.0086, γ = 2.0, ρ = 0.073, κs = 0.116, b = 0.011 * 4)
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)

# Long Run Risk Models
## Bansal Yaron (2004)
m = BansalYaronModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)
using Plots
plotly()
surface(grid[:μ], grid[:σ], result[:p])
## Bansal, Kiku, Yaron (2009) calibration
m = BansalYaronModel(μbar = 0.018, νc = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189, ρ = 0.0132, γ = 7.5, ψ = 1.5)
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)

# Heterogeneous Agent Models
## Garleanu Panageas (2015)
m = GarleanuPanageasModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)
using Plots
plotly()
plot(grid[:x], result[:p])
## DiTella (2016)
using EconPDEs
m = DiTellaModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
is_algebraic = fill(false, size(y0)...)
is_algebraic[:, :, 3] = true
result, distance = solve(m, grid, y0, is_algebraic = is_algebraic)
using Plots
plotly()
surface(grid[:x], grid[:ν], result[:p])

# Consumption - saving problem with idiosyncratic income risk
## Wang Wang Yang (2016)
ap = WangWangYangModel()
grid = StateGrid(ap)
y0 = initialize(ap, grid)
result, distance = solve(ap, grid, y0)

using Plots
plotly()
plot(grid[:w], result[:p])
```

