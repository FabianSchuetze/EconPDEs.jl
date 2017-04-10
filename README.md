[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)

# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

This package proposes a new, fast, and robust algorithm to solve PDEs associated with economic models. I discuss in details this algorithm [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf).

# Solving  PDEs
To `solve` a PDE, the user needs to define a type and three functions. I go through these definitions for the PDE associated to the Campbell Cochrane (1999) model.
1. A type that stores the parameters of the PDEs. For the case of Campbell Cochrane (1999),
	```julia
	type CampbellCochraneModel  <: EconPDEModel
	    # consumption process parameters
	    μ::Float64 
	    σ::Float64

	    # utility parameters
	    γ::Float64
	    ρ::Float64

	    # habit parameters
	    κs::Float64
	end
	# initialize
	function CampbellCochraneModel(;μ = 0.0189, σ = 0.015, γ = 2.0, ρ = 0.116, κs = 0.13)
	    CampbellCochraneModel(μ, σ, γ, ρ, κs)
	end
	```
2. a `state_grid` function that returns the state space on which the PDE must be solved. For the case of Campbell Cochrane (1999), there is only one state variable, the habit.
	```julia
	function state_grid(m::CampbellCochraneModel; smin = -300.0, n = 1000)
	    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs
	    Sbar = σ * sqrt(γ / κs)
	    sbar = log(Sbar)
	    smax = sbar + 0.5 * (1 - Sbar^2)
	    s = logspace(- 5, smax, n)
	    @NT(s = s)
	end
	```
3. an `initialize` function that returns an initial guess for the solution. For the case of Campbell Cochrane (1999) there is only one function to solve for, the price-dividend ratio. I simply take a vector of ones:

	```julia
	function initialize(m::CampbellCochraneModel, grid::StateGrid)
	    @NT(p = ones(grid.s))
	end
	```
4. a `pde` function that encodes the PDE. The function takes as argument the model `m`, the value of state variables `state`, and a current guess for the solution `solution`. It must return  a tuple of two terms.
	1. A tuple corresponding to the value of the system of PDEs at this grid point.
	2. A tuple corresponding to the drift of state variables at this grid point (used for upwinding).

	This is the `pde` function for the Campbell Cochrane (1999) model:
	```julia
	function pde(m::CampbellCochraneModel, state, solution)
	    μ = m.μ ; σ = m.σ ; γ = m.γ ; ρ = m.ρ ; κs = m.κs
	    s = state.s
	    p, ps, pss = solution.p, solution.ps, solution.pss
	    
	    # drift and volatility of state variable s
	    Sbar = σ * sqrt(γ / κs)
	    sbar = log(Sbar)
	    λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
	    μs = - κs * (s - sbar)
	    σs = λ * σ

	    # market price of risk κ
	    κ = γ * (σ + σs)

	    # risk free rate  r
	    r = ρ + γ * μ - γ * κs / 2

	    # drift and volatility of p
	    σp = ps / p * σs
	    μp = ps / p * μs + 0.5 * pss / p * σs^2

	    # PDE
	    out = p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
	    return out, μs
	end
	```


# Examples

I have coded a variety of well-known economic models in continuous time in  `src/models`.
- Asset pricing model with time varying habit (Campbell Cochrane (1999), Wachter (2005))
- Asset pricing model with long run risk (Bansal Yaron (2004), Bansal, Kiku, Yaron (2009))
- Asset pricing model with heterogeneous agents (Garleanu Panageas (2015), DiTella (2016))

```julia
# Habit Model
## Campbell Cochrane (1999)
m = CampbellCochrane()
grid = state_grid(m)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)

## Wachter (2005) calibration
m = CampbellCochraneModel(μ = 0.022, σ = 0.0086, γ = 2.0, ρ = 0.073, κs = 0.116, b = 0.011 * 4)
grid = state_grid(m)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)


# Long Run Risk Models
## Bansal Yaron (2004)
m = BansalYaronModel()
grid = state_grid(m)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)
## Bansal, Kiku, Yaron (2009) calibration
m = BansalYaronModel(μbar = 0.018, νc = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189, ρ = 0.0132, γ = 7.5, ψ = 1.5)
grid = state_grid(m)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)

# Heterogeneous Agent Models
## Garleanu Panageas (2015)
m = GarleanuPanageasModel()
grid = state_grid(m)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)
## DiTella (2016)
m = DiTellaModel()
grid = state_grid(m)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)

# Consumption - saving problem with idiosyncratic income risk
## Wang Wang Yang (2016)
ap = WangWangYangModel()
grid = state_grid(ap)
y0 = initialize(ap, grid)
result, distance = pde_solve(ap, grid, y0)
```



# Internal Functions
`pde_solve` internally calls `nl_solve` to solve the non linear system associated with finite difference schemes. You can also call this function directly.

 The solver `nl_solve` has the following syntax. Denote `F` the finite difference scheme corresponding to a PDE. The goal is to find `y` such that `F(y) = 0`.
 - The first argument is a function `F!(y, out)` which transforms `out = F(y)` in place.
 - The second argument is an array of arbitrary dimension for the initial guess for `y`
 - The option `is_algebraic` (defaults to an array of `false`) is an array indicating the eventual algebraic equations (typically market clearing conditions).

 Some options control the algorithm:
 - The option `Δ` (default to 1.0) specifies the initial time step. 
 - The option `inner_iterations` (default to `10`) specifies the number of inner Newton-Raphson iterations. 
 - The option `autodiff` (default to `true`) specifies that the Jacobian is evaluated using automatic differentiation.



