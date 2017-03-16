[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)

# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

Economic models can be hard to solve numerically. This package implements a fast and robust algorithm to solve them. I demonstrate the interest of this algorithm by solving a variety of well-known economic models in continuous time. This package was developed while writing my PhD thesis.



More precisely, the package includes 
1. a fast and robust function `Ψtc` to solve systems of PDEs + algebraic equations
2. a higher-level function `solve` to solve  economics models. The function directly solves
	- Asset pricing model with *time varying habit* (Campbell Cochrane (1999), Wachter (2005))
	- Asset pricing model with *long run risk* (Bansal Yaron (2004), Bansal, Kiku, Yaron (2009))
	- Asset pricing model with *heterogeneous agents* (Garleanu Panageas (2015), DiTella (2016))


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

The function `solve` is a higher-level function to solve economic models. In the background, the function still relies on the PDE solver `Ψtc` in the background. The goal of this function is to reduce the boilerplate when solving new models.

To `solve` a economic model, the user only needs to define three functions.
- a `Stategrid` function that creates the state space grid specific to the economic problem.
- an `initialize` function that returns an initial guess for the functions that solve the system of PDEs.
- a `pde` function that returns the system of PDEs. More precisely, the function takes as argument a current guess and a grid position. It returns  a tuple of three terms.
	1. The first term is a tuple corresponding to the system of equations evaluated at this grid point.
	2. the second term is a tuple for the drift of state variables (used to upwind the scheme).
	3. the third term is a dictionary from symbols to value. It is not used to solve the model, but it allows to compute interesting quantities beyond the solution of the PDEs.

Examples of these functions can be found in the folder `src/models`. 	Each model is coded as a system of PDEs, where each PDE corresponds to the no-arbitrage condition for an asset.

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

