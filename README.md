[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)

# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

The package includes (i) a fast and robust PDE solver (ii) a higher lever solver for economics models in continuous time. 

The PDE solver is based on implicit time-stepping, as explained in details [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf). The function handles systems of PDEs + eventual algebraic equations.




# PDE Solver

 The solver `Ψtc` has the following syntax. 
 - The first argument is the model, given as a function `F!(y, out)`.
 - The second argument is an initial value for `y`
 - The option `is_algebraic` (defaults to an array of `false`) is an array indicating the eventual algebraic equations (typically market clearing conditions).

 Some options control the algorithm:
 - The option `Δ` (default to 1.0) specifies the initial time step. 
 - The option `inner_iterations` (default to `10`) specifies the number of inner Newton-Raphson iterations. 
 - The option `autodiff` (default to `true`) specifies that the Jacobian is evaluated using automatic differentiation.


# Higher-level Solver for Economics Models
I use the PDE solver `Ψtc` to solve several models in continous time.

Each model only needs to define four elements
- a type that includes the model parameters
- a `Stategrid` function that creates the state space grid
- an `initialize` function that corresponds to the initial solution of the problem
- a `pde` function that codes the PDE of the model.

Models are coded in the folder `src/models`.

### Campbell Cochrane (1999)
Asset pricing model with time varying habit
```julia
using EconPDEs
m = CampbellCochraneModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)

# plot results
using Plots
plotly()
plot(exp(grid[:s]), result[:p])
```

Wachter (2005) calibration:
```julia
m = CampbellCochraneModel(μ = 0.022, σ = 0.0086, γ = 2.0, ρ = 0.073, κs = 0.116, b = 0.011 * 4)
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)
```



### Bansal Yaron (2004)
Asset pricing model with time varying drift and volatility

```julia
using EconPDEs
m = BansalYaronModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)

# plot results
using Plots
plotly()
surface(grid[:μ], grid[:σ], result[:p])
```

Bansal, Kiku, Yaron (2009) calibration:
```julia
using EconPDEs
m = BansalYaronModel(μbar = 0.018, νc = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189, ρ = 0.0132, γ = 7.5, ψ = 1.5)
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)
```

### Garleanu Panageas (2015)
Asset pricing model with heterogeneous agents
```julia
using EconPDEs
m = GarleanuPanageasModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)

# plot results
using Plots
plotly()
plot(grid[:x], result[:p])
```

### DiTella (2016)
Asset pricing model with heterogeneous agents and time varying volatility

```julia
using EconPDEs
m = DiTellaModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
is_algebraic = fill(false, size(y0)...)
is_algebraic[:, :, 3] = true
result, distance = fullsolve(m, grid, y0, is_algebraic = is_algebraic)

# plot results
using Plots
plotly()
surface(grid[:x], grid[:ν], result[:p])
```


# Macro Models

### Wang Wang Yang (2016)
Consumption - saving problem with idiosyncratic income risk
```julia
using EconPDEs
ap = WangWangYangModel()
grid = StateGrid(ap)
y0 = initialize(ap, grid)
result, distance = fullsolve(ap, grid, y0)

# plot results
using Plots
plotly()
plot(grid[:w], result[:p])
```

