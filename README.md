# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

The package introduces a novel, fast and robust solver for economics models in continuous time. 

The package can solves PDEs, systems of PDEs, or systems of PDEs + algebraic equations.

The underlying algorithm is explained in details [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf).



# Syntax

 The solver `Ψtc` has the following syntax. 
 - The first argument is the model, given as a function `F!(y, out)`.
 - The second argument is an initial value for `y`
 - The option `is_algebraic` (defaults to an array of `false`) is an array indicating the eventual algebraic equations (typically market clearing conditions).

 Some options control the algorithm:
 - The option `Δ` (default to 1.0) specifies the initial time step. 
 - The option `inner_iterations` (default to `10`) specifies the number of inner Newton-Raphson iterations. 
 - The option `autodiff` (default to `true`) specifies that the Jacobian is evaluated using automatic differentiation.


# Models
The packages solves several asset pricing models to show how `Ψct` works. Models are coded in the folder `src/models`. I set the initial condition of each model to the constant `1.0` to show the robustness of the solution method.

Open an issue or send me an email if you spot mistakes. Pull requests are welcome.

### Bansal Yaron (2004)
```julia
using EconPDEs
ap = BansalYaronModel()
grid = StateGrid(ap)
y0 = initialize(ap, grid)
result, distance = fullsolve(ap, grid, y0)

# plot results
using Plots
plotly()
surface(grid[:μ], grid[:σ], result[:p])
```

### Garleanu Panageas (2015)

```julia
using EconPDEs
ap = GarleanuPanageasModel()
grid = StateGrid(ap)
y0 = initialize(ap, grid)
result, distance = fullsolve(ap, grid, y0)

# plot results
using Plots
plotly()
plot(grid[:x], result[:p])
```

### DiTella (2016)

```julia
using EconPDEs
ap = DiTellaModel()
grid = StateGrid(ap)
y0 = initialize(ap, grid)
is_algebraic = fill(false, size(y0)...)
is_algebraic[:, :, 3] = true
result, distance = fullsolve(ap, grid, y0, is_algebraic = is_algebraic)

# plot results
using Plots
plotly()
surface(grid[:x], grid[:ν], result[:p])
```
