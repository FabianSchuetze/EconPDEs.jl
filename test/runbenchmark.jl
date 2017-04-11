using EconPDEs, Base.Test


m = CampbellCochraneModel()
grid = state_grid(m)
y0 = initialize(m, grid)
@time pde_solve(m, grid, y0)
@time pde_solve(m, grid, y0)
#0.45s, 93MB

m = BansalYaronModel()
grid = state_grid(m)
y0 = initialize(m, grid)
@time pde_solve(m, grid, y0)
@time pde_solve(m, grid, y0)
#1.13s, 150MB


m = GarleanuPanageasModel()
grid = state_grid(m)
y0 = initialize(m, grid)
@time pde_solve(m, grid, y0)
@time pde_solve(m, grid, y0)
#1.9s, 150MB

m = DiTellaModel()
grid = state_grid(m)
y0 = initialize(m, grid)
@time pde_solve(m, grid, y0)
@time pde_solve(m, grid, y0)
# 30.820142s, 10.2 GB
@time pde_solve(m, grid, y0, is_algebraic = (false, false, true))

m = WangWangYangModel()
grid = state_grid(m)
y0 = initialize(m, grid)
@time pde_solve(m, grid, y0)
@time pde_solve(m, grid, y0)
# 0.08s, 3.0MB
