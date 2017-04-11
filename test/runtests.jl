using EconPDEs, Base.Test


m = CampbellCochraneModel()
grid = state_grid(m; n = 1000)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)
@test distance <= 1e-5


m = BansalYaronModel()
grid = state_grid(m; μn = 5, σn = 5)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)
@test distance <= 1e-5



m = GarleanuPanageasModel()
grid = state_grid(m; n = 10)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)
@test distance <= 1e-5


m = DiTellaModel()
grid = state_grid(m ; xn = 10, νn = 3)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)
@time pde_solve(m, grid, y0, is_algebraic = (false, false, true))
@test distance <= 1e-5


m = WangWangYangModel()
grid = state_grid(m; n = 10)
y0 = initialize(m, grid)
result, distance = pde_solve(m, grid, y0)
@test distance <= 1e-5


