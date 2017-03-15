using EconPDEs, Base.Test


m = CampbellCochraneModel()
grid = StateGrid(m; n = 1000)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)
EconPDEs.simulate(grid, result, Dict(:Z => randn(10)))
@test distance <= 1e-5


m = BansalYaronModel()
grid = StateGrid(m; μn = 5, σn = 5)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)
EconPDEs.simulate(grid, result, Dict(:Zμ => randn(10), :Zσ => randn(10)))
@test distance <= 1e-5



m = GarleanuPanageasModel()
grid = StateGrid(m; n = 10)
y0 = initialize(m, grid)
result, distance = solve(m, grid, y0)
@test distance <= 1e-5


m = DiTellaModel()
grid = StateGrid(m ; xn = 10, νn = 3)
y0 = initialize(m, grid)
is_algebraic = fill(false, size(y0)...)
is_algebraic[:, :, 3] = true
result, distance = solve(m, grid, y0, is_algebraic = is_algebraic)
@test distance <= 1e-5


ap = WangWangYangModel()
grid = StateGrid(ap; n = 10)
y0 = initialize(ap, grid)
result, distance = solve(ap, grid, y0)
@test distance <= 1e-5


