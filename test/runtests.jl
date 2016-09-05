using EconPDEs, Base.Test


m = CampbellCochraneModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)
@test distance <= 1e-5


m = BansalYaronModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)
@test distance <= 1e-5


m = BansalYaronModel(μbar = 0.018, νD = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189, ρ = 0.0132, γ = 7.5, ψ = 1.5)
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)
@test distance <= 1e-5


m = GarleanuPanageasModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
result, distance = fullsolve(m, grid, y0)
@test distance <= 1e-5


m = DiTellaModel()
grid = StateGrid(m)
y0 = initialize(m, grid)
is_algebraic = fill(false, size(y0)...)
is_algebraic[:, :, 3] = true
result, distance = fullsolve(m, grid, y0, is_algebraic = is_algebraic)
@test distance <= 1e-5


ap = WangWangYangModel()
grid = StateGrid(ap)
y0 = initialize(ap, grid)
result, distance = fullsolve(ap, grid, y0)
@test distance <= 1e-5

