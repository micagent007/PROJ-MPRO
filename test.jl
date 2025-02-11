include("coupes.jl")
include("dual.jl")
include("problem_static.jl")
include("Heuristic.jl")

path = "data/n_8-euclidean_true"

Cutting_planes(path)
#Dual_solve(path)
#Static_problem(path)
Heuristic(path)