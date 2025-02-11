include("coupes.jl")
include("dual.jl")
include("problem_static.jl")
include("Heuristic.jl")

path = "data/n_7-euclidean_false"

Cutting_planes(path)
Dual_solve(path)
Static_problem(path)
#cross_validate()
Heuristic(path, 0.5, 2) #parameters fine_tuned on average gap