include("coupes.jl")
include("dual.jl")
include("problem_static.jl")
include("Heuristic.jl")
include("branch_and_cut.jl")

path = "data/n_5-euclidean_false"

#Cutting_planes(path)
#Dual_solve(path)
#Static_problem(path)
#cross_validate()
start_time = time()
#T, V = Heuristic(path, 0.5, 2) #parameters fine_tuned on average gap
T, V, y, exec= Cutting_planes(path, 20)
#runtime, pb_value = Cutting_planes(path, 20)
println(T - start_time)
println(V)
println(y)
println(exec)