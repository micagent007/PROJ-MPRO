include("coupes.jl")
include("branch_and_cut.jl")
include("dual.jl")
include("problem_static.jl")
include("Heuristic.jl")

"""
Résout les instances de data et renvoie la liste des temps de résolution des instances
"""
function testing(data="./data", timeout=10)
    # start_time = time()
    # list_time = []
    fout = open("./results/output.txt", "w")
    println(fout, "Instances,Static_sol,Static_val,Static_bound,Static_time,BC_sol,BC_val,BC_bound,BC_time,Dual_sol,Dual_val,Dual_bound,Dual_time,Heuristic_sol,Heuristic_val,Heuristic_bound,Heuristic_time")
    for inst in readdir(data)
    # inst = "n_5-euclidean_false"
        instance = data * "/" * inst
        print(fout, instance, ",")
        static = Static_problem(instance, timeout)
        print(fout, static, ",")
        # pc = Cutting_planes(instance, timeout)
        # print(fout, pc, ",")
        bb = solve_vrp_with_callbacks(instance, timeout)
        print(fout, bb, ",")
        dual = Dual_solve(instance, timeout)
        print(fout, dual)
            # end_time = solve(instance,timeout)
            # push!(list_time, end_time - start_time)
        print(fout, "\n")
    end
    # return list_time
    close(fout)
end

testing()