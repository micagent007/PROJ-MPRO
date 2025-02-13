using JuMP
using CPLEX

# Charger les données depuis le fichier instance_n5.txt
# include("data/instance_n5.txt")


function solve_main_problem(nb_nodes::Int, time_matrix::Matrix{Int64}, node_weights::Vector{Int64}, max_time::Int64, demand::Vector{Int64}, capacity::Int64, delta1_list, delta2_list, cost_list)
    is_feasible = false
    main_model = Model(CPLEX.Optimizer)

    # Définition des variables de décision
    @variable(main_model, y[1:nb_nodes, 1:nb_nodes], Bin)
    @variable(main_model, flow[1:nb_nodes] >= 0)
    @variable(main_model, total_cost)

    # Contraintes
    @constraint(main_model, [i in 2:nb_nodes], sum(y[i, j] for j in 1:nb_nodes if j != i) == 1)
    @constraint(main_model, [j in 2:nb_nodes], sum(y[i, j] for i in 1:nb_nodes if i != j) == 1)
    @constraint(main_model, [i in 2:nb_nodes], flow[i] <= capacity - demand[i])
    @constraint(main_model, [i in 2:nb_nodes], flow[i] <= capacity * (1 - y[1, i]))
    @constraint(main_model, [i in 2:nb_nodes, j in 2:nb_nodes, i != j], flow[j] - flow[i] >= demand[i] - capacity * (1 - y[i, j]))
    @constraint(main_model, total_cost >= sum(time_matrix[i, j] * y[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j))
    
    for idx in 1:length(delta1_list)
        delta1 = delta1_list[idx]
        delta2 = delta2_list[idx]
        @constraint(main_model, total_cost >= sum((time_matrix[i, j] + delta1[i, j] * (node_weights[i] + node_weights[j]) + delta2[i, j] * (node_weights[i] * node_weights[j])) * y[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j))
    end
    
    @objective(main_model, Min, total_cost)
    optimize!(main_model)

    solution_status = termination_status(main_model)
    if solution_status == MOI.OPTIMAL
        return value.(y), objective_value(main_model)
    elseif solution_status in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        println("Modèle infaisable ou non borné.")
        return false
    else
        println("Pas de solution optimale trouvée. État : ", solution_status)
        return false
    end
end


function solve_subproblem(ys, nb_nodes::Int, time_matrix::Matrix{Int64}, node_weights::Vector{Int64}, max_time::Int64 = 10^6)
    subproblem_model = Model(CPLEX.Optimizer)

    # Définition des variables
    @variable(subproblem_model, delta1[1:nb_nodes, 1:nb_nodes], lower_bound=0, upper_bound=1)
    @variable(subproblem_model, delta2[1:nb_nodes, 1:nb_nodes], lower_bound=0, upper_bound=2)

    @objective(subproblem_model, Max, sum((time_matrix[i, j] + delta1[i, j] * (node_weights[i] + node_weights[j]) + delta2[i, j] * (node_weights[i] * node_weights[j])) * ys[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j))

    # Contraintes
    @constraint(subproblem_model, sum(delta1[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j) <= max_time)
    @constraint(subproblem_model, sum(delta2[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j) <= max_time * max_time)

    optimize!(subproblem_model)
    obj_value = objective_value(subproblem_model)
    delta1_values = value.(delta1)
    delta2_values = value.(delta2)
    return delta1_values, delta2_values, obj_value
end


function Cutting_planes(path = "data/instance_n5.txt", max_runtime = 10^6)
    
    include(path)

    """println("n = ", n)
    println("t = ", t)
    println("th = ", th)
    println("T = ", T)
    println("d = ", d)
    println("C = ", C)"""

    nb_nodes = n
    time_matrix = t
    node_weights = th
    max_time = T
    demand = d
    capacity = C

    delta1_list = []
    delta2_list = []
    cost_values = []
    ys, z_value = solve_main_problem(nb_nodes, time_matrix, node_weights, max_time, demand, capacity, delta1_list, delta2_list, cost_values)
    println("Solution optimale initiale trouvée :")
    println("y = ", ys)
    println("Coût total = ", z_value)
    delta1_values, delta2_values, obj_cost = solve_subproblem(ys, nb_nodes, time_matrix, node_weights, max_time)
    push!(delta1_list, delta1_values)
    push!(delta2_list, delta2_values)
    push!(cost_values, obj_cost)
    
    println("delta1 = ", delta1_values)
    println("delta2 = ", delta2_values)
    println("Coût total = ", obj_cost)
    
    start_time = time()
    while obj_cost >= z_value + 1 - 1e-6 && time() - start_time < max_runtime
        ys, z_value = solve_main_problem(nb_nodes, time_matrix, node_weights, max_time, demand, capacity, delta1_list, delta2_list, cost_values)
        println("Nouvelle solution optimale trouvée :")
        println("y = ", ys)
        println("Coût total = ", z_value)
        delta1_values, delta2_values, obj_cost = solve_subproblem(ys, nb_nodes, time_matrix, node_weights, max_time)
        push!(delta1_list, delta1_values)
        push!(delta2_list, delta2_values)
        push!(cost_values, obj_cost)
        println("delta1 = ", delta1_values)
        println("delta2 = ", delta2_values)
        println("Coût total = ", obj_cost)
    end
    
    if time() - start_time >= max_runtime
        println("Temps limite atteint pour la méthode des plans coupants")
    else
        println("Solution optimale finale trouvée.")
    end
    
    end_time = time()
    exec_time = end_time - start_time
    println("Temps de résolution : ", exec_time)
    println("Valeur optimale", z_value)
    
    return time(), z_value
end
