using JuMP
using CPLEX

# Charger les données depuis le fichier instance_n5.txt
# include("data/instance_n5.txt")

function solve_main_problem(nb_nodes::Int, time_matrix::Matrix{Int64}, node_weights::Vector{Int64}, max_time::Int64, demand::Vector{Int64}, capacity::Int64, delta1_list, delta2_list, cost_list, max_runtime)
    main_model = Model(CPLEX.Optimizer)
    set_time_limit_sec(main_model, max_runtime)
    set_silent(main_model)

    @variable(main_model, y[1:nb_nodes, 1:nb_nodes], Bin)
    @variable(main_model, flow[1:nb_nodes] >= 0)
    @variable(main_model, total_cost)

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

    if has_values(main_model)
        return value.(y), objective_value(main_model)
    else
        return "Pas de solution trouvée", Inf
    end
end

function solve_subproblem(ys, nb_nodes::Int, time_matrix::Matrix{Int64}, node_weights::Vector{Int64}, max_time::Int64, max_runtime)
    subproblem_model = Model(CPLEX.Optimizer)
    set_time_limit_sec(subproblem_model, max_runtime)
    set_silent(subproblem_model)

    @variable(subproblem_model, delta1[1:nb_nodes, 1:nb_nodes], lower_bound=0, upper_bound=1)
    @variable(subproblem_model, delta2[1:nb_nodes, 1:nb_nodes], lower_bound=0, upper_bound=2)
    
    @objective(subproblem_model, Max, sum((time_matrix[i, j] + delta1[i, j] * (node_weights[i] + node_weights[j]) + delta2[i, j] * (node_weights[i] * node_weights[j])) * ys[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j))
    
    @constraint(subproblem_model, sum(delta1[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j) <= max_time)
    @constraint(subproblem_model, sum(delta2[i, j] for i in 1:nb_nodes, j in 1:nb_nodes if i != j) <= max_time * max_time)
    
    optimize!(subproblem_model)
    
    if has_values(subproblem_model)
        return value.(delta1), value.(delta2), objective_value(subproblem_model)
    else
        return "Pas de solution trouvée", "Pas de solution trouvée", Inf
    end
end

function Cutting_planes(path = "data/instance_n5.txt", max_runtime = 60)
    include(path)
    
    nb_nodes = n
    time_matrix = t
    node_weights = th
    max_time = T
    demand = d
    capacity = C

    delta1_list = []
    delta2_list = []
    cost_values = []
    ys, z_value = solve_main_problem(nb_nodes, time_matrix, node_weights, max_time, demand, capacity, delta1_list, delta2_list, cost_values, max_runtime)
    
    if z_value == Inf
        println("Aucune solution trouvée pour le problème maître dans le temps imparti.")
        return time(), "Inf", "Pas de solution trouvée", 0
    end
    
    delta1_values, delta2_values, obj_cost = solve_subproblem(ys, nb_nodes, time_matrix, node_weights, max_time, max_runtime)
    push!(delta1_list, delta1_values)
    push!(delta2_list, delta2_values)
    push!(cost_values, obj_cost)
    
    start_time = time()
    while obj_cost >= z_value + 1 - 1e-6 && time() - start_time < max_runtime
        new_ys, new_z_value = solve_main_problem(nb_nodes, time_matrix, node_weights, max_time, demand, capacity, delta1_list, delta2_list, cost_values, max_runtime)
        
        if new_z_value == Inf
            break
        else
            z_value = new_z_value
            ys = new_ys
        end
        
        delta1_values, delta2_values, obj_cost = solve_subproblem(ys, nb_nodes, time_matrix, node_weights, max_time, max_runtime)
        push!(delta1_list, delta1_values)
        push!(delta2_list, delta2_values)
        push!(cost_values, obj_cost)
        println(obj_cost)
    end
    
    if time() - start_time >= max_runtime
        println("Temps limite atteint pour la méthode des plans coupants")
    else
        println("Solution optimale finale trouvée.")
    end
    
    end_time = time()
    exec_time = end_time - start_time
    println("Temps de résolution : ", exec_time)
    println("Valeur optimale : ", z_value)
    
    return end_time, z_value, ys, exec_time
end
