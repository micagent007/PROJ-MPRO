using JuMP
using CPLEX

include("coupes.jl")

function Heuristic(path="data/instance_n5.txt", max_runtime = 60, alpha = 0.5, beta = 2)
    include(path)
    start_time = time()

    # Initialize the solution matrix X
    X = zeros(Int, n, n)

    # Compute the reward matrix based on demand, travel time, and variation time
    reward = zeros(n, n)
    for i in 1:n
        for j in 1:n
            if i != j
                reward[i, j] = d[j] / (alpha * t[i, j] + beta * th[j] + 1e-6)  # Avoid division by zero
            end
        end
    end

    # Track visited nodes and remaining demand
    visited = falses(n)
    visited[1] = true  # Warehouse is always visited
    remaining_demand = sum(d)  # Ensure it's defined in the global scope

    while remaining_demand > 0
        route = [1]  # Start from warehouse
        capacity_used = 0
        current_node = 1

        while true
            # Select the best next node based on reward, ensuring capacity is not exceeded
            next_node = argmax([(reward[current_node, j] * !visited[j] * (capacity_used + d[j] <= C)) for j in 1:n])

            if next_node == 0 || visited[next_node] || capacity_used + d[next_node] > C
                break  # No valid next node, return to warehouse
            end

            # Update tracking variables
            push!(route, next_node)
            X[current_node, next_node] = 1
            visited[next_node] = true
            capacity_used += d[next_node]
            remaining_demand -= d[next_node]

            current_node = next_node
        end

        # Return to warehouse
        push!(route, 1)
        X[current_node, 1] = 1
    end


    _, _, H_value = solve_subproblem(X, n, t, th, T, max_runtime)

    end_time = time()
    exec_time = end_time - start_time
    
    return end_time, H_value, X, exec_time
end

function cross_validate()
    alpha_values = [0.1, 0.5, 1, 2, 5, 10, 50]
    beta_values = [0.1, 0.5, 1, 2, 5, 50, 200, 500, 1000]
    optimal_values = Dict("n_5-euclidean_false" => 2786, "n_6-euclidean_false" => 3307, "n_7-euclidean_false" => 2066, 
                          "n_8-euclidean_false" => 3270, "n_9-euclidean_false" => 3113, "n_10-euclidean_false" => 2684, 
                          "n_6-euclidean_true" => 3087, "n_7-euclidean_true" => 3132, "n_8-euclidean_true" => 3776, 
                          "n_9-euclidean_true" => 4181, "n_10-euclidean_true" => 5423)
    
    best_alpha, best_beta = 1, 1
    min_avg_gap = Inf

    for alpha in alpha_values
        for beta in beta_values
            total_gap = 0
            count = 0
            
            for (instance, optimal) in optimal_values
                _, heuristic_value = Heuristic("data/$(instance)", alpha, beta)
                gap = (heuristic_value - optimal) / optimal * 100
                total_gap += gap
                count += 1
            end
            
            avg_gap = total_gap / count
            if avg_gap < min_avg_gap
                min_avg_gap = avg_gap
                best_alpha, best_beta = alpha, beta
            end
        end
    end

    println("Best alpha: ", best_alpha, ", Best beta: ", best_beta)
    println("Minimum average gap: ", min_avg_gap)
end

function cross_validate_max_gap()
    alpha_values = [0.1, 0.5, 1, 2, 5, 10, 50]
    beta_values = [0.1, 0.5, 1, 2, 5, 50, 200, 500, 1000]
    optimal_values = Dict("n_5-euclidean_false" => 2786, "n_6-euclidean_false" => 3307, "n_7-euclidean_false" => 2066, 
                          "n_8-euclidean_false" => 3270, "n_9-euclidean_false" => 3113, "n_10-euclidean_false" => 2684, 
                          "n_6-euclidean_true" => 3087, "n_7-euclidean_true" => 3132, "n_8-euclidean_true" => 3776, 
                          "n_9-euclidean_true" => 4181, "n_10-euclidean_true" => 5423)
    
    best_alpha, best_beta = 1, 1
    min_max_gap = Inf

    for alpha in alpha_values
        for beta in beta_values
            max_gap = -Inf
            
            for (instance, optimal) in optimal_values
                _, heuristic_value = Heuristic("data/$(instance)", alpha, beta)
                gap = (heuristic_value - optimal) / optimal * 100
                max_gap = max(max_gap, gap)
            end
            
            if max_gap < min_max_gap
                min_max_gap = max_gap
                best_alpha, best_beta = alpha, beta
            end
        end
    end

    println("Best alpha: ", best_alpha, ", Best beta: ", best_beta)
    println("Minimum maximum gap: ", min_max_gap)
end
