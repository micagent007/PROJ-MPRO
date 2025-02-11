using JuMP
using CPLEX

include("coupes.jl")


function Heuristic(path="data/instance_n5.txt", alpha = 1, beta = 1)
    include(path)

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


    println("Optimized tour matrix X:")
    println(X)

    H_value = solve_subproblem(X, n, t, th)

    return time(), H_value
end