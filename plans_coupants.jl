using JuMP
using CPLEX

# Charger les données à partir du fichier data_loader.jl


# if length(ARGS) == 0
#     println("Veuillez spécifier le chemin du fichier de données en argument. Donnez dans cet ordre : le nombre de noeuds, si les données sont euclidiennes ou non (true/false).")
#     exit(1)
# end

# number_of_nodes = parse(Int, ARGS[1])
# is_euclidean = ARGS[2] == "true"
# path = "./data/n_$(number_of_nodes)-euclidean_$(is_euclidean)"

include("data/instance_n5.txt")

# # Charger les données depuis un fichier texte
# n, t, th, T, d, C = load_data(path)



function resolution_master_pb(n::Int, t::Matrix{Int64}, th::Vector{Int64}, T::Int64, d::Vector{Int64}, C::Int64, list_delta_1, list_delta_2, list_zs)
    verif = false
    master_model = Model(CPLEX.Optimizer)


    # Définir les variables binaires et continues
    @variable(master_model, x[1:n, 1:n], Bin)  # Variables binaires de décision
    @variable(master_model, u[1:n] >= 0)       # Variables de flux u
    @variable(master_model, z)            # Variable pour le coût total

   # Définir les contraintes
    # 1. Chaque noeud (sauf le premier) doit avoir une arrivée
    @constraint(master_model, [i in 2:n], sum(x[i, j] for j in 1:n if j != i) == 1)

    # 2. Chaque noeud (sauf le premier) doit avoir un départ
    @constraint(master_model, [j in 2:n], sum(x[i, j] for i in 1:n if i != j) == 1)

    # 3. Conservation de flux au noeud de départ : pas nécessaire car implicite
    # @constraint(model, sum(x[1, j] for j in 2:n) == sum(x[i, 1] for i in 2:n))


    # 7. Contraintes sur les limites des variables u
    @constraint(master_model, [i in 2:n], u[i] <= C - d[i])
    @constraint(master_model, [i in 2:n], u[i] <= C * (1 - x[1, i]))
    @constraint(master_model, [i in 2:n, j in 2:n, i != j], u[j] - u[i] >= d[i] - C * (1 - x[i, j]))
    @constraint(master_model, z >= sum(t[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j))
    for i in 1:length(list_delta_1)
        delta_1 = list_delta_1[i]
        delta_2 = list_delta_2[i]
        @constraint(master_model, z >= sum((t[i, j] + delta_1[i, j] * (th[i] + th[j]) + delta_2[i, j] * (th[i] * th[j])) * x[i, j] for i in 1:n, j in 1:n if i != j))
    end
    @objective(master_model, Min, z)
    optimize!(master_model)

    status = termination_status(master_model)
    if status == MOI.OPTIMAL
        return value.(x), objective_value(master_model)
    elseif status in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        println("Le modèle est infaisable ou non borné.")
        return false
    else
        println("Aucune solution optimale trouvée. État : ", status)
        return false
    end
end




function resolution_sous_probleme(xs::Matrix{Float64}, n::Int, t::Matrix{Int64}, th::Vector{Int64}, T::Int64)
    # Créer un nouveau modèle avec CPLEX
    model_sous_probleme = Model(CPLEX.Optimizer)

    # Définir les variables binaires et continues
    @variable(model_sous_probleme, delta_1[1:n, 1:n], lower_bound=0, upper_bound=1)  # Variables continues entre 0 et 1
    @variable(model_sous_probleme, delta_2[1:n, 1:n], lower_bound=0, upper_bound=2)  # Variables continues entre 0 et 2

    @objective(model_sous_probleme, Max, sum((t[i, j] + delta_1[i, j] * (th[i] + th[j]) + delta_2[i, j] * (th[i] * th[j])) * xs[i, j] for i in 1:n, j in 1:n if i != j))

    # Définir les contraintes
    @constraint(model_sous_probleme, sum(delta_1[i, j] for i in 1:n, j in 1:n if i != j) <= T)
    @constraint(model_sous_probleme, sum(delta_2[i, j] for i in 1:n, j in 1:n if i != j) <= T*T)

    optimize!(model_sous_probleme)
    zs = objective_value(model_sous_probleme)
    delta_1s = value.(delta_1)
    delta_2s = value.(delta_2)
    return delta_1s, delta_2s, zs
end


function plans_coupants(n::Int, t::Matrix{Int64}, th::Vector{Int64}, T::Int64, d::Vector{Int64}, C::Int64, max_time::Float64)

    list_delta_1 = []
    list_delta_2 = []
    list_zs = []
    xs, z = resolution_master_pb(n, t, th, T, d, C, list_delta_1, list_delta_2, list_zs)
    println("Solution optimale trouvée :")
    println("x = ", xs)
    println("Coût total = ", z)
    delta_1s, delta_2s, zs = resolution_sous_probleme(xs, n, t, th, T)
    list_delta_1 = [delta_1s]
    list_delta_2 = [delta_2s]
    list_zs = [zs]

    println("delta_1 = ", delta_1s)
    println("delta_2 = ", delta_2s)
    println("Coût total = ", zs)
    time_start = time()
    while zs >= z + 1 - 1e-6 && time() - time_start < max_time
        xs, z = resolution_master_pb(n, t, th, T, d, C, list_delta_1, list_delta_2, list_zs)
        println("Solution optimale trouvée :")
        println("x = ", xs)
        println("Coût total = ", z)
        delta_1s, delta_2s, zs = resolution_sous_probleme(xs, n, t, th, T)
        push!(list_delta_1, delta_1s)
        push!(list_delta_2, delta_2s)
        push!(list_zs, zs)
        println("delta_1 = ", delta_1s)
        println("delta_2 = ", delta_2s)
        println("Coût total = ", zs)
    end
    if time() - time_start >= max_time
        println("Temps maximum atteint pour la méthode : plans_coupants")
    else 
        println("Solution optimale trouvée :")
    end
    time_end = time()
    t = time_end - time_start
    println("Temps de résolution : ", t)
    
    return z, t 
end


plans_coupants(n, t, th, T, d, C, 1000.0)