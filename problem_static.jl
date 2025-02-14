using JuMP
using CPLEX

function Static_problem(path="data/instance_n5.txt", max_runtime=60.0)
    include(path)

    m = Model(CPLEX.Optimizer)
    set_silent(m)

    # Définir la limite de temps
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", max_runtime)

    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)

    ### Variables de décision
    @variable(m, x[1:n, 1:n], Bin)
    @variable(m, u[1:n] >= 0, Int)

    ### Fonction objectif
    @objective(m, Min, sum(t[i, j] * x[i, j] for i in 1:n, j in 1:n if j != i))

    ### Contraintes
    @constraint(m, [i in 2:n], sum(x[j, i] for j in 1:n if j != i) == 1)
    @constraint(m, [i in 2:n], sum(x[i, j] for j in 1:n if j != i) == 1)
    @constraint(m, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n))
    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n, j in 2:n; (i != j)], u[j] - u[i] >= d[i] - C * (1 - x[i, j]))
    @constraint(m, [j in 2:n], u[j] <= C * (1 - x[1, j]))

    start_time = time()
    optimize!(m)
    end_time = time()
    exec_time = end_time - start_time

    best_obj = nothing
    best_x = nothing
    best_u = nothing

    if has_values(m)
        best_obj = objective_value(m)
        best_x = value.(x)
        best_u = value.(u)
        #println("Solution trouvée en ", max_runtime, " secondes.")
        #println("x = ", best_x)
        #println("u = ", best_u)
        println("Valeur de l'objectif : ", best_obj)
    else
        best_obj = "Inf"
        best_x = "Pas de solution trouvée"
        println("Aucune solution trouvée dans le temps imparti.")
    end

    return end_time, best_obj, best_x, exec_time
end