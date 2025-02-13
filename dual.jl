using JuMP
using CPLEX

function Dual_solve(path="data/instance_n5.txt", max_runtime=60.0)
    include(path)

    m = Model(CPLEX.Optimizer)
    set_time_limit_sec(m, max_runtime)  # Définir la limite de temps
    set_silent(m)

    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", timeout)

    ### Variables de décision
    @variable(m, x[1:n, 1:n], Bin)
    @variable(m, u[2:n] >= 0, Int)
    @variable(m, alpha1 >= 0)
    @variable(m, alpha2 >= 0)
    @variable(m, beta1[1:n, 1:n] >= 0)
    @variable(m, beta2[1:n, 1:n] >= 0)

    ### Fonction objectif
    @objective(m, Min, alpha1 * T + alpha2 * T^2 + sum(t[i, j] * x[i, j] + beta1[i, j] + 2 * beta2[i, j] for i in 1:n, j in 1:n if j != i))

    ### Contraintes
    @constraint(m, [i in 2:n], sum(x[j, i] for j in 1:n if j != i) == 1)
    @constraint(m, [i in 2:n], sum(x[i, j] for j in 1:n if j != i) == 1)
    @constraint(m, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n))
    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n, j in 2:n; (i != j)], u[j] - u[i] >= d[i] - C * (1 - x[i, j]))
    @constraint(m, [j in 2:n], u[j] <= C * (1 - x[1, j]))

    for i in 1:n, j in 1:n
        if i != j
            @constraint(m, alpha1 + beta1[i, j] >= (th[i] + th[j]) * x[i, j])
            @constraint(m, alpha2 + beta2[i, j] >= th[i] * th[j] * x[i, j])
        end
    end

    start_time = time()
    optimize!(m)
    end_time = time()
    exec_time = end_time - start_time

    best_objective = nothing
    if has_values(m)
        best_objective = objective_value(m)
        vX = JuMP.value.(x)
        vU = JuMP.value.(u)
        #println("x = ", vX)
        #println("u = ", vU)
        println("Meilleure valeur de l'objectif : ", best_objective)
    else
        best_objective = "Inf"
        vX = "Pas de solution trouvée"
        println("Aucune solution trouvée dans le temps imparti.")
    end

    return time(), best_objective, vX, exec_time
end
