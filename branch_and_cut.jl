using JuMP
using CPLEX

# Fonction pour résoudre le problème avec callbacks et une limite de temps
function solve_vrp_with_callbacks(path="data/n_6-euclidean_false", max_runtime=60.0)
    include(path)

    model = Model(CPLEX.Optimizer)
    MOI.set(model, MOI.NumberOfThreads(), 1)  # Limite à un seul thread pour le callback
    set_optimizer_attribute(model, "CPXPARAM_TimeLimit", max_runtime)  # Définir la limite de temps
    set_silent(model)

    # Désactiver certaines options de CPLEX pour éviter les optimisations automatiques
    set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(model, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    
    # Déclaration des variables
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, u[1:n] >= 0, Int)
    @variable(model, z >= 0)

    ### Contraintes
    @constraint(model, [i in 2:n], sum(x[j, i] for j in 1:n if j != i) == 1)
    @constraint(model, [i in 2:n], sum(x[i, j] for j in 1:n if j != i) == 1)
    @constraint(model, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n))
    @constraint(model, [i in 2:n], u[i] <= C - d[i])
    @constraint(model, [i in 2:n, j in 2:n; (i != j)], u[j] - u[i] >= d[i] - C * (1 - x[i, j]))
    @constraint(model, [j in 2:n], u[j] <= C * (1 - x[1, j]))
    @constraint(model, sum(t[i, j] * x[i, j] for i in 1:n, j in 1:n) <= z)

    # Objectif initial (nominal)
    @objective(model, Min, z)

    # Callback pour ajouter les coupes robustes
    function robust_cut_callback(cb_data, context_id)
        if isIntegerPoint(cb_data, context_id)
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_vals = callback_value.(cb_data, x)

            model_sp = Model(CPLEX.Optimizer)
            @variable(model_sp, 0 <= δ1[1:n, 1:n] <= 1)
            @variable(model_sp, 0 <= δ2[1:n, 1:n] <= 2)

            @constraint(model_sp, sum(δ1[i, j] for i in 1:n, j in 1:n if i != j) <= T)
            @constraint(model_sp, sum(δ2[i, j] for i in 1:n, j in 1:n if i != j) <= T^2)

            @objective(model_sp, Max, sum((t[i, j] + δ1[i, j] * (th[i] + th[j]) + δ2[i, j] * th[i] * th[j]) * x_vals[i, j] for i in 1:n, j in 1:n if i != j))

            optimize!(model_sp)

            worst_case_cost = objective_value(model_sp)
            z_val = callback_value(cb_data, z)
            if z_val < worst_case_cost - 1e-4
                cstr = @build_constraint(z >= sum((t[i, j] + value(δ1[i, j]) * (th[i] + th[j]) + value(δ2[i, j]) * th[i] * th[j]) * x[i, j] for i in 1:n, j in 1:n if i != j))
                MOI.submit(model, MOI.LazyConstraint(cb_data), cstr)
            end
        end
    end

    MOI.set(model, CPLEX.CallbackFunction(), robust_cut_callback)
    start_time = time()
    optimize!(model)
    end_time = time()
    exec_time = end_time - start_time

    # Vérification de l'état de la solution
    if has_values(model)
        best_value = objective_value(model)
        best_solution = value.(x)
    else
        best_value = "Inf"
        best_solution = "Pas de solution trouvée"
    end
    
    println("Temps de fin : ", end_time, " Meilleure valeur trouvée : ", best_value)
    return end_time, best_value, best_solution, exec_time
end

function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)
    return ret == 0 && ispoint_p[] == 1
end

# Lancer la résolution avec une limite de temps
solve_vrp_with_callbacks("data/n_6-euclidean_false", 300.0)
