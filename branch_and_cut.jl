using JuMP
using CPLEX

# Données de l'instance
include("data/instance_n5.txt")

# Fonction pour résoudre le problème avec callbacks
function solve_vrp_with_callbacks()
    model = Model(CPLEX.Optimizer)
    MOI.set(model, MOI.NumberOfThreads(), 1)  # Limite à un seul thread pour le callback

    # Désactiver certaines options de CPLEX pour éviter les optimisations automatiques
    set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(model, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_FPHeur", -1)

    # Déclaration des variables
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, u[1:n] >= 0, Int)
    @variable(model, z >= 0)

    ### Contraintes

    # 1 véhicule atteint le client i
    @constraint(model, [i in 2:n], sum(x[j, i] for j in 1:n if j != i) == 1)

    # 1 véhicule quitte le client i
    @constraint(model, [i in 2:n], sum(x[i, j] for j in 1:n if j != i) == 1)

    # Autant de véhicules quittent et atteignent l'entrepôt
    @constraint(model, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n))

    # En arrivant en i le véhicule doit contenir suffisamment 
    # de vaccins pour satisfaire la demande de i
    @constraint(model, [i in 2:n], u[i] <= C - d[i])

    # Si i est le prédécesseur de j, au moins d_i vaccins
    # de plus ont été livrés en arrivant en j qu'en i
    @constraint(model, [i in 2:n, j in 2:n; (i != j)], u[j] - u[i] >= d[i] - C * (1 - x[i, j]))

    # Si i est le premier client de la tournée, aucun
    # vaccin n'a été livré en y arrivant
    @constraint(model, [j in 2:n], u[j] <= C * (1 - x[1, j]))


    # Ensemble de génération égal aux temps non perturbés
    @constraint(model, sum(t[i, j] * x[i, j] for i in 1:n, j in 1:n) <= z)

    # Objectif initial (nominal)
    @objective(model, Min, z)

    # Callback pour ajouter les coupes robustes
    function robust_cut_callback(cb_data, context_id)
        println("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")
        if isIntegerPoint(cb_data, context_id)

            # Charger les valeurs des variables de la relaxation
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_vals = [callback_value(cb_data, x[i, j]) for i in 1:n, j in 1:n]

            # Résolution du sous-problème pour générer un scénario adverse
            model_sp = Model(CPLEX.Optimizer)
            @variable(model_sp, δ1[1:n, 1:n] >= 0)
            @variable(model_sp, δ2[1:n, 1:n] >= 0)

            for i in 1:n, j in 1:n
                @constraint(model_sp, δ1[i, j] <= 1)
                @constraint(model_sp, δ2[i, j] <= 2)
            end

            @constraint(model_sp, sum(δ1[i, j] for i in 1:n, j in 1:n if i != j) <= T)
            @constraint(model_sp, sum(δ2[i, j] for i in 1:n, j in 1:n if i != j) <= T^2)

            @objective(model_sp, Max, sum((t[i, j] + δ1[i, j] * (th[i] + th[j]) + δ2[i, j] * th[i] * th[j]) * x_vals[i, j] for i in 1:n, j in 1:n if i != j))

            optimize!(model_sp)

            δ1_vals = value.(δ1)
            δ2_vals = value.(δ2)
            worst_case_cost = objective_value(model_sp)
            println("value pb slave :", worst_case_cost)

            # Ajout de la coupe si nécessaire
            z_val = callback_value(cb_data, z)
            println("value courante :", z_val)
            if z_val < worst_case_cost
                cstr = @build_constraint(z >= worst_case_cost)
                MOI.submit(model, MOI.LazyConstraint(cb_data), cstr)
                println("Ajout d'une coupe robuste avec coût ", worst_case_cost)
            end
        end
    end

    # Ajouter le callback
    MOI.set(model, CPLEX.CallbackFunction(), robust_cut_callback)

    # Résolution du modèle
    optimize!(model)

    # Affichage des résultats
    if termination_status(model) == MOI.OPTIMAL
        println("Solution optimale trouvée avec coût : ", objective_value(model))
        println("value de x :", value.(x))
        println("value de u :", value.(u))
    else
        println("Solution non optimale trouvée.")
    end
end


function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)
    # context_id == CPX_CALLBACKCONTEXT_CANDIDATE si le callback est
    # appelé dans un des deux cas suivants :
    # cas 1- une solution entière a été obtenue; ou
    # cas 2- une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end
    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)
    # S’il n’y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end


# Lancer la résolution
solve_vrp_with_callbacks()
