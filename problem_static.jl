using JuMP
using CPLEX


function Static_problem(path="data/n_5-euclidean_false", timeout=10)

    include(path)

    start_time = time()

    """println("n = ", n)
    println("t = ", t)
    println("th = ", th)
    println("T = ", T)
    println("d = ", d)
    println("C = ", C)"""


    m = Model(CPLEX.Optimizer)

    set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur", -1)
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", timeout)

    ### Variables de décision

    @variable(m, x[1:n, 1:n], Bin)

    @variable(m, u[1:n] >= 0, Int)

    ### Fonction objectif

    @objective(m, Min, sum(t[i, j] * x[i, j] for i in 1:n, j in 1:n if j != i))

    ### Contraintes

    # 1 véhicule atteint le client i
    @constraint(m, [i in 2:n], sum(x[j, i] for j in 1:n if j != i) == 1)

    # 1 véhicule quitte le client i
    @constraint(m, [i in 2:n], sum(x[i, j] for j in 1:n if j != i) == 1)

    # Autant de véhicules quittent et atteignent l'entrepôt
    @constraint(m, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n))

    # En arrivant en i le véhicule doit contenir suffisamment 
    # de vaccins pour satisfaire la demande de i
    @constraint(m, [i in 2:n], u[i] <= C - d[i])

    # Si i est le prédécesseur de j, au moins d_i vaccins
    # de plus ont été livrés en arrivant en j qu'en i
    @constraint(m, [i in 2:n, j in 2:n; (i != j)], u[j] - u[i] >= d[i] - C * (1 - x[i, j]))

    # Si i est le premier client de la tournée, aucun
    # vaccin n'a été livré en y arrivant
    @constraint(m, [j in 2:n], u[j] <= C * (1 - x[1, j]))

    optimize!(m)

    # if primal_status(m) == MOI.FEASIBLE_POINT
    #     vX = JuMP.value.(x)
    #     vU = JuMP.value.(u)
    #     println("x = ", vX)
    #     println("u = ", vU)
    #     println("Valeur de l'objectif : ", JuMP.objective_value(m))
    # end

    return JuMP.value.(x), JuMP.objective_value(m), JuMP.objective_bound(m), time() - start_time
end

Static_problem()