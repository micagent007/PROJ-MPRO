using JuMP
using CPLEX

path = "/Users/guillaume/Documents/MPRO/1.Tronc Commun/PROJ/data/instance_n5.txt"

include(path)

println("n = ", n)
println("t = ", t)
println("th = ", th)
println("T = ", T)
println("d = ", d)
println("C = ", C)


m = Model(CPLEX.Optimizer)

### Variables de décision

@variable(m, x[1:n, 1:n], Bin)

@variable(m, u[2:n] >= 0, Int)

### Fonction objectif

@objective(m, Min, sum(t[i,j] * x[i,j] for i in 1:n, j in 1:n if j != i))

### Contraintes

# 1 véhicule atteint le client i
@constraint(m, [i in 2:n], sum(x[j,i] for j in 1:n if j != i) == 1)

# 1 véhicule quitte le client i
@constraint(m, [i in 2:n], sum(x[i,j] for j in 1:n if j != i) == 1)

# Autant de véhicules quittent et atteignent l'entrepôt
@constraint(m, sum(x[1,j] for j in 2:n) == sum(x[j,1] for j in 2:n))

# En arrivant en i le véhicule doit contenir suffisamment 
# de vaccins pour satisfaire la demande de i
@constraint(m, [i in 2:n], u[i] <= C - d[i])

# Si i est le prédécesseur de j, au moins d_i vaccins
# de plus ont été livrés en arrivant en j qu'en i
@constraint(m, [i in 2:n, j in 2:n; (i != j)], u[j] - u[i] >= d[i] - C * (1 - x[i,j]))

# Si i est le premier client de la tournée, aucun
# vaccin n'a été livré en y arrivant
@constraint(m, [j in 2:n], u[j] <= C * (1 - x[1,j]))

optimize!(m)

if primal_status(m) == MOI.FEASIBLE_POINT
    vX = JuMP.value.(x)
    vU = JuMP.value.(u)
    println("x = ", vX)
    println("u = ", vU)
    println("Valeur de l'objectif : ", JuMP.objective_value(m))
end
