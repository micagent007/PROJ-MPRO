using DataFrames, XLSX
using Dates
using ProgressBars

include("coupes.jl")
include("dual.jl")
include("problem_static.jl")
include("Heuristic.jl")
include("branch_and_cut.jl")

directory = "data/"
max_runtime = 60  # Temps maximum en secondes

# Récupération des instances
instances = readdir(directory)

# Dictionnaire pour stocker les résultats par méthode
results = Dict(
    "Cutting Planes" => [],
    "Branch & Cut" => [],
    "Dual" => [],
    "Heuristic" => [],
    "Static" => []
)

data = []
for instance in ProgressBar(instances)
    path = joinpath(directory, instance)
    println("Résolution de l'instance : ", instance)

    #_, cp_best, cp_sol, cp_time = Cutting_planes(path, max_runtime)
    #_, bc_best, bc_sol, bc_time = solve_vrp_with_callbacks(path, max_runtime)
    #_, dual_best, dual_sol, dual_time = Dual_solve(path, max_runtime)
    _, heur_best, heur_sol, heur_time = Heuristic(path, max_runtime, 0.1, 1, 5)
    #_, static_best, static_sol, static_time = Static_problem(path, max_runtime)

    # Stockage des résultats dans le dictionnaire
    #push!(results["Cutting Planes"], "$instance : $cp_sol")
    #push!(results["Branch & Cut"], "$instance : $bc_sol")
    #push!(results["Dual"], "$instance : $dual_sol")
    push!(results["Heuristic"], "$instance : $heur_sol")
    #push!(results["Static"], "$instance : $static_sol")

    #Calcul de la meilleure sol robuste
    """Best_list = [cp_best, dual_best, bc_best, Inf]
    true_list = []
    for sol in Best_list
        if typeof(sol) != String
            push!(true_list, sol)
        end
    end
    println(true_list)
    best_sol = minimum(true_list)

    # Calcul du gap
    gap = static_best != 0 ? ((best_sol - static_best) / static_best) * 100 : Inf"""

    push!(data, (instance=instance, heur_time=heur_time, heur_best=heur_best))
end

# Écriture des fichiers texte
for (method, lines) in results
    open("Solutions/$(method).txt", "w") do file
        write(file, join(lines, "\n"))
    end
end

# Création du DataFrame pour les résultats
df = DataFrame(data)  # Utilisation de NamedTuples pour créer le DataFrame

# Écriture des résultats dans un fichier Excel avec plusieurs feuilles
output_file = "Solutions/resultats.xlsx"
XLSX.writetable(output_file, "Résultats" => df)

println("Fichier Excel généré: ", output_file)
println("Fichiers texte générés pour chaque méthode.")