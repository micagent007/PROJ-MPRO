"""
Résout les instances de data et renvoie la liste des temps de résolution des instances
"""
function testing(method, data="./data", timeout=3600)
    start_time = time()
    list_time = []
    for instance in readdir(data)
        end_time = time()
        push!(list_time, end_time - start_time)
    end
    return list_time
end