# using DP method solving the EV charging problem with fixed visiting routing

function data()
    N = 1:10
    D = [45, 35, 40, 55, 20, 50, 30, 40, 30, 50, 0]
    d1 = [3, 6, 8, 10, 1, 48, 28, 25, 15, 5]
    d2 = [46, 34, 39, 54, 21, 5, 8, 25, 20, 48]
    L = 100   
    return N, D, d1, d2, L
end

function get_solution(routeDict, g, N)
    min_dis = routeDict[0]["h"]
    path = [routeDict[0]["q"]]
    start_arc = path[end]
    while g[start_arc] < length(N)
        push!(path, routeDict[start_arc]["q"])
        start_arc = path[end]
    end
    return min_dis, path
end

function dp_evtour(N, D, d1, d2, L)
    routeDict = Dict()
    currentArc = length(N)
    routeDict[currentArc + 1] = Dict(
        "h" => 0,
        "q" => 0
    )

    #define g
    g = zeros(Int, length(N))
    for i in N
        if i < length(N)
            current_dis = d2[i]        
            current_i = i + 1
            while current_dis <= L && current_i < length(N)
                if current_dis + D[current_i] + d1[current_i+1] <= L
                    current_dis += D[current_i]
                else
                    break
                end
                current_i += 1
            end
            g[i] = current_i
            
            if d2[i] + sum(D[ii] for ii in (i+1):length(N)) <= L
                g[i] = length(N) + 1
            end
        elseif i == length(N)
            g[i] = length(N) + 1
        end
    end

    while currentArc >= 0
        if currentArc >= 1
            if g[currentArc] > length(N)
                dis = d1[currentArc] + d2[currentArc] + (currentArc == length(N) ? 0 : sum(D[ii] for ii in (currentArc+1):length(N)))
                routeDict[currentArc] = Dict(
                    "h" => dis,
                    "q" => currentArc
                )
            else
                dis = [d1[currentArc] + d2[currentArc] + routeDict[i+1]["h"] + (i == currentArc ? 0 : sum(D[ii] for ii in (currentArc+1):i)) for i in currentArc:g[currentArc]-1]
                routeDict[currentArc] = Dict(
                    "h" => minimum(dis),
                    "q" => findfirst(x -> x == minimum(dis), dis) + currentArc
                )
            end
        elseif currentArc == 0
            current_node = currentArc + 1
            dis = [routeDict[i]["h"] + (i == current_node ? 0 : sum(D[ii] for ii in 1:i-1)) for i in current_node:g[current_node]]
            routeDict[currentArc] = Dict(
                "h" => minimum(dis),
                "q" => findfirst(x -> x == minimum(dis), dis)
            )
        end
        currentArc -= 1
    end
    return get_solution(routeDict, g, N)
end

function main()
    N, D, d1, d2, L = data()
    distance, path = dp_evtour(N, D, d1, d2, L)
end
main()
# solution 
# min_distance = 412
# charging_node = [3, 5, 6, 9]
