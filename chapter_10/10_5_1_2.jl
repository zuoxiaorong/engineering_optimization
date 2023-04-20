# using DP method solving network design with relays

function data()
    N = 1:9
    r = [0, 7, 3, 5, 4, 1, 3, 2, 4]
    d = [2, 1, 1, 3, 1, 2, 2, 2, 0]
    L = 5
    return N, r, d, L
end

function far_node(N, d, L)
    g = []
    for i in 1:length(N)
        cur_i = i
        cur_dis = L
        while cur_dis >= 0 && cur_i <= length(N)
            cur_dis -= d[cur_i]
            cur_i += 1
            if cur_i == length(N) && cur_dis > 0
                cur_i =  length(N) + 2
            end
        end
        if i == length(N)
            push!(g, i + 1)
        else
            push!(g, cur_i - 1)
        end
    end
    return g
end

function DP_networkdesign(N, r, d, L, g, startNode)
    f = Dict()
    f[length(N)+1] = Dict(
        "v" => 0,
        "q" => nothing
    )
    currentNode = length(N)

    while currentNode >= 1
        relay = [r[currentNode] + (currentNode == i ? 0 : f[i]["v"]) for i in currentNode+1:g[currentNode]]
        f[currentNode] = Dict(
            "v" => minimum(relay),
            "q" => findall(x -> x == minimum(relay), relay) .+ currentNode
        )
        currentNode -=1
    end
    total_min_cost = f[startNode]["v"]
    solution = get_solution(f, startNode, N, [])
    num = findall(x -> x == startNode, solution)
    solutions = [i < length(num) ? collect(solution[num[i]:num[i+1]-1]) : collect(solution[num[i]:length(solution)]) for i in 1:length(num)]
    return solutions, total_min_cost
end

function get_solution(f, startNode, N, set)
    push!(set, startNode)

    sets = []
    for i in f[startNode]["q"]
        if !(i in set)
            if i == length(N) + 1
                return set
                break
            end
            new_sets = get_solution(f, i, N, copy(set))
            append!(sets, new_sets)
        end
    end
    return sets

end
    
function main()
    N, r, d, L = data()
    g = far_node(N, d, L)
    startNode = 1
    solutions, total_min_cost = DP_networkdesign(N, r, d, L, g, startNode)
end
main()
# solution 
# total_min_cost = 6
# solutions = [1, 3, 6, 8]