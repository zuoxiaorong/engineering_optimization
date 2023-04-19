# Dynamic Programming from back to front

function data()
    T = 1:6
    s = 10
    h = 0.5
    d = [10, 15, 8, 12, 20, 16]
    return T, s, h, d
end

function get_solution(prodDict, T, d)
    min_cost = prodDict[1]["v"]
    current_stage = 1
    prod_stage = [current_stage]
    while prodDict[current_stage]["q"] <= length(T)
        push!(prod_stage, prodDict[current_stage]["q"])
        current_stage = prod_stage[end]
    end
    prod = zeros(length(T))
    for i in 1:length(prod_stage)
        if i < length(prod_stage)
            prod[prod_stage[i]] = sum(d[t] for t in prod_stage[i]:prod_stage[i+1]-1)
        else
            prod[prod_stage[i]] = sum(d[t] for t in prod_stage[i]:length(T))
        end
    end
   return min_cost, prod_stage, prod
end

function dp_lotsize(T, s, h, d)
    prodDict = Dict()
    currentStage = length(T)
    prodDict[currentStage + 1] = Dict(
        "v" => 0,
        "q" => 0
    ) 

    while currentStage >= 1
        cost = [s + sum(h*d[tt]*(tt - currentStage) for tt in currentStage:t) + prodDict[t+1]["v"] for t in currentStage:length(T)]
        prodDict[currentStage] = Dict(
            "v" => minimum(cost),
            "q" => findfirst(x -> x == minimum(cost), cost) + currentStage
        )
        currentStage -= 1 
    end
    cost, stage, prod = get_solution(prodDict, T, d)
end

function main()
    T, s, h, d = data()
    cost, stage, prod = dp_lotsize(T, s, h, d)
end
main()
# solution
# cost = 51.5
# order_stage = [1, 3, 5]
# lotsizing = [25.0, 0.0, 20.0, 0.0, 36.0, 0.0]