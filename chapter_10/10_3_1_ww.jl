# W-W Algorithm (proposed by H. M. Wagner and T. M. Whitin in 1958)

function data()
    T = 1:6
    s = 10
    h = 0.5
    d = [10, 15, 8, 12, 20, 16]
    transfer_t = [Int(floor((s/d[i])/h)) for i in T]
    return T, s, h, d, transfer_t
end

function get_solution(orderDict, T, d)
    min_cost = minimum(values(filter(x -> x[1][1] == collect(T)[end], orderDict)))
    prod = zeros(length(T))
    start_month = length(T)
    while prod[begin] == 0
        index = [key[2] for key in keys(orderDict) if key[1] == start_month && orderDict[key] == minimum(values(filter(x -> x[1][1] == start_month, orderDict)))]
        for i in index
            prod[i] = sum(d[i:start_month])
        end
        start_month = minimum(index) - 1
    end
    min_cost, prod 
end

function W_W(T, s, h, d, pervT)
    orderDict = Dict()

    startMonth = collect(T) .- pervT
    monthSet = []
    for i in T
        if i == 1
            push!(monthSet, [i])
        elseif i >= 2 && startMonth[i] <= 0
            push!(monthSet, collect(minimum(monthSet[i-1]):i))
        elseif i >= 2 && startMonth[i] > 0
            push!(monthSet, collect(max(startMonth[i], minimum(monthSet[i-1])):i))
        end
    end

    currentMonth = 1
    while currentMonth <= length(T)
        if currentMonth == 1
            orderDict[(currentMonth, 1)] = s[currentMonth]
        else
            for t in monthSet[currentMonth]
                if t == currentMonth
                    orderDict[(t,t)] = minimum(values(filter(x -> x[1][1] == currentMonth - 1, orderDict))) + s
                else
                    orderDict[(currentMonth, t)] = orderDict[(currentMonth-1, t)] + d[currentMonth]*(currentMonth - t)*h
                end
            end
        end
        currentMonth += 1
    end
    cost, solution = get_solution(orderDict, T, d)
    return cost, solution
end

function main()
    T, s, h, d, pervT = data()
    cost, solution = W_W(T, s, h, d, pervT)
end
main()

# cost = 51.5
# solution = [25.0, 0.0, 20.0, 0.0, 36.0, 0.0]