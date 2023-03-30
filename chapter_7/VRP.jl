using JuMP
using SCIP
using Distances
using Plots

function cal_distance(C,O)
    D = pairwise(Euclidean(), C, O)
    return D
end

function allnodes()
    XY = [
        50 61 19 92 96 11 71 55 26 35 41 67 58 94 36 5 76 45 17 79 84
        50 59 20 25 53 17 13 33 81 64 75 49 8 65 48 31 60 39 96 59 76
    ]
    return XY
end

function data(XY)
    NODE = 1:size(XY,2)
    n = length(NODE) - 1
    D = cal_distance(XY, XY)
    m = 3
    return NODE, n, D, m
end

function milpmodel(model, NODE, n, D, m)
    @variable(model, x[NODE,NODE], Bin)
    @variable(model, u[NODE] >= 0)

    @objective(model, Min, sum(x[i,j]*D[i,j] for i in NODE for j in NODE if i != j))
    @constraint(model, sum(x[1,j] for j in NODE if j > 1) == m)
    @constraint(model, [i in NODE; i > 1], sum(x[i,j] for j in NODE if i != j) == 1)
    @constraint(model, [i in NODE; i > 1], sum(x[j,i] for j in NODE if i != j) == 1)
    @constraint(model, [i in NODE, j in NODE; i != j && j > 1], u[j] - u[i] >= 1 - n*(1 - x[i,j]))
    # @constraint(model, [i in NODE, j in NODE; i != j && j > 1], u[j] - u[i] <= 1 + n*(1 - x[i,j]))
    @constraints(model, begin
        [i in NODE; i > 1], u[i] <= ceil(n/m)
        u[1] == 0
    end)
    return model
end

function plotresult(path, XY)
    X = XY[1,:]
    Y = XY[2,:]
    x_order = map(x -> X[x], path)
    y_order = map(x -> Y[x], path)
    de_x = [[x_order[i][j+1] - x_order[i][j] for j in 1:length(x_order[i])-1] for i in 1:length(path)]
    map(x -> pop!(x_order[x]), 1:length(path))
    # x_order = [x_order[i][1:end-1] for i in 1:length(path)]
    de_y = [[y_order[i][j+1] - y_order[i][j] for j in 1:length(y_order[i])-1] for i in 1:length(path)]
    map(x -> pop!(y_order[x]), 1:length(path))
    map(x -> quiver!(x_order[x], y_order[x], quiver=(de_x[x],de_y[x])), 1:length(path))
    s = plot!(X, Y, seriestype=:scatter, series_annotation = text.(collect(0:length(X) - 1), :left, 11), legend=false, aspect_ratio = 1)
    png(s, "vrpout")
end

function startmodel()
    model = Model(SCIP.Optimizer)
    XY = allnodes()
    NODE, n, D, m = data(XY)
    model = milpmodel(model, NODE, n, D, m)
    optimize!(model)
    z = collect(value.(model[:x]))[NODE, NODE]
    firstnodes = findall(x -> x > 0.8, z[1,:])
    startvehicle = 1
    path = []
    for h in firstnodes
        tour = []
        cur_node = h
        push!(tour, 1, cur_node)
        while cur_node != 1
            ii = cur_node
            cur_node = findfirst(x->x>0.5, z[ii,:])
            push!(tour, cur_node)
        end
        push!(path, tour)
        startvehicle += 1
    end
    plotresult(path, XY)
end
startmodel()

# solution 
# Vehicel: 1
# 1  2  21  14  5  20  17  12  1    
# Vehicel: 2
# 1  8  4  7  13  3  6  16  1    
# Vehicel: 3
# 1  18  15  10  9  19  11  1  
