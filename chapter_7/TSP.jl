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
        0 25.1 13.9 3.5  5.2  7.5  35.6 37.8 1.6  11.4 32.6
        0 36.3 29.2 10.4 38.5 21.7 6.3  2.4  37.3 36.3 12.5
    ]
    return XY
end

function data(XY)
    NODE = 1:size(XY,2)
    n = length(NODE) - 1
    D = cal_distance(XY, XY)
    return NODE, n, D
end

function milpmodel(model, NODE, n, D)
    @variable(model, x[NODE,NODE], Bin)
    @variable(model, u[NODE] >= 0)

    @objective(model, Min, sum(x[i,j]*D[i,j] for i in NODE for j in NODE if i != j))
    @constraint(model, [i in NODE], sum(x[i,j] for j in NODE if i != j) == 1)
    @constraint(model, [i in NODE], sum(x[j,i] for j in NODE if i != j) == 1)
    @constraint(model, [i in NODE, j in NODE; i != j && j > 1], u[j] - u[i] >= 1 - n*(1 - x[i,j]))
    @constraint(model, [i in NODE; i > 1], u[i] <= n)
    return model
end

function plotresult(u, XY)
    x = XY[1,:]
    y = XY[2,:]
    order = sortperm(u)
    x_order = push!(x[order], x[1])
    y_order = push!(y[order], y[1])
    de_x = [x_order[i+1] - x_order[i] for i in 1:length(x)]
    de_y = [y_order[i+1] - y_order[i] for i in 1:length(y)]
    quiver(x_order, y_order, quiver=(de_x,de_y))
    s = plot!(x, y, seriestype=:scatter, series_annotation = text.(collect(0:10), :left, 11), legend=false)
    png(s, "tspout")
end

function startmodel()
    model = Model(SCIP.Optimizer)
    XY = allnodes()
    NODE, n, D = data(XY)
    model = milpmodel(model, NODE, n, D)
    optimize!(model)
    u_v = [value.(model[:u])[i] for i in NODE]
    plotresult(u_v, XY)
end
startmodel()
