# Uncertainty damaged in facility siting problem
# 面向不确定性损毁的设施选址问题
using JuMP
using SCIP
using COPT
using Distances
using Plots

function cal_distance(C,O)
    D = pairwise(Euclidean(), C, O)
    return D
end

function candidates_nodes()
    params = [
        15 62 0.1
        29 24 0.2
        39 39 0.15
        21 48 0.2
        67 62 0.5
        47 64 0.4
        89 49 0.2
        81 79 0.3
        52 75 0
        80 20 0.1
    ]
    XY = params[:,1:2]
    p = params[:,3]
    return XY', p
end

function demand_nodes()
    params = [
        57 69 66
        88 35 34
        68 14 68
        31 45 76
        94 8 85
        98 74 84 
        63 96 90
        84 62 48
        38 23 69
        66 61 39
        53 77 92
        85 45 92
        97 68 67
        68 6  84
        64 34 94
        63 54 93
        22 61 86
        12 38 38
        55 26 52
        37 10 42
        28 76 99
        10 2 87
        47 43 43
        17 19 12
        70 83 67
        28 84 93 
        3 52 15
        81 23 80
        58 60 69 
        9 87 69
    ]
    XY = params[:,1:2]
    a = params[:,3]
    return XY', a
end

function data()
    can_xy, p = candidates_nodes()      #保障点的坐标及损毁概率
    de_xy, a = demand_nodes()           #被保障点的坐标及需求
    N = 1:length(a)                     #被保障点的集合     
    K = 1:length(p)                     #保障设施的候选地址集合
    k = 4                               #保障设施的数量
    d = cal_distance(can_xy, de_xy)     #地址j与被保障点i之间的距离
    M = 999                             #一个大数
    return can_xy, de_xy, N, K, p, a, k, d, M
end

function lpmodel(model, N, K, p, a, k, d, M)
    @variable(model, x[K], Bin)             #0/1 变量，表示选择建造位置
    @variable(model, y[K,N], Bin)           #0/1 变量，表示保障服务的初始分配
    
    @objective(model, Min, sum(p[j]*a[i]*y[j,i] for i in N for j in K))
    @constraint(model, sum(x[j] for j in K) == k)
    @constraint(model, [j in K, i in N], y[j,i] <= x[j])
    @constraint(model, [i in N], sum(y[j,i] for j in K) == 1)
    @constraint(model, [i in N, j in K, j1 in K], d[j,i] <= d[j1,i] + M*(3 - y[j,i] + y[j1,i] - x[j] - x[j1]))
    return model
end

function plotresult(model, can_xy, de_xy)
    for p in 1:result_count(model)
        x_bin =  Int.(round.(collect(value.(model[:x]; result = p))))
        y_bin =  Int.(round.(collect(value.(model[:y]; result = p))))
        index = findall(x -> x > 0, x_bin)
        scatter(can_xy[1,:], can_xy[2,:], markersize= 5, markershape=:circle, markercolor=:purple, label="candidatesnodes", title="Solution$(p)")
        for i in index
            for j in findall(x -> x > 0, y_bin[i,:])
                sx, sy = can_xy[1,:][i], can_xy[2,:][i]
                ex, ey = de_xy[1,j], de_xy[2,j]
                plot!([sx, ex], [sy, ey], seriestype =:line, color=:black, linewidth=1, legend= false)
            end
        end
        scatter!(can_xy[1,:][index], can_xy[2,:][index], markersize= 6, markershape=:star, markercolor=:blue, label="chosen nodes")
        s = scatter!(de_xy[1,:], de_xy[2,:], markersize= 3, markershape=:cross, markercolor=:red, label="demandnodes")
        png(s, "newout$(p)")
    end
end

function startmodel()
    model = Model(SCIP.Optimizer)
    can_xy, de_xy, N, K, p, a, k, d, M = data()
    model = lpmodel(model, N, K, p, a, k, d, M)
    optimize!(model)
    plotresult(model, can_xy, de_xy)
end

@time startmodel()

# SCIP
# obj = 1.44850000000000e+02
# solved_time = 0.14
# time = 11.427696 seconds (23.00 M allocations: 1.367 GiB, 4.80% gc time, 11.01% compilation time)
# COPT
# obj = 144.850000000
# solved_time = 0.06
# time = 6.747458 seconds (22.44 M allocations: 1.342 GiB, 5.39% gc time, 12.57% compilation time)