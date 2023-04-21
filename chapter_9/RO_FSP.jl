# facility siting problem with robust optimization
# 不确定攻击下设施选址健壮性优化问题

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
        15 62 350
        29 24 400
        39 39 450
        21 48 500
        67 62 250
        47 64 350
        89 49 250
        81 79 300
        52 75 350
        80 20 450
    ]
    XY = params[:,1:2]
    C = params[:,3]
    return XY', C
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
    can_xy, C = candidates_nodes()      #保障点的坐标及容量上限
    de_xy, a = demand_nodes()           #被保障点的坐标及需求
    N = 1:length(a)                     #被保障点的集合     
    K = 1:size(can_xy,2)                #保障设施的候选地址集合
    k = 6                               #保障设施的数量
    d = cal_distance(can_xy, de_xy)     #地址j与被保障点i之间的距离
    r = 3                               #被攻击的设施数量
    M = 9999                            #一个大数
    return can_xy, de_xy, N, K, a, d, k, r, C, M
end

function lpmodel(model, N, K, a, d, k, r, C, M)
    @variable(model, x[K], Bin)             #0/1 变量，表示选择建造位置
    @variable(model, y[K,N], Bin)           #0/1 变量，表示保障服务的初始分配
    @variable(model, s[K], Bin)             #0/1 变量，表示设施是否被攻击
    @variable(model, y1[K,N], Bin)          #0/1 变量，表示被攻击后中断的服务关系
    
    @objective(model, Min, sum(a[i]*d[j,i]*y1[j,i] for j in K for i in N))

    @constraint(model, sum(x[j] for j in K) == k)
    @constraint(model, [i in N], sum(y[j,i] for j in K) == 1)
    @constraint(model, [j in K, i in N], y[j,i] <= x[j])
    @constraint(model, [j in K], sum(a[i]*y[j,i] for i in N) <= C[j]*x[j])
    @constraint(model, sum(s[j] for j in K) == r)
    @constraint(model, [j in K], s[j] <= x[j])
    @constraint(model, [i in N, j in K], y1[j,i] <= y[j,i])
    @constraint(model, [i in N, j in K], y1[j,i] <= s[j])
    @constraint(model, [i in N, j in K], y1[j,i] >= 1 - (2 - y[j,i] - s[j]))
    @constraint(model, [j in K, j1 in K], sum(a[i]*d[j,i]*y[j,i] for i in N) >= sum(a[i]*d[j1,i]*y[j1,i] for i in N) - M*(3 - s[j] + s[j1] - x[j] - x[j1]))   
    return model 
end

function plotresult(model, can_xy, de_xy)
    x_bin = Int.(round.(collect(value.(model[:x]))))
    y_bin = Int.(round.(collect(value.(model[:y]))))
    s_bin = Int.(round.(collect(value.(model[:s]))))
    index = findall(x -> x > 0, x_bin)
    index_b = findall(x -> x > 0, s_bin)
    scatter(can_xy[1,:], can_xy[2,:], markersize= 8, markershape=:circle, markercolor=:pink, label="candidatesnodes", title="Solution")
    for i in index
        for j in findall(x -> x > 0, y_bin[i,:])
            sx, sy = can_xy[1,:][i], can_xy[2,:][i]
            ex, ey = de_xy[1,j], de_xy[2,j]
            plot!([sx, ex], [sy, ey], seriestype =:line, color=:black, linewidth=1, legend= false)
        end
    end
    scatter!(can_xy[1,:][index], can_xy[2,:][index], markersize= 8, markershape=:star, markercolor=:red, label="chosen nodes")
    scatter!(can_xy[1,:][index_b], can_xy[2,:][index_b], markersize= 3, markershape=:circle, markercolor=:blue, label="broken nodes")
    s = scatter!(de_xy[1,:], de_xy[2,:], markersize= 3, markershape=:cross, markercolor=:red, label="demandnodes")
    png(s, "RO_FSP_out")
end

function startmodel()
    model = Model(COPT. Optimizer)
    can_xy, de_xy, N, K, a, d, k, r, C, M = data()
    model = lpmodel(model, N, K, a, d, k, r, C, M)
    optimize!(model)
    plotresult(model, can_xy , de_xy)
end
startmodel()
# solution
# obj = +1.97104944595954e+04