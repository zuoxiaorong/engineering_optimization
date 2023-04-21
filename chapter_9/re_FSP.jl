# facility siting problem with system elastic recovery
# 面向系统弹性恢复的设施选址优化问题
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
    return XY'
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
    can_xy = candidates_nodes()         #保障点的坐标及损毁概率
    de_xy, a = demand_nodes()           #被保障点的坐标及需求
    N = 1:length(a)                     #被保障点的集合     
    K = 1:size(can_xy,2)                #保障设施的候选地址集合
    k = 4                               #保障设施的数量
    d = cal_distance(can_xy, de_xy)     #地址j与被保障点i之间的距离
    M = 999                             #一个大数
    R = 1:3                             #保障点发生损毁的数量状况集合
    p = [0.3, 0.1, 0.05]                #R集合中每个情况发生的概率
    return can_xy, de_xy, N, K, p, a, k, d, M, R
end

function lpmodel(model, N, K, p, a, k, d, M, R)
    @variable(model, x[K], Bin)             #0/1 变量，表示选择建造位置
    @variable(model, y[K,N], Bin)           #0/1 变量，表示保障服务的初始分配
    @variable(model, z[R,K], Bin)           #0/1 变量，表示在场景r下j是否发生故障或损毁
    @variable(model, y1[R,K,N], Bin)        #0/1 变量，表示在场景r下是否由j向i提供服务

    @objective(model, Min, (1 - sum(p[r] for r in R))*sum(a[i]*d[j,i]*y[j,i] for i in N for j in K) + sum(p[r]*a[i]*d[j,i]*y1[r,j,i] for r in R for j in K for i in N))
    @constraint(model, sum(x[j] for j in K) == k)                                                                                           #建造数量约束
    @constraint(model, [i in N], sum(y[j,i] for j in K) == 1)                                                                               #每一个需求点都有一个点为它提供服务
    @constraint(model, [j in K, i in N], y[j,i] <= x[j])                                                                                    #只有i点有建设保障点才可以为其他点提供服务
    @constraint(model, [i in N, j in K, j1 in K], d[j,i] <= d[j1,i] + M*(3 - y[j,i] + y[j1,i] - x[j] - x[j1]))                              #需求点选择建造的保障点中最近的一个提供服务
    @constraint(model, [r in R], sum(z[r,j] for j in K) == r)                                                                               #场景r下发生损毁的数量约束
    @constraint(model, [r in R, j in K], z[r,j] <= x[j])                                                                                    #只有i点有建设保障点才可以被损毁
    @constraint(model, [r in R, j in K, j1 in K], sum(a[i]*y[j,i] for i in N) >= sum(a[i]*y[j1,i] for i in N) - M*(1 - z[r,j] + z[r,j1]))   #需求最高的几个点发生故障
    @constraint(model, [r in R, i in N], sum(y1[r,j,i] for j in K) == 1)                                                                    #损毁后每一个需求点仍然被保障
    @constraint(model, [r in R, i in N, j in K], y1[r,j,i] <= x[j] - z[r,j])                                                                #损毁后，重新分配保障关系后的约束
    return model
end

function plotresult(model, can_xy, de_xy, R)
    P = plot()
    x_bin =  Int.(round.(collect(value.(model[:x]))))
    y_bin =  Int.(round.(collect(value.(model[:y]))))
    index = findall(x -> x > 0, x_bin)
    scatter(can_xy[1,:], can_xy[2,:], markersize= 5, markershape=:circle, markercolor=:purple, label="candidatesnodes", title="R = 0")
    for i in index
        for j in findall(x -> x > 0, y_bin[i,:])
            sx, sy = can_xy[1,:][i], can_xy[2,:][i]
            ex, ey = de_xy[1,j], de_xy[2,j]
            plot!([sx, ex], [sy, ey], seriestype =:line, color=:black, linewidth=1, legend= false)
        end
    end
    scatter!(can_xy[1,:][index], can_xy[2,:][index], markersize= 6, markershape=:star, markercolor=:blue, label="chosen nodes")
    s = scatter!(de_xy[1,:], de_xy[2,:], markersize= 3, markershape=:cross, markercolor=:red, label="demandnodes")
    figures = [s]
    for r in R
        z_bin = Int.(round.(collect(value.(model[:z]))[r,:]))
        new_x_bin = x_bin .- z_bin
        new_y_bin =  Int.(round.(collect(value.(model[:y1]))[r,:,:]))
        index = findall(x -> x > 0, new_x_bin)
        scatter(can_xy[1,:], can_xy[2,:], markersize= 5, markershape=:circle, markercolor=:purple, label="candidatesnodes", title="R = $(r)")
        for i in index
            for j in findall(x -> x > 0, new_y_bin[i,:])
                sx, sy = can_xy[1,:][i], can_xy[2,:][i]
                ex, ey = de_xy[1,j], de_xy[2,j]
                plot!([sx, ex], [sy, ey], seriestype =:line, color=:black, linewidth=1, legend= false)
            end
        end
        scatter!(can_xy[1,:][index], can_xy[2,:][index], markersize= 6, markershape=:star, markercolor=:blue, label="chosen nodes")
        n = scatter!(de_xy[1,:], de_xy[2,:], markersize= 3, markershape=:cross, markercolor=:red, label="demandnodes")
        push!(figures, n)
    end
    S = plot(figures..., layout=(length(R)+1))
    png(S, "re_FSP_out")
end

function startmodel()
    model = Model(COPT.Optimizer)
    can_xy, de_xy, N, K, p, a, k, d, M, R = data()
    model = lpmodel(model, N, K, p, a, k, d, M, R)
    optimize!(model)
    plotresult(model, can_xy, de_xy, R)
end
startmodel()
# solution
# obj = +4.80891687496860e+04