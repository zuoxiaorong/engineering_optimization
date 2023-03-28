using JuMP
using SCIP

function data()
    N = 1:10
    params_N = [
        1.1 4
        1.2 1
        1.4 4
        0.9 2
        1.5 2
        0.8 5
        1   4
        0.9 2
        0.8 7
        3   2
    ]
    M = 9999
    return N, params_N[:,1], params_N[:,2], M
end

function data_new()
    N = 1:10
    params_N = [
        1.1 4
        1.2 1
        1.4 4
        0.9 2
        1.5 2
        0.8 5
        1   4
        0.9 2
        0.8 7
        3   2
    ]
    M = 9999
    d = [10, 5, 10, 3, 5, 10, 7, 8, 9, 15]
    return N, params_N[:,1], params_N[:,2], M, d
end

function mlpmodel(model, N, w, p, M)
    @variable(model, s[N,N], Bin)                   #s[i,j] = 1 代表i订单排在j订单前被处理
    @variable(model, c[N] >= 0)                     #订单的完成时间

    @objective(model, Min, sum(c[i]*w[i] for i in N))
    @constraint(model, [i in N, j in N; i != j], s[i,j] + s[j,i] == 1)
    @constraint(model, [i in N], c[i] >= p[i])
    @constraint(model, [i in N, j in N; i != j], c[j] - c[i] >= p[j] - M*(1 - s[i,j]))
    return model
end

function mlpmodel_new(model, N, w, p, M, d)             #带惩罚时间的
    @variable(model, s[N,N], Bin)                   #s[i,j] = 1 代表i订单排在j订单前被处理
    @variable(model, c[N] >= 0)                     #订单的完成时间
    @variable(model, t[N] >= 0)                     #订单延期交货的时间

    @objective(model, Min, sum(t[i]*w[i] for i in N))
    @constraint(model, [i in N, j in N; i != j], s[i,j] + s[j,i] == 1)
    @constraint(model, [i in N], c[i] >= p[i])
    @constraint(model, [i in N, j in N; i != j], c[j] - c[i] >= p[j] - M*(1 - s[i,j]))
    @constraint(model, [i in N], c[i] - d[i] <= t[i])
    # @constraint(model, [i in N], c[i] <= d[i])
    return model
end

function startmodel()
    model = Model(SCIP.Optimizer)
    N, w, p, M = data()
    model = mlpmodel(model, N, w, p, M)
    optimize!(model)
    s = value.(model[:s])
    order = sortperm([sum(s[i,j] for j in N) for i in N]; rev=true)    
end

function startmodel_new()
    model = Model(SCIP.Optimizer)
    N, w, p, M, d = data_new()
    model = mlpmodel_new(model, N, w, p, M, d)
    optimize!(model)
    s = value.(model[:s])
    order = sortperm([sum(s[i,j] for j in N) for i in N]; rev=true)    
end
startmodel_new()

# solution
# obj = 1.36599999999995e+02
# order = [10, 2, 5, 4, 8, 3, 1, 7, 6, 9]
# solution_new 延迟惩罚
# obj = +5.51000000000000e+01
# order = [2, 4, 5, 8, 3, 10, 1, 7, 6, 9]
# solution_new 必须早于交货时间完成
# infeasible 交货时间的设置太紧，没办法完成
