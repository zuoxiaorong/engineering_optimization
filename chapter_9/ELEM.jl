# Optimization problem of equipment life extension maintenance based on reliability
# 基于可靠性的设备延寿维修优化问题

using JuMP
using SCIP

function data()
    T = 1:6                                         # 设备延寿后的使用周期集合
    h = [1000, 1000, 1000, 1500, 1500, 1500]        # 设备在期间t内的工作量
    D = [3, 3, 3, 5, 5, 5]                          # 设备在期间t内的日均收益
    F = [10, 10, 10, 15, 15, 20]                    # 设备发生故障后产生的损失，如修复费用、停机摄失等
    λ1 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]        # 以设备当前现状(无任何维修活动)在期间t的故障率
    p1 = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5]             # 以设备当前现状(无任何维修活动)在期间t的运行成本
    R = 1:4                                         # 设备的计划性维修方案的集合(包括:现状维修,增强维修,设备大修,更换核心机)
    C = [5, 10, 20, 15]                             # 方案工的成本
    d = [1, 1.5, 2, 2]                              # 执行方案工的停机天数
    λ = [                                           # 设备在期间t执行方案工后，在后续期间的故障率
        0.5 0.5 0.6 0.6 0.8 0.8
        0.8 0.8 0.9 0.9 0.1 0.1
        0.2 0.2 0.8 1   1.2 1.4
        0.3 0.4 0.5 0.7 0.9 1.0
    ]
    p = [                                           # 设备在期间t执行方案工后，在后续期间的运行成本
        2 2 2 2 2 2
        5 5 5 5 5 5
        3 3 3 3 3 3
        4 4 4 4 4 4
    ]
    M = 9999
    return T, h, D, F, λ1, p1, R, C, d, λ, p, M 
end

function lpmodel(model, T, h, D, F, λ1, p1, R, C, d, λ, p, M )
    @variable(model, x[R,T], Bin)                   # 0/1 变量,表示是否在期间t执行了方案r
    @variable(model, y[R,T], Bin)                   # 0/1 变量，表示期间t的有效方案是否为方案r
    @variable(model, f[T] >= 0)                     # 非负连续变量，表示设备在期间t内发生的期望故障次数
    @variable(model, P[T] >= 0)                     # 非负连续变量，表示设备在期间t内的运行成本
    @variable(model, z[T] >= 0)                     # 非负连续变量，表示设备在期间t内的维修停机天数

    @objective(model, Min, sum(C[r]*x[r,t] for r in R for t in T) + sum(D[t]*z[t] + F[t]*f[t] + P[t] for t in T))
    @constraint(model, [t in T], sum(x[r,t] for r in R) <= 1)
    @constraint(model, [r in R, t in T; t > 1], 2*y[r,t] >= y[r,t-1] + x[r,t] - sum(x[rl,t] for rl in R if rl != r))
    @constraint(model, [r in R, t in T], y[r,t] <= 1 - sum(x[r1,t] for r1 in R if r1 != r))
    @constraint(model, [t in T], sum(y[r,t] for r in R) <= 1)
    @constraint(model, [r in R, t in T], y[r,t] <= sum(x[r,t1] for t1 in T if t1 != t))
    @constraint(model, [r in R, t in T], y[r,t] >= x[r,t])
    @constraint(model, [r in R, t in T, τ in T; τ <= t], f[t] >= h[t]*λ[r,t] - M*(2-x[r,τ] - y[r,t]))
    @constraint(model, [r in R, t in T, τ in T; τ <= t], f[t] <= h[t]*λ[r,t] + M*(2 - x[r,τ] - y[r,t]))
    @constraint(model, [t in T], f[t] >= h[t]*λ1[t] - M*sum(x[r,t1] for r in R for t1 in T if t1 <= t))
    @constraint(model, [t in T], f[t] <= h[t]*λ1[t] + M*sum(x[r,t1] for r in R for t1 in T if t1 <= t))
    @constraint(model, [r in R, t in T, τ in T; τ <= t], P[t] >= p[r,t] - M*(2 - x[r,τ] - y[r,t]))
    @constraint(model, [r in R, t in T, τ in T; τ <= t], P[t] <= p[r,t] + M*(2 - x[r,τ] - y[r,t]))
    @constraint(model, [t in T], P[t] >= p1[t] - M*sum(x[r,t1] for r in R for t1 in T if t1 <= t))
    @constraint(model, [t in T], P[t] <= p1[t] + M*sum(x[r,t1] for r in R for t1 in T if t1 <= t))
    @constraint(model, [t in T], z[t] >= sum(d[r]*x[r,t] for r in R))
    return model
end

function startmodel()
    model = Model(SCIP.Optimizer)
    T, h, D, F, λ1, p1, R, C, d, λ, p, M  = data()
    model = lpmodel(model, T, h, D, F, λ1, p1, R, C, d, λ, p, M )
    optimize!(model)
end
startmodel()
