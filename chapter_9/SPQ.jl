# Optimization problems based on the spare parts quantity based on an exponential distribution
# 基于指数分布的备件数量优化问题

using JuMP
using SCIP
using COPT
using Plots

function data(x)
    θ = 1000                            # 表示平均无故障间隔时间(MTBF), hour
    T = 24000                           # 年平均工作量
    p = 0.5                             # 期初单位价格
    u = 1.5                             # 单位保管费用
    λ = T/θ                             # 备件发生故障的次数, 即 λ = T/θ
    o = 3.5                             # 发生故障后的故障损失, 包括备件订购成本和停机损失成本
    exp_y = sum(k*(λ^(x + k))*(ℯ^(-λ))/factorial(big(x + k)) for k in 1:50)                     #这里计算不是特别准确
    return p, u, o, exp_y
end

function lpmodel(model, p, u, o, exp_y, x)
    @variable(model, y >= 0)            # 备件缺货的期望次数

    @objective(model, Min, x*(u + p) + y*o)
    
    @constraint(model, y >= exp_y)
    return model
end

function plotresult(num, order_num, cost)
    order_fig = plot(num, order_num, color = :orange, line = (:dot, 3), marker = (:hexagon, 5, 0.8), label="Y Expection")
    cost_fig = plot(num, cost, color = :blue, line = (:dot, 3), marker = (:d, 5, 0.8), label="Cost Expection")
    F = plot(order_fig, cost_fig, layout=(2,1))
    png(F, "SPQ_out")
end

function startmodel()
    num = []
    cost = []
    order_num = []
    for x in 0:50
        model = Model(COPT. Optimizer) 
        p, u, o, exp_y = data(x)
        model = lpmodel(model, p, u, o, exp_y, x)
        optimize!(model)
        push!(num, x)
        push!(order_num, value.(model[:y]))
        push!(cost, objective_value(model))
    end
    plotresult(num, order_num, cost)
end
startmodel()
# solution
# x = 23
# cost = 54.661752190372