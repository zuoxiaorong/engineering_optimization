using JuMP
using SCIP

function data()
    T = 1:6
    s = 10
    h = 0.5
    d = [10, 15, 8, 12, 20, 16]
    return T, s, h, d
end

function lpmodel(model, T, s, h, d)
    M = 999
    @variable(model, y[T], Bin)
    @variable(model, x[T], lower_bound = 0)
    @variable(model, I[T], lower_bound = 0)
    
    @constraint(model, I[1] == x[1] - d[1])
    @constraint(model, [i in T; i > 1], I[i] == I[i-1] + x[i] - d[i])
    @constraint(model, [i in T], x[i] - M*y[i] <= 0)

    @objective(model, Min, sum(h*I[i] + s*y[i] for i in T))
    return model
end

function startmodel()
    T, s, h, d = data()
    model = Model(SCIP.Optimizer)
    model = lpmodel(model, T, s, h, d)

    optimize!(model)
end
startmodel()
# solution
# obj = 5.15000000000000e+01
# y = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
# x = [25.0, 0.0, 20.0, 0.0, 36.0, 0.0]