using JuMP
using SCIP

function data()
    N = 1:10
    D = [45, 35, 40, 55, 20, 50, 30, 40, 30, 50]
    d1 = [3, 6, 8, 10, 1, 48, 28, 25, 15, 5]
    d2 = [46, 34, 39, 54, 21, 5, 8, 25, 20, 48]
    L = 100   
    return N, D, d1, d2, L
end

function lpmodel(model, N, D, d1, d2, L)
    @variable(model, y[N], Bin)
    @variable(model, r[N] >= 0)

    @objective(model, Min, sum((1 - y[i])*D[i] + y[i]*(d1[i] + d2[i]) for i in N))
    @constraint(model, L >= d1[1]*y[1])
    @constraint(model, L >= r[1] + d2[1]*y[1])
    @constraint(model, L - r[1] >= D[1]*(1 - y[1]))
    @constraint(model, [i in N; i >= 2], r[i-1] >= d1[i] + L*(y[i] - 1))
    @constraint(model, [i in N; i >= 2], L - r[i] >= d2[i]*y[i])
    @constraint(model, [i in N; i >= 2], r[i-1] >= r[i] + D[i]*(1 - y[i]) - L*y[i])
    return model
end

function startmodel()
    N,  D, d1, d2, L = data()
    model = Model(SCIP.Optimizer)
    model = lpmodel(model, N,  D, d1, d2, L)

    optimize!(model)
end
startmodel()

# solution
# obj = +4.12000000000000e+02
# y = [-0.0, -0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0]