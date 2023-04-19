using JuMP
using SCIP

function data()
    N = 1:9
    r = [0, 7, 3, 5, 4, 1, 3, 2, 4]
    d = [2, 1, 1, 3, 1, 2, 2, 2]
    L = 5
    return N, r, d, L
end

function lpmodel(model, N, r, d, L)
    @variable(model, x[N], Bin)
    @variable(model, y[N] >= 0)

    @objective(model, Min, sum(x[i]*r[i] for i in N))
    @constraint(model, y[1] == L)
    @constraint(model, [i in N; i < length(N)], y[i] - y[i+1] >= d[i]*(1 - x[i]) - x[i]*L)
    @constraint(model, [i in N; i < length(N)], y[i] - y[i+1] <= d[i]*(1 - x[i]) + x[i]*L)
    @constraint(model, [i in N; i < length(N)], y[i+1] >= (L - d[i])*x[i] - L*(1 - x[i]))
    @constraint(model, [i in N; i < length(N)], y[i+1] <= (L - d[i])*x[i] + L*(1 - x[i]))
    return model
end

function startmodel()
    N, r, d, L = data()
    model = Model(SCIP.Optimizer)
    model = lpmodel(model, N, r, d, L)

    optimize!(model)
end
startmodel()
# solution 
# obj = +6.00000000000000e+00
# plan = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0]