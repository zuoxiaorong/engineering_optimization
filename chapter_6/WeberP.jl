using JuMP
using SCIP

function cal_accuracy_outer(epsilon)
    cita = acos(1 - 4*epsilon + 2*epsilon^2)
    n = Int(ceil(pi/(2*cita)))
    return n, cita
end

function cal_accuracy_inner(epsilon)
    cita = acos(1 - 4*epsilon + 2*epsilon^2)
    n = Int(ceil(pi/(2*cita)))
    a = [(sin(k*cita) - sin(k*cita - cita))/sin(cita) for k in 1:n]
    b = [(cos(k*cita - cita) - cos(k*cita))/sin(cita) for k in 1:n]
    return n, a, b
end

function data()
    XY = [
        75 54 80 23 47 36  9 53 47 44 73 58 38  3 70 99 50 73 42 25
        63 47 63  1 69 95 73 64 66 54 57 74 73 94 13 26 75 95 52 76
    ]
    N = 1:size(XY,2)
    X = XY[1,:]
    Y = XY[2,:]
    return N, X, Y
end

function lpmodel(model, N, X, Y, cita, n)
    @variable(model, x >= 0)
    @variable(model, y >= 0)
    @variable(model, dx[N])
    @variable(model, dy[N])
    @variable(model, d[N])

    @objective(model, Min, sum(d[i] for i in N))
    @constraint(model, [i in N], dx[i] >= x - X[i])
    @constraint(model, [i in N], dx[i] >= X[i] - x)
    @constraint(model, [i in N], dy[i] >= y - Y[i])
    @constraint(model, [i in N], dy[i] >= Y[i] - y)
    @constraint(model, [i in N, p in 0:n-1], d[i] >= dx[i]*cos(p*cita) + dy[i]*sin(p*cita))
    return model
end

function lpmodel_inner(model, N, X, Y, m, a, b)
    @variable(model, x >= 0)
    @variable(model, y >= 0)
    @variable(model, dx[N])
    @variable(model, dy[N])
    @variable(model, d[N])

    @objective(model, Min, sum(d[i] for i in N))
    @constraint(model, [i in N], dx[i] >= x - X[i])
    @constraint(model, [i in N], dx[i] >= X[i] - x)
    @constraint(model, [i in N], dy[i] >= y - Y[i])
    @constraint(model, [i in N], dy[i] >= Y[i] - y)
    @constraint(model, [i in N, p in 1:m], d[i] >= dx[i]*a[p] + dy[i]*b[p])
    return model
end

function startmodel_outer()
    model = Model(SCIP.Optimizer)
    epsilon = 0.001
    n, cita = cal_accuracy_outer(epsilon)
    N, X, Y = data()
    model = lpmodel(model, N, X, Y, cita, n)
    optimize!(model)
    println("obj = $(objective_value(model))")
    println("x = $(value.(model[:x]))")
    println("y = $(value.(model[:y]))")
end

function startmodel_inner()
    model = Model(SCIP.Optimizer)
    epsilon = 0.001
    n, a, b = cal_accuracy_inner(epsilon)
    N, X, Y = data()
    model = lpmodel_inner(model, N, X, Y, n, a, b)
    optimize!(model)
    println("obj = $(objective_value(model))")
    println("x = $(value.(model[:x]))")
    println("y = $(value.(model[:y]))")
end
startmodel_outer()

# 1 solution_outter
# obj = 556.0952666848525
# x = 49.09935061293565
# y = 65.53483210350612
# 2 solution_inner 
# obj = 556.6346140674315
# x = 49.32624035860518
# y = 65.5018679436858