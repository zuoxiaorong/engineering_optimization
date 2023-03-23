using JuMP
using SCIP

function x_squared(epsilon, x_min, x_max)
    miu = 1 + 2*epsilon + 2*sqrt(epsilon + epsilon*epsilon)
    eta = round(0.49999 + (log(x_max) - log(x_min))/log(miu))
    k = [(miu + 1)*(miu^(p - 1))*x_min for p in 1:Int(eta)]
    b = [-(miu^(2*p-1))*x_min^2 for p in 1:Int(eta)]
    return Int(eta), k, b
end

function x_derivative(epsilon, x_min, x_max)
    miu = 1 + 2*epsilon + 2*sqrt(epsilon + epsilon*epsilon)
    eta = round(0.49999 + (log(x_max) - log(x_min))/log(miu))
    k = [-1/((miu^(2*p - 1))*x_min^2) for p in 1:Int(eta)]
    b = [(miu + 1)/((miu^p)*x_min) for p in 1:Int(eta)]
    return Int(eta), k, b
end   

function x_root(epsilon, x_min, x_max)
    e = 1 - epsilon
    miu = (8 + e^4 - 8*e*e + (8 - 4*e*e)*(1 - e*e)^0.5)/e^4
    eta = round(0.49999 + (log(x_max) - log(x_min))/log(miu))
    k = [1/((sqrt(miu) + 1)*(sqrt(miu))^(p - 1)*sqrt(x_min)) for p in 1:Int(eta)]
    b = [1/((1 + miu^(-0.5))*((miu^(-0.5))^(p - 1))*x_min^(-0.5)) for p in 1:Int(eta)]
    return Int(eta), k, b
end 

function nonlpmodel(model, k1, k2, k3, b1, b2, b3, n1, n2, n3, x_min, x_max)
    @variable(model, x1, lower_bound = x_min, upper_bound = x_max)
    @variable(model, x2, lower_bound = x_min, upper_bound = x_max)
    @variable(model, x3, lower_bound = x_min, upper_bound = x_max)
    @variable(model, y1)
    @variable(model, y2)
    @variable(model, y3)

    @objective(model, Min, x1 - x2 - x3)
    @constraint(model, [p in 1:n1], y1 >= k1[p]*x1 + b1[p])
    @constraint(model, [p in 1:n2], y2 >= k2[p]*x2 + b2[p])
    @constraint(model, [p in 1:n3], y3 <= k3[p]*x3 + b3[p])
    @constraint(model, y1 == 0.5*x2 + 0.5*x3 + 4.5)
    @constraint(model, y2 == 2*x1 - x3 + 4)
    @constraint(model, y3 == 1 - 0.5*x1 + 0.5*x2)
    return model
end

function startmodel()
    model = Model(SCIP.Optimizer)
    epsilon = 0.001
    x_min = 0.1
    x_max = 10

    n1, k1, b1 = x_squared(epsilon, x_min, x_max)
    n2, k2, b2 = x_derivative(epsilon, x_min, x_max)
    n3, k3, b3 = x_root(epsilon, x_min, x_max)
   
    model = nonlpmodel(model, k1, k2, k3, b1, b2, b3, n1, n2, n3, x_min, x_max)
    optimize!(model)
    println("obj = $(objective_value(model))")
    println("x1 = $(value.(model[:x1]))")
    println("x2 = $(value.(model[:x2]))")
    println("x3 = $(value.(model[:x3]))")
end
startmodel()
# solution 
# obj = -14.31960347011567
# x1 = 3.0676932249927713
# x2 = 7.38729669510844
# x3 = 10.0