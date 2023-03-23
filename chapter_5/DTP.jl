using JuMP
using SCIP

function data()
   V = 1:18
   arc = [
        (1, 18)  1.8
        (1, 10)  1.2
        (1, 4)   2.6
        (1, 5)   4.8
        (1, 11)  3.3
        (1, 12)  3.9
        (2, 7)   4.5
        (2, 9)   1.6
        (2, 12)  1
        (2, 18)  3.1
        (3, 6)   2.3
        (3, 16)  1.6
        (3, 17)  2.9
        (4, 11)  4.8
        (4, 8)   4.1
        (4, 5)   2.5
        (5, 8)   3
        (5, 12)  2.9
        (5, 14)  3.2
        (6, 7)   4.6
        (6, 12)   2.6
        (6, 17)   3.2
        (7, 9)    1.7
        (7, 12)   3
        (8, 13)   4.7
        (8, 14)   2.2
        (8, 15)   2.5
        (9, 10)   4.3
        (9, 18)   1.8
        (10, 11)  1.5
        (10, 18)  4.2
        (12, 17)  4.1
        (12, 18)  2.9
        (13, 15)  2.7
        (14, 15)  3.5
        (14, 16)  2.4
        (14, 17)  2.8
        (15, 16)  3
        (16, 17)  3
   ] 
   E = vcat([(i,j) for i in V for j in V if (i,j) in arc[:,1]], [(i,j) for i in V for j in V if (j,i) in arc[:,1]])
   c = [arc[:,2][findfirst(x -> (x == E[k] || x == reverse(E[k])), arc[:,1])] for k in 1:length(E)]
   s = 3
   t = 11
   a = 1024
   return V, E, c, s, t, a
end

function lpmodel(model, V, E, c, s, t, a)
    @variable(model, x[i in V, j in V; (i,j) in E], lower_bound = 0)                #弧上数据量
    @variable(model, T, lower_bound = 0)                                            #最大传输时间

    @objective(model, Min, T)
    @constraint(model, [i in V; (i, s) in E], x[i,s] == 0)
    @constraint(model, [i in V; (t, i) in E], x[t,i] == 0)
    @constraint(model, sum(x[s,j] for j in V if (s,j) in E) == a)
    @constraint(model, sum(x[i,t] for i in V if (i,t) in E) == a)
    @constraint(model, [i in V; i != s && i != t], sum(x[j,i] for j in V if (j,i) in E && j != t) == sum(x[i,j] for j in V if (i,j) in E && j != s))
    @constraint(model, [i in V, j in V, k in 1:length(c); (i,j) == E[k]], T >= x[i,j]/c[k])
    return model
end
    
function startmodel()
    model = Model(SCIP.Optimizer)

    V, E, c, s, t, a = data()
    model = lpmodel(model, V, E, c, s, t, a)
    println(model)
    optimize!(model)
    dt = value.(model[:x])
    display([("[$(i),$(j)]", dt[i,j]) for i in V for j in V if (i,j) in E && dt[i,j] > 0])
end
    
startmodel()
