using JuMP
using SCIP
using Distances

function xy()
    params = [
        2.0   83.3  1 
        8.8   110.0 1
        -9.5  58.8  1
        48.6  9.0   3
        70.5  20.0  3
        93.2  27.3  3
        95.5  8.9   3
        93.5  68.0  3
        78.0  109.5 2
        20.6  86.9  2
        69.5  85.7  2
        23.0  30.5  2
        61.5  80.5  2
        37.0  94.8  2
        47.4  56.4  2
        17.2  66.7  2
        62.0  95.0  2
        20.8  68.8  2
        3.4   25.5  2
        36.9  28.4  2
        21.2  39.9  2
        0.5   40.0  2
        22.0  97.3  2
        46.1  83.2  2
        55.1  68.8  2
        10.1  51.9  2
        38.3  61.5  2
        38.0  72.8  2
        35.0  81.5  2
        6.5   9.0   2
    ]
    return params[:,1:2]', params[:,3]
end

function cal_distance(C,O)
    D = pairwise(Euclidean(), C, O)
    return D
end

function data()
    Node, a = xy()
    D = cal_distance(Node, Node)
    T = findall(x -> x == 3, a)                 #保障目标集合
    S = findall(x -> x == 1, a)                 #机场集合
    K = findall(x -> x == 2, a)                 #待选中继点集合
    N = 1:length(a)                             #union(T,S,K)
    H = 1:2                                     #机型集合
    Link = [
        (1, 2, 4),
        (1, 1, 4),
        (2, 3, 4),
        (1, 2, 5),
        (2, 1, 5),
        (1, 1, 5),
        (1, 3, 6),
        (2, 1, 6),
        (1, 3, 7),
        (2, 1, 7),
        (2, 1, 8),
        (1, 3, 8),
        (2, 2, 8)
    ]
    w = [5, 3, 3, 6, 3, 5, 4, 2, 6, 2, 6, 2, 4]
    r = [185, 150] 
    p = [0.7, 0.9, 0.8, 0.6, 0.5]
    e = 4 
    M = 9999
    return T, S, K, N, H, Link, w, r, p, D, a, e, M
end

function milpmodel(model, T, S, K, N, H, Link, w, r, p, D, a, e, M)
    @variable(model, z[K], Bin)                         #是否选择该点建立中继点
    @variable(model, y[N], Bin)                         #是否允许飞机起降
    @variable(model, x[h in H, s in S, t in T, i in N, j in N; (h, s, t) in Link && i != j], Bin)                     #路径的选择变量

    @objective(model, Min, sum(x[h,s,t,i,j]*D[i,j]*w[k]*p[findfirst(x -> x == t, T)] for h in H for s in S for t in T for i in N for j in N for k in 1:length(w) if (h, s, t) == Link[k] && i != j))
    @constraint(model, [h in H, s in S, t in T; (h, s, t) in Link], sum(x[h,s,t,s,j] for j in N if j != s) == 1)
    @constraint(model, [h in H, s in S, t in T, i in N; (h, s, t) in Link && i != s && i != t], sum(x[h,s,t,j,i] for j in N if i != j) == sum(x[h,s,t,i,j] for j in N if i != j))
    @constraint(model, [h in H, s in S, t in T; (h, s, t) in Link], sum(x[h,s,t,j,t] for j in N if j != t) == 1)
    @constraint(model, [i in N; a[i] == 1], y[i] == 1)
    @constraint(model, sum(y[i] for i in N if a[i] == 3) == 0)
    @constraint(model, [i in K], y[i] == z[i])
    @constraint(model, sum(y[i] for i in K) <= e)
    @constraint(model, [h in H, s in S, t in T, i in N, j in N; (h, s, t) in Link && i != j], x[h,s,t,i,j] <= y[i])
    @constraint(model, [h in H, s in S, t in T, i in N, j in N; (h, s, t) in Link && i != j && a[j] != 3], x[h,s,t,i,j]*D[i,j] <= r[h])
    @constraint(model, [h in H, s in S, t in T, i in N, j in N; (h, s, t) in Link && i != j && a[j] == 3], x[h,s,t,i,j]*D[i,j] <= 0.5*r[h])
    return model
end

function  startmodel()
    optimizer = SCIP.Optimizer(display_verblevel=4 , limits_gap=0.00, parallel_maxnthreads=12)
    model = JuMP.direct_model(optimizer)
    set_time_limit_sec(model,100)

    T, S, K, N, H, Link, w, r, p, D, a, e, M = data()
    model = milpmodel(model, T, S, K, N, H, Link, w, r, p, D, a, e, M)
    optimize!(model)
    c = value.(model[:z])
    candidate = [i for i in K if c[i] > 0]
    path_bin = value.(model[:x])
    path = unique([("Path$(k):($(h),$(s),$(t))", s, ifelse(isempty(setdiff(union(s,i,j,t), union(s,t))), nothing, sum(setdiff(union(i,j,s,t), union(s,t)))), t) for k in 1:length(w) for h in H for s in S for t in T for i in N for j in N if (h, s, t) == Link[k] && i != j && path_bin[h,s,t,i,j] > 0])
    display(path)
end
startmodel()

# solution
# obj = +3.64987170724828e+03
# candidate = [13, 15, 18, 21]
# path = [
#     ("Path1:(1,2,4)", 2, 18, 4),
#     ("Path2:(1,1,4)", 1, missing, 4),
#     ("Path3:(2,3,4)", 3, 21, 4),
#     ("Path4:(1,2,5)", 2, 15, 5),
#     ("Path5:(2,1,5)", 1, 18, 5),
#     ("Path6:(1,1,5)", 1, 18, 5),
#     ("Path7:(1,3,6)", 3, 21, 6),
#     ("Path8:(2,1,6)", 1, 15, 6),
#     ("Path9:(1,3,7)", 3, 21, 7),
#     ("Path10:(2,1,7)", 1, 15, 7),
#     ("Path11:(2,1,8)", 1, 13, 8),
#     ("Path12:(1,3,8)", 3, 15, 8),
#     ("Path13:(2,2,8)", 2, 13, 8)
# ]