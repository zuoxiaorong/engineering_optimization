using JuMP
using SCIP

function cal_rou(z,x,p)
    rou = (p^x)*((1 - p)^(z-x))*factorial(big(z))/(factorial(big(z-x))*factorial(big(x)))
    return rou
end

function cal_sigma(z,x,p)
    return sum(cal_rou(z,i,p) for i in x:z)
end

function data()
    L = 1:4                             #部署位置集合
    T = 1:3                             #目标集合
    N = 1:4                             #导弹型号集合          
    namda = [                           #导弹平均命中概率
        0.9 0.7 0.8
        0.7 0.7 0.7
        0.8 0.8 0.85
        0.5 0.6 0.65
    ]
    A = [800, 500, 400, 700]            #导弹射程下限
    B = [1200, 800, 600, 900]           #导弹射程上限
    m = [20, 40, 20, 30]                #导弹可用数量
    c = [60, 15, 10, 8]                 #导弹成本单价
    cita = [0.1, 0.3, 0.2, 0.4]         #导弹可靠度
    P = [1, 2, 1.5]                     #目标重要程度（或发生概率）
    K = [3, 2, 2]                       #目标t的打击方案总数
    PI = [                              #打击方案pi：目标t的第k个打击方案采用导弹i，命中n枚完成吕为p（t,k,i,n,p）
        (1, 1, 1),
        (1, 2, 2),
        (1, 3, 4),
        (2, 1, 1),
        (2, 2, 3),
        (3, 1, 3),
        (3, 1, 4),
        (3, 2, 3),
        (3, 2, 2)
    ]
    n = [4, 10, 6, 6, 14, 5, 10, 5, 8]                              #打击方案的最低命中数
    p = [1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 0.4, 0.6, 0.4]               #打击方案的命中概率
    D = [                               #部署位置与目标之间的距离
        870 720 840
        480 580 460
        790 870 460
        620 860 810
    ]
    C = 800
    M = 99999
    return L, T, N, A, B, m, c, P, K, PI, D, C, M, n, p, namda
end

function milpmodel(model, L, T, N, A, B, m, c, P, K, PI, D, C, M, n, p, namda)
    @variable(model, x[L,N,T], Bin)                         #是否将导弹i部署在位置l用于打击目标t
    @variable(model, y[L,N,T] >= 0, Int)                    #表示部署在位置l用于打击目标t的导弹i的数量
    @variable(model, z[N,T] >= 0, Int)                      #表示部署用于打击目标t的导弹i的总数量
    @variable(model, f[t in T, k in 1:K[t]], Bin)           #是否用方案k打击目标t
    @variable(model, o[i in N, T, a in 1:m[i]], Bin)        #表示zit >= a 是否成立
    @variable(model, e[t in T, k in 1:K[t], i in N; (t,k,i) in PI] >= 0)
    @variable(model, E[T] >= 0)                             #打击目标t的完成概率

    @objective(model, Max, sum(P[t]*E[t] for t in T))

    @constraint(model, [j in L, i in N, t in T], y[j,i,t] >= x[j,i,t])
    @constraint(model, [j in L, i in N, t in T], y[j,i,t] <= M*x[j,i,t])
    @constraint(model, [i in N, t in T], sum(y[j,i,t] for j in L) == z[i,t])
    @constraint(model, [i in N], sum(z[i,t] for t in T) <= m[i])
    @constraint(model, sum(y[j,i,t]*c[i] for j in L for i in N for t in T) <= C)
    @constraint(model, [j in L, i in N, t in T], D[j,t] >= A[i] - M*(1 - x[j,i,t]))
    @constraint(model, [j in L, i in N, t in T], D[j,t] <= B[i] + M*(1 - x[j,i,t]))
    @constraint(model, [t in T], sum(f[t,k] for k in 1:K[t]) == 1)
    @constraint(model, [t in T, k in 1:K[t], i in N, g in 1:length(PI); (t,k,i) == PI[g]], z[i,t] >= n[g] - M*(1 - f[t,k]))
    @constraint(model, [i in N, t in T, a in 1:m[i]], o[i,t,a]*M >= a - z[i,t] + 1)
    @constraint(model, [i in N, t in T, a in 1:m[i]], (1 - o[i,t,a])*M >= z[i,t] - a - 1)
    @constraint(model, [t in T, k in 1:K[t], i in N, a in 1:m[i], g in 1:length(PI); (t,k,i) == PI[g] && a >= n[g]], e[t,k,i] <= p[g]*cal_sigma(a,n[g],namda[i,t]) + (2 - f[t,k] - o[i,t,a]))
    @constraint(model, [t in T, k in 1:K[t], i in N, g in 1:length(PI); (t,k,i) == PI[g]], e[t,k,i] <= f[t,k]*p[g])
    @constraint(model, [t in T], E[t] == sum(e[t,k,i] for k in 1:K[t] for i in N if (t,k,i) in PI))
    return model
end

function  startmodel()
    optimizer = SCIP.Optimizer(display_verblevel=4, limits_gap=0.00, parallel_maxnthreads=12)
    model = JuMP.direct_model(optimizer)
    set_time_limit_sec(model,100)

    L, T, N, A, B, m, c, P, K, PI, D, C, M, n, p, namda = data()
    model = milpmodel(model, L, T, N, A, B, m, c, P, K, PI, D, C, M, n, p, namda)
    optimize!(model)
    build_bin = value.(model[:x])
    build_num = value.(model[:y])
    fin_rate = value.(model[:f])
    num = value.(model[:z])
    fin_prob = value.(model[:e])
    total_prob = value.(model[:E])
    dispose_plan = [("L = $(j), DD = $(i), TT = $(t), num = $(build_num[j,i,t])") for j in L for i in N for t in T if build_bin[j,i,t] == 1]
    task_prob = [("t = $(t), k = $(k), i = $(i), n = $(n[g]), z = $(num[i,t]), finishprob = $(fin_prob[t,k,i])") for t in T for k in 1:K[t] for i in N for g in 1:length(PI) if (t,k,i) == PI[g] && fin_rate[t,k] == 1]
    total_prob = [(" t  = $(t), E = $(total_prob[t])") for t in T]
    cost = sum(num[i,t]*c[i] for i in N for t in T)
end
startmodel()

# solution
# obj = +3.18524310521767e+00 (14 solutions)
# dispose_plan
# "L = 1, DD = 4, TT = 1, num = 14.0"
# "L = 1, DD = 4, TT = 3, num = 1.0"
# "L = 3, DD = 3, TT = 3, num = 8.0"
# "L = 4, DD = 1, TT = 2, num = 8.0"
# "L = 4, DD = 4, TT = 3, num = 15.0"
# task_prob
# "t = 1, k = 3, i = 4, n = 6, z = 14.0, finishprob = 0.78802490234375"
# "t = 2, k = 1, i = 1, n = 6, z = 8.0, finishprob = 0.5517738099999998"
# "t = 3, k = 1, i = 3, n = 5, z = 8.0, finishprob = 0.5871885179296874"
# "t = 3, k = 1, i = 4, n = 10, z = 16.0, finishprob = 0.27525853731959593"
# total_prob
# " t  = 1, E = 0.78802490234375"
# " t  = 2, E = 0.5517738099999998"
# " t  = 3, E = 0.8624470552492833"
# cost = 800.0