using JuMP
using SCIP
using Distances

function cal_distance(C,O)
    D = pairwise(Euclidean(), C, O)
    return D
end

function cal_d()
    param_N = [
        1 25 26.4 87.6 
        2 30 23.5 91.9 
        3 45 27.8 93.6 
        4 25 34.6 70.5 
        5 15 40.5 65.2 
        6 45 48.6 19.3 
        7 20 58.5 21.1 
        8 40 57.2 24.5 
        9 10 72.2 75.1 
        10 25 21.7 44.2 
        11 10 24.6 51.3 
        12 30 27.2 39.6
        13 45 22.8 22.5
        14 45 89.5 90.9 
        15 35 95.7 90.6
        16 10 96.8 83.5 
        17 30 78.4 28.9 
        18 35 60.7 80.8 
        19 10 62.7 90.1 
        20 50 45.7 51.4 
    ]
    param_LP = [
        77.5 28.1
        4.50 34.60
        38.5 48.2
        71.2 64.5
        49.4 77.9
    ]
    d = cal_distance(param_N[:,3:4]', param_N[:,3:4]')
    d1 = cal_distance(param_LP', param_N[:,3:4]')
    return param_N[:,2], d, d1
end

function data()
    S = 1:10
    N = 1:20
    P = 1:5
    v = 60
    C = [1.2, 1.3, 0.9, 0.7, 1.8, 1.4, 0.9, 2.2, 3.8, 1.7]
    tao, d, d1 = cal_d() 
    r = [
        1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
        0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
        0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1
        0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 1 
        0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
    ]
    h = 120
    M = 9999
    return S, N, P, v, C, tao, d, d1, r, h, M
end

function lpmodel(model, S, N, P, v, C, tao, d, d1, r, h, M)
    @variable(model, x[N], Bin)                         #设备是否维修
    @variable(model, y[N,P], Bin)                       #设备N是否由P维修
    @variable(model, y1[N,P], Bin)                      #是否为维修路径第一站
    @variable(model, y2[N,P], Bin)                      #是否为维修路径最后一站
    @variable(model, z[N,N], Bin)                       #维修的先后顺序关系
    @variable(model, a[N] >= 0)                         #到达时间
    @variable(model, o[S], Bin)                         #系统是否能恢复
    @variable(model, e[S] >= 0)                         #系统恢复时间
    @variable(model, e1[S] >= 0)                        #系统恢复后的时效能力值

    @objective(model, Max, sum(e1[s] for s in S))
    @constraint(model, [i in N], sum(y[i,p] for p in P) == x[i])                #若维修则分配一个人
    @constraints(model, begin
        #维修第一站
        [i in N, p in P], y1[i,p] <= y[i,p]
        [p in P], sum(y1[i,p] for i in N) <= 1
        #维修最后一站
        [i in N, p in P], y2[i,p] <= y[i,p]
        [p in P], sum(y2[i,p] for i in N) <= 1        
    end)
    @constraints(model, begin
        [i in N], sum(y1[i,p] for p in P) + sum(z[j,i] for j in N if i != j) == x[i]            #最多进站一次
        [i in N], sum(y2[i,p] for p in P) + sum(z[i,j] for j in N if i != j) == x[i]            #最多出站一次
    end)
    @constraint(model, [i in N, j in N, p in P; i != j], y[i,p] - y[j,p] <= 1 - z[i,j])
    @constraint(model, [i in N, j in N; i != j], x[i] + x[j] >= 2*z[i,j])
    @constraint(model, [i in N, p in P], a[i] >= 60*d1[p,i]/v - M*(1 - y1[i,p]))
    @constraint(model, [i in N, j in N; i != j], a[j] >= a[i] + tao[i] + 60*d[i,j]/v - M*(1 - z[i,j]))
    @constraint(model, [i in N], a[i] + tao[i] <= h + M*(1 - x[i]))
    @constraint(model, [i in N], a[i] <= M*x[i])
    @constraint(model, [s in S, i in N; r[s,i] == 1], o[s] <= x[i])
    @constraint(model, [s in S, i in N; r[s,i] == 1], e[s] >= a[i] + tao[i] - M*(1 - x[i]))
    @constraints(model, begin
        [s in S], e1[s] <= C[s]*(h - e[s])
        [s in S], e1[s] <= M*o[s]
    end) 
    return model
end

function  startmodel()
    optimizer = SCIP.Optimizer(display_verblevel=4 , limits_gap=0.00, parallel_maxnthreads=12)
    model = JuMP.direct_model(optimizer)
    set_time_limit_sec(model,1000)

    S, N, P, v, C, tao, d, d1, r, h, M = data()
    model = lpmodel(model, S, N, P, v, C, tao, d, d1, r, h, M)
    optimize!(model)
end
startmodel()