using JuMP
using SCIP
using Plots

function tangent_linearize_method(l_min, l_max)
    epsilon = 0.06
    mu = (1 + sqrt(epsilon))/(1 - sqrt(epsilon))
    eta = log(mu,l_max/l_min) + 1
    k = [-(mu^(2 - 2p))*(l_min^(-2)) for p in 1:eta]
    b = [2*mu^(1 - p)*(l_min^(-1)) for p in 1:eta]
    return k,b
end

function data()
    DIM = 1:2                                                                           #维度集合,即x和y维度 index by e
    DEP = 1:10                                                                          #部门集合 index by i and j
    DDF = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 10), (4, 6), (4, 10), (5, 10), (6, 7), (6, 9), (7, 9), (8, 10), (9, 10)]   #具有物料流动关系的部门对; #部门之间的物料流量 
    f = [1, 16, 1, 4, 1, 4, 4, 4, 1, 4, 4, 1, 4, 64, 16, 4, 16, 64, 16]
    MET = 1:3                                                                           #距离类型:1-矩形距离;2-欧氏距离;3-天车距离 
    L = [90, 95]                                                                        #车间的长和宽 
    l_min = [                                                                           #部门最短边要求
        20 20
        10 10
        10 10
        10 10
        10 10
        15 15
        10 10
        20 20
        10 10
        20 20
        30 30
    ]   
    l_max = [                                                                           #部门最长边要求 
        90 90
        20 20
        40 40
        50 50
        80 80
        40 40
        50 50
        70 70
        70 70
        100 100
    ]
    a = [1200, 150, 300, 400, 600, 300, 900, 600, 1000, 3000]                           #部门最小面积要求 
    ar = [5,5,5,5,5,5,5,5,5,5]                                                          #部门最大长宽比要求
    #利用切线线性化切分幂函数
    TL_k = [tangent_linearize_method(l_min[i], l_max[i])[1] for i in DEP]               #上述切线的斜率 
    TL_b = [tangent_linearize_method(l_min[i], l_max[i])[2] for i in DEP]               #上述切线的截距 
    TL_n = [length(TL_k[i]) for i in DEP]                                               #部门面积线性化的切线数量,即n 
    TP_cita = 0.0872                                                                    #欧氏距离线性化的0值 
    TP_n = 18                                                                           #欧氏距离线性化的切平面数量 
    w = [                                                                               #部门之间的物料搬运距离类型
        1  0  0
        0  1  0
        0  0  1
        1  0  0
        1  0  0
        0  1  0
        0  1  0
        0  1  0
        1  0  0
        0  0  1
        0  1  0
        1  0  0
        0  1  0
        1  0  0
        0  1  0
        0  0  1
        0  1  0
        1  0  0
        0  1  0
    ] 
    M = 9999                                                                            #一个大数 
    return DIM, DEP, DDF, f, MET, L, l_min, l_max, a, ar, TL_b, TL_k, TL_n, TP_cita, TP_n, w, M
end

function lpmodel(model, DIM, DEP, DDF, f, MET, L, l_min, l_max, a, ar, TL_b, TL_k, TL_n, TP_cita, TP_n, w, M)

    @variable(model, c[DEP,DIM] >= 0)                                                   #部门的中心 
    @variable(model, l[DEP,DIM] >= 0)                                                   #部门的边长 
    @variable(model, s[DEP,DEP,DIM], Bin)                                               #部门之间的相对位置:（1代表i和j在x轴向分开，即i在j的左边; i和j在y轴向分开，即i在j的下边）
    @variable(model, o[DEP,DEP,DIM], Bin)                                               #部门之间在轴向上是否重叠 
    @variable(model, dxy[i in DEP, j in DEP,DIM; (i,j) in DDF] >=0)                     #部门之间的轴向距离 
    @variable(model, dis[i in DEP, j in DEP,MET; (i,j) in DDF] >=0)                     #部门之间的距离 

    @objective(model, Min, sum(f[k]*dis[i,j,d]*w[k,d] for i in DEP for j in DEP for d in MET for k in 1:length(DDF) if (i,j) == DDF[k]))      #minimize objective Cost:


    @constraints(model, begin
        [i in DEP, j in DEP, e in DIM; i != j], c[i,e] + 0.5*l[i,e] <= c[j,e] - 0.5*l[j,e] + M*(1 - s[i,j,e])           #部门面积不能重叠
        [i in DEP, j in DEP; i != j], sum(s[i,j,e] + s[j,i,e] for e in DIM) == 1                                        #i和j在四个方向必须有一个关系
    end)

    @constraints(model, begin
        [i in DEP, j in DEP, e in DIM; (i,j) in DDF], dxy[i,j,e] >= c[i,e] - c[j,e]                                     #部门之间的轴向距离
        [i in DEP, j in DEP, e in DIM; (i,j) in DDF], dxy[i,j,e] >= c[j,e] - c[i,e]                                     #部门之间的轴向距离
    end)

    @constraints(model, begin
        [i in DEP, j in DEP, k in 1:length(DDF); (i,j) == DDF[k] && w[k,1] == 1], dis[i,j,1] >= dxy[i,j,1] + dxy[i,j,2]                      #矩形面积
        [i in DEP, j in DEP, p in 0:TP_n, k in 1:length(DDF); (i,j) == DDF[k] && w[k,2] == 1], dis[i,j,2] >= dxy[i,j,1]*cos(p*TP_cita) + dxy[i,j,2]*sin(p*TP_cita)      #欧氏距离
        [i in DEP, j in DEP, e in DIM, k in 1:length(DDF); (i,j) == DDF[k] && w[k,3] == 1], dis[i,j,3] >= dxy[i,j,e]                         #天车距离    
    end)

    @constraints(model, begin           #采用天车距离的部门对，至少在某一个轴向是重叠的（天车距离）
        [i in DEP, j in DEP, e in DIM, k in 1:length(DDF); (i,j) == DDF[k] && w[k,3] == 1], 0.5*(l[i,e] + l[j,e]) >= dxy[i,j,e] - M*(1 - o[i,j,e])
        [i in DEP, j in DEP, k in 1:length(DDF); (i,j) == DDF[k] && w[k,3] == 1], sum(o[i,j,e] for e in DIM) >= 1
    end)

    @constraints(model,begin           #满足部门面积要求和长宽比要求
        [i in DEP, p in 1:TL_n[i]], l[i,2] >= a[i]*TL_k[i][p]*l[i,1] + a[i]*TL_b[i][p]
        [i in DEP; ar[i] > 0], l[i,1] <= ar[i]*l[i,2]
        [i in DEP; ar[i] > 0], l[i,2] <= ar[i]*l[i,1]        
    end)

    @constraints(model, begin           #部门布局在车间内部
        [i in DEP, e in DIM], c[i,e] + 0.5*l[i,e] <= L[e]
        [i in DEP, e in DIM], c[i,e] - 0.5*l[i,e] >= 0        
    end)

    @constraints(model, begin
        [i in DEP, e in DIM; l_min[i,e] > 0], l[i,e] >= l_min[i,e]
        [i in DEP, e in DIM; l_max[i,e] > 0], l[i,e] <= l_max[i,e]
    end)
    return model
end
function plotresult(c,l)
    rectangle(cx,cy,lx,ly) = Shape(cx .+ [-0.5lx,0.5lx,0.5lx,-0.5lx], cy .+ [-0.5ly,-0.5ly,0.5ly,0.5ly])
    s = plot([rectangle(c[i,1],c[i,2],l[i,1],l[i,2]) for i in 1:10], legend = false, aspect_ratio = 1, alpha = 0.2, grid = 0.2)
    png(s, "out")
end

function startmodel()
    optimizer = SCIP.Optimizer(display_verblevel=4 , limits_gap=0.00, parallel_maxnthreads=12)
    model = JuMP.direct_model(optimizer)
    set_time_limit_sec(model,7200)
    # model = Model(SCIP.Optimizer)
    DIM, DEP, DDF, f, MET, L, l_min, l_max, a, ar, TL_b, TL_k, TL_n, TP_cita, TP_n, w, M = data()
    model = lpmodel(model, DIM, DEP, DDF, f, MET, L, l_min, l_max, a, ar, TL_b, TL_k, TL_n, TP_cita, TP_n, w, M)
    optimize!(model)
    cxy = value.(model[:c])
    lxy = value.(model[:l])
    plotresult(cxy, lxy)    
end

startmodel()
