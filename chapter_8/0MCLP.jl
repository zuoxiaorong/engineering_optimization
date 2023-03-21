using JuMP
using SCIP

function xy_d(N)
    params = [   
        1 (63.7, 17.3) (59.2, 47.1) (57.7, 49.9) (51.2, 62.0) (24.8, 89.4) 
        2 (58.6, 10.0) (52.6, 40.7) (96.7, 6.5) (9.7, 39.8) (42.6, 97.6) 
        3 (34.2, 18.6) (37.8, 41.8) (16.5, 48.6) (55.7, 69.9) (96.8, 86.0) 
        4 (54.2, 10.3) (77.1, 39.9) (74.2, 57.3) (13.0, 75.0) (68.0, 85.8) 
        5 (30.7, 11.0) (15.4, 25.6) (1.3, 19.4) (46.9, 60.5) (25.4, 80.7) 
        6 (57.5, 20.3) (73.1, 38.2) (75.3, 46.0) (55.4, 84.3) (15.2, 78.1) 
        7 (42.9, 16.0) (48.1, 48.2) (63.0, 41.2) (20.0, 79.7) (77.0, 71.2) 
        8 (48.4, 20.5) (83.6, 30.5) (37.2, 42.2) (91.9, 64.4) (67.9, 96.3) 
        9 (60.5, 0.4) (14.0, 12.5) (66.3, 43.8) (11.0, 63.9) (81.2, 93.2) 
        10 (53.9, 19.7) (34.7, 30.3) (10.9, 12.9) (22.0, 61.8) (81.4, 87.5) 
        11 (39.3, 16.3) (20.0, 33.0) (15.9, 55.5) (93.1, 52.9) (18.1, 68.3) 
        12 (60.1, 17.7) (11.9, 30.1) (89.9, 46.1) (24.8, 50.5) (18.9, 91.9) 
        13 (27.5, 16.9) (69.3, 24.9) (24.8, 56.1) (47.9, 80.9) (0.9, 94.2) 
        14 (24.0, 5.1) (23.7, 24.5) (93.9, 42.8) (33.7, 55.5) (2.9, 87.9) 
        15 (50.5, 10.2) (92.2, 6.9) (7.3, 27.0) (15.7, 61.5) (31.0, 74.2) 
        16 (32.1, 21.1) (25.1, 74.1) (32.9, 92.8) (40.6, 45.7) (7.0, 51.9)
        17 (38.1, 0.5) (13.3, 32.8) (29.2, 46.2) (83.8, 53.9) (23.3, 67.8)
        18 (75.3, 12.2) (44.4, 37.4) (20.9, 59.2) (76.7, 70.3) (44.1, 72.7)
        19 (62.4, 5.7) (39.0, 31.8) (89.1, 21.4) (70.1, 81.8) (6.4, 99.4)
        20 (43.0, 5.3) (21.1, 27.0) (39.1, 49.3) (74.6, 63.9) (49.4, 85.0)
        21 missing (68.1, 16.5) (54.8, 59.9) (17.1, 47.7) (15.3, 88.4)
        22 missing (24.0, 32.5) (90.2, 15.5) (49.1, 76.5) (12.8, 85.0)
        23 missing (74.9, 11.4) (9.3, 47.1) (98.2, 29.4) (23.9, 85.6)
        24 missing (9.5, 0.2) (7.6, 4.8) (20.0, 58.6) (47.1, 95.0)
        25 missing (10.7, 26.4) (23.2, 40.8) (16.2, 78.4) (74.2, 80.0)
        26 missing missing (29.3, 37.8) (90.2, 41.9) (53.5, 99.8)
        27 missing missing (9.0, 18.7) (12.1, 42.7) (69.2, 75.9)
        28 missing missing (14.8, 19.5) (22.5, 66.1) (89.6, 79.1)
        29 missing missing (22.5, 53.7) (89.6, 62.2) (27.5, 98.2)
        30 missing missing (32.2, 51.9) (45.6, 74.9) (23.8, 94.7)
        31 missing missing missing (49.6, 84.6) (59.5, 94.2)
        32 missing missing missing (12.2, 51.3) (49.1, 82.0)
        33 missing missing missing (89.5, 57.2) (80.8, 80.2)
        34 missing missing missing (84.2, 63.8) (72.0, 86.0)
        35 missing missing missing (38.6, 60.3) (47.6, 88.4)
        36 missing missing missing missing (65.3, 79.3)
        37 missing missing missing missing (2.0, 60.9)
        38 missing missing missing missing (91.6, 61.5)
        39 missing missing missing missing (82.1, 68.4)
        40 missing missing missing missing (80.0, 99.0)
    ]
    x = zeros(size(params, 1), length(N))
    y = zeros(size(params, 1), length(N))
    for t in 1:5
        for n in 1:N[t]
            x[n,t] = params[n, 1+t][1]
            y[n,t] = params[n, 1+t][2]
        end
    end
    return x, y
end

function data()
    T = 1:5                                         #周期的下标
    N = [20, 25, 30, 35, 40]                        #周期t内的行动单元的数量, i，k 行动单元的下标 
    V = 1:3                                         #移动通信车的集合, j 通信车的下标
    R = 15                                          #同通信车j的最大覆盖半径
    L = 30                                          #通信车j的最大移动距离(受限于移动速度)
    C = 18                                          #移动通信车的信道容量
    M = 9999
    X0, Y0 = 50, 0                                  #行动单元初始坐标位置
    x, y = xy_d(N)                                  #行动单元在周期t的坐标位置
    q = 18					                        #given parameter for Euclidean distance linearization
    cita = 0.0872				                    #given parameter for Euclidean distance linearization
    return T, N, V, R, L, C, M, X0, Y0, x, y, q, cita
end

function setvariables(model, T, V, N)
    X = @variable(model, [T,V], base_name = "X", lower_bound = 0)				        #coordinate X
    Y = @variable(model, [T,V], base_name = "Y", lower_bound = 0)				        #coordinate Y
    dx = @variable(model, [t in T,1:N[t],V], base_name = "dx", lower_bound = 0)           #distance in x_d
    dy = @variable(model, [t in T,1:N[t],V], base_name = "dy", lower_bound = 0)           #distance in y_d
    d = @variable(model, [t in T,1:N[t],V], base_name = "d", lower_bound = 0)             #distance between CV and RT
    DX = @variable(model, [T,V], base_name = "DX", lower_bound = 0)                     #distance in x_d
    DY = @variable(model, [T,V], base_name = "DY", lower_bound = 0)                     #distance in y_d
    D = @variable(model, [T,V], base_name = "D", lower_bound = 0)			            #distance between t-1 and t
    E = @variable(model, [t in T,1:N[t],V], base_name = "E", Bin)
    e = @variable(model, [t in T,1:N[t]], base_name = "e", Bin)
    return X, Y, dx, dy, d, DX, DY, D, E, e
end

function obj_maximize_Signal_Coverage(model, T, N, e)
    @objective(model, Max, sum(sum(e[t,i] for i in 1:N[t])/N[t] for t in T))
    return model
end

function obj_minimize_Total_Distance(model, T, V, D)
    @objective(model, Min, sum(D[t,j] for t in T for j in V)) 
    return model
end

function milpmodel(model, T, N, V, R, L, C, M, X0, Y0, x, y, q, cita, X, Y, dx, dy, d, DX, DY, D, E, e)
    @constraint(model, [t in T, i in 1:N[t], j in V], dx[t,i,j] >= X[t,j] - x[i,t])
    @constraint(model, [t in T, i in 1:N[t], j in V], dx[t,i,j] >= x[i,t] - X[t,j])
    @constraint(model, [t in T, i in 1:N[t], j in V], dy[t,i,j] >= Y[t,j] - y[i,t])
    @constraint(model, [t in T, i in 1:N[t], j in V], dy[t,i,j] >= y[i,t] - Y[t,j])
    @constraint(model, [t in T, i in 1:N[t], j in V, p in 0:q-1], d[t,i,j] >= dx[t,i,j]*cos(p*cita) + dy[t,i,j]*sin(p*cita))

    @constraint(model, [j in V], DX[1,j] >= X[1,j] - X0)
    @constraint(model, [j in V], DX[1,j] >= X0 - X[1,j])
    @constraint(model, [t in T, j in V; t>1], DX[t,j] >= X[t-1,j] - X[t,j])
    @constraint(model, [t in T, j in V; t>1], DX[t,j] >= X[t,j] - X[t-1,j])
    @constraint(model, [j in V], DY[1,j] >= Y[1,j] - Y0)
    @constraint(model, [j in V], DY[1,j] >= Y0 - Y[1,j])
    @constraint(model, [t in T, j in V; t>1], DY[t,j] >= Y[t-1,j] - Y[t,j])
    @constraint(model, [t in T, j in V; t>1], DY[t,j] >= Y[t,j] - Y[t-1,j])
    @constraint(model, [t in T, j in V, p in 0:q-1], D[t,j] >= DX[t,j]*cos(p*cita) + DY[t,j]*sin(p*cita))

    @constraint(model, [t in T, j in V], C >= sum(E[t,i,j] for i in 1:N[t]))
    @constraint(model, [t in T, j in V], D[t,j] <= L)
    @constraint(model, [t in T, i in 1:N[t], j in V], M*(E[t,i,j] - 1) <= R - d[t,i,j])
    @constraint(model, [t in T, i in 1:N[t], j in V], e[t,i] <= E[t,i,j])
    return model
end

function  startmodel()
    optimizer = SCIP.Optimizer(display_verblevel=4 , limits_gap=0.00, parallel_maxnthreads=12)
    model = JuMP.direct_model(optimizer)
    set_time_limit_sec(model,1000)

    T, N, V, R, L, C, M, X0, Y0, x, y, q, cita = data()
    X, Y, dx, dy, d, DX, DY, D, E, e = setvariables(model, T, V, N)
    model = milpmodel(model, T, N, V, R, L, C, M, X0, Y0, x, y, q, cita, X, Y, dx, dy, d, DX, DY, D, E, e)
    model = obj_maximize_Signal_Coverage(model, T, N, e)
    optimize!(model)
    cover= value.(e)
    display(cover)
    println(sum(cover))
    for t in 1:5
        a = sum(cover[t,i] for i in 1:N[t])
        println(a/N[t])
    end
    best_obj_SCR = objective_value(model)
    if result_count(model) > 1
        @constraint(model, sum(e[t,i] for t in T for i in 1:N[t]) >= best_obj_SCR)              #把第一个目标函数的求解结果作为约束进行第二轮求解
        model = obj_minimize_Total_Distance(model,  T, V, D)                                    #导入第二个目标函数
        optimize!(model)
        cover= value.(e)
        display(cover)
        println(sum(cover)
        )
    end
end
startmodel()