using Distances
using JuMP
using SCIP
using COPT
using Cbc

function cal_distance(C,O)
    D = pairwise(Euclidean(), C, O)
    return D
end

function allnodes()
    xy = [
            0	0	0	0	2	5	5	5	5	5	5	5	8	9	10	11	12	12	15	15	15	15	15	15	15	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	25	25	25	25	25	25	25	25	25	28	28	28	28	28	28	28	28	32	32	33	33	33	33	34	34	34	34	34	34	34	35	35	38	38	38	38	40	41	41	41	41	41	48	48	48	51	51	56	57	57	57	61	61	63	64	71	71	71	71	71	74	74	74	74	74	74	74	74	77	78	78	78	78	79	79	79	80	80	80	81	84	84	84	84	84	84	107
            13	26	27	39	0	13	19	25	31	37	43	8	0	10	10	10	10	5	13	19	25	31	37	43	8	11	13	15	17	19	21	23	25	27	29	31	33	35	37	39	41	42	44	45	11	15	22	23	24	26	28	29	9	16	20	28	30	34	40	43	47	26	31	15	26	29	31	15	26	29	31	38	41	5	17	31	16	20	30	34	22	23	32	34	35	36	22	27	6	45	47	25	12	25	44	45	47	6	22	11	13	16	45	47	12	16	20	24	29	35	39	6	21	10	32	35	39	10	33	37	10	41	5	17	20	24	29	34	38	6	27
        ]
    return xy
end

function data(XY)
    NODE = 1:size(XY,2)
    n = length(NODE) - 1
    D = cal_distance(XY, XY)
    return NODE, n, D
end

function setvariables(model, NODE)
    x = @variable(model, [NODE,NODE], base_name = "x", Bin)
    u = @variable(model, [NODE], base_name = "u", lower_bound = 0)
    return x,u
end

function milpmodel(model, NODE, n, D, x, u)
    @objective(model, Min, sum(x[i,j]*D[i,j] for i in NODE for j in NODE if i != j))
    @constraint(model, [i in NODE], sum(x[i,j] for j in NODE if i != j) == 1)
    @constraint(model, [i in NODE], sum(x[j,i] for j in NODE if i != j) == 1)
    @constraint(model, [i in NODE, j in NODE; i != j && j > 1], u[j] - u[i] >= 1 - (n+1)*(1 - x[i,j]))
    @constraint(model, [i in NODE; i > 1], u[i] <= n)
    return model
end

function initial_solution(model, NODE, n, x)
    cons1 = @constraint(model, [i in NODE; i >1], x[i-1,i] == 1, base_name = "cons1")
    cons2 = @constraint(model, x[n+1,1] == 1, base_name = "cons2")
    return cons1, cons2
end

function ini_solution(NODE, n)
    a = zeros(length(NODE), length(NODE))
    for i in NODE
        if i > 1
            a[i-1,i] = 1
        end
    end
    a[n+1,1] = 1
    return a
end

function TSP_OPT_algo()
    # optimizer = SCIP.Optimizer(display_verblevel=4, limits_gap=0.00, lp_threads=5)
    # model = JuMP.direct_model(optimizer)
    # set_time_limit_sec(model, 3600)
    # model = Model(SCIP.Optimizer)
    model = Model(COPT.Optimizer)

    XY = allnodes()
    NODE, n, D = data(XY)
    x, u = setvariables(model, NODE)
    model = milpmodel(model, NODE, n, D, x, u)
    #构造初始解
    initial_path = ini_solution(NODE, n)
    fix.(x, initial_path; force = true)
    optimize!(model)
    #偏优化
    last_obj = objective_value(model)                                   #上次解
    x_value = collect(value.(x))[NODE, NODE]
    wd = 30                                                             #偏优化范围
    icount = 0                                                          #计数器
    while icount <= 50
        fix.(x, x_value; force = true)
        #选择分布图的wd*wd范围区域内进行偏优化
        cur_i = rand(1:n+1)     
        for i in NODE
            for j in NODE
                if i != j && D[i,cur_i] <= wd && D[j,cur_i] <= wd
                    unfix(x[i,j])
                end
            end
        end
        for i in NODE
            for j in NODE
                if i != j && D[i,cur_i] <= wd && D[j,cur_i] > wd && (x_value[i,j] == 1 || x_value[j,i] == 1)
                    unfix(x[i,j])
                end
            end
        end
        optimize!(model)
        if solve_time(model) < 0.2
            wd += 2
        end
        if solve_time(model) > 2 && wd  >= 5
            wd -= 2
        end
        if objective_value(model) < last_obj
            last_obj = objective_value(model)
            icount = 0
            x_value = collect(value.(x))[NODE, NODE]
        else
            icount += 1
        end
    end

end

@time TSP_OPT_algo()

# COPT
# objective = 673.45
# time = 176.528221 seconds (53.11 M allocations: 2.444 GiB, 0.32% gc time, 0.42% compilation time)
# objective = 709.83
# time =  542.313119 seconds (103.99 M allocations: 4.125 GiB, 0.54% gc time, 0.34% compilation time)
# objective = 752.02
# time = 361.945157 seconds (65.26 M allocations: 2.856 GiB, 0.66% gc time, 0.56% compilation time)
# SCIP
# objective = +7.17422163853608e+02
# time = 491.309024 seconds (82.27 M allocations: 2.734 GiB, 0.11% gc time, 0.17% compilation time)