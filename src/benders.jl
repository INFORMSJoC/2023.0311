#########################################################
################ MASTER PROBLEM
#########################################################

"""
CONSTRUCTOR
Master Problem for Benders decomposition, without any cuts

### Keywords
* `D` - DayBefore
* `time_limit_sec` - time limit for Gurobi optimizer in seconds
* `mip_gap` - optimality tolerance for solver
* `min_workers` - whether to determine min required workers
* `min_cost` - whether to compute min cost
* `no_recourse` - whether to compute min cost without scenario input (straightforward itinerary costs)
### Returns
* JuMP Model
"""
function MasterProblem(
        D::DayBefore,
        time_limit_sec::Int64=MP_TIME_LIMIT,
        mip_gap::Float64=MIP_GAP,
        min_workers::Bool=false,
        min_cost::Bool=false,
        relaxed::Bool=false,
        no_recourse::Bool=false
    )
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB),
    "TimeLimit" => time_limit_sec, "MIPGap" => mip_gap, "OutputFlag" => 0))

    if relaxed
        JuMP.@variable(m, x[D.S] >= 0)
    else
        JuMP.@variable(m, x[D.S], Bin)
    end

    JuMP.@constraints m begin
        [r in D.R], sum(x[i] for i in D.S_r[r]) == 1
    end

    if min_workers
        JuMP.@objective(m, Min, sum(x))
    elseif min_cost
        JuMP.@objective(m, Min, sum(D.l[i] * x[i] for i in D.S))
        JuMP.@constraint(m, sum(x) <= D.W)
    elseif no_recourse
        JuMP.@objective(m, Min, sum(D.l_nr[i] * x[i] for i in D.S))
        JuMP.@constraint(m, sum(x) <= D.W)
    else
        JuMP.@variable(m, theta[D.O])
        JuMP.@constraints m begin
            sum(x) <= D.W
            theta .>= - LARGE_NUMBER    # prevent unboundedness
        end
        JuMP.@objective(m, Min, sum(D.l[i] * x[i] for i in D.S) + sum(D.p[w] * theta[w] for w in D.O))
    end

    return m
end

"""
Solve master problem
### Keywords
* `MP` - JuMP Model
* `max_time` - extended time limit if no solution is found
### Returns
* solve time (seconds)
* objective value
* x values
"""
function solve_MP(
        MP::JuMP.Model,
        max_time::Int64=MP_TIME_LIMIT_MAX
    )
    time = @elapsed JuMP.optimize!(MP)

    # time out -> try one more time
    if (JuMP.termination_status(MP) == MOI.TIME_LIMIT) & (JuMP.primal_status(MP) == MOI.NO_SOLUTION)
        JuMP.set_time_limit_sec(MP, max_time)
        time += @elapsed JuMP.optimize!(MP)
    end

    # cannot proceed without x
    local x_val, obj
    if JuMP.primal_status(MP) == MOI.FEASIBLE_POINT
        x_val = JuMP.value.(MP[:x])
        obj = JuMP.objective_value(MP)
    else
        error("Master problem did not solve correctly. Termination status: " * string(JuMP.termination_status(MP)) * ", primal status: " * string(JuMP.primal_status(MP)))
    end

    return time, obj, x_val
end


################################################
###################### SUBPROBLEMS
################################################

"""
Given an incumbent first-stage solution, retrieve primal second-stage solution for a given scenario

### Keywords
* `w` - Scenario ID
* `D` - DayBefore
* `S_xr` - dictionary of trips to their selected itineraries
* `swapping` - whether to implement swapping recourse formulation
* `time_limit_sec` - time limit in seconds to solve
* `mip_gap` - optimality gap
* `relaxed` - whether to use binary or cts variables
### Returns
* `vals` - tuple of decision variable values
* objective value
"""
function SecondStageScenario(
        w::String,
        D::DayBefore,
        S_xr::Dict{String,String},
        swapping::Bool,
        time_limit::Int64=SP_TIME_LIMIT,
        mip_gap::Float64=MIP_GAP,
        relaxed::Bool=false
    )
    # gather all itineraries
    S_x = intersect(unique([S for S in values(S_xr)]), D.S_w[w])

    # build model
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB), "TimeLimit" => time_limit, "MIPGap" => mip_gap, "OutputFlag" => 0))

    if relaxed
        JuMP.@variables m begin
            y[D.R_w[w]] >= 0
            z[D.A_w[w]] >= 0
        end
    else
        JuMP.@variables m begin
            y[D.R_w[w]], Bin
            z[D.A_w[w]], Bin
        end
    end

    # outsourcing
    JuMP.@constraints m begin
        adhoc_out[r in D.R_w[w]], sum(z[(r, s)] for (q, s) in D.A_w[w] if q == r) <= y[r]
        adhoc_in[r in D.R_w[w]], sum(z[(q, r)] for (q, s) in D.A_w[w] if s == r) <= y[r]
    end

    if swapping
        # swappable requests depend on activated itin. DNS aren't swapped: can be insourced, though
        R_swap = setdiff(unique([r for ((v, i), R_a) in D.R_accept for r in R_a if (v == w) & (i in S_x)]), D.R_w_DNS[w])

        # which itineraries can accept requests (DNS or not)a?
        inswappable = intersect(D.inswappable[w], S_x)

        # scenario amenable to insourcing?
        swappable = length(inswappable) > 0

        # which itineraries have requests that are possible to insource?
        outswappable = [i for i in S_x if length(intersect(R_swap, D.S_list[i])) > 0]

        # itineraries irrelevant to insourcing
        notswappable = setdiff(S_x, union(inswappable, outswappable))

        if swappable
            if relaxed
                JuMP.@variables m begin
                    v[r in R_swap] >= 0
                    a[S in inswappable, r in D.R_accept[(w, S)]] >= 0
                end
            else
                JuMP.@variables m begin
                    v[r in R_swap], Bin
                    a[S in inswappable, r in D.R_accept[(w, S)]], Bin
                end
            end
        end

        # single modifications
        for S in setdiff(inswappable, outswappable)
            JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[S], D.R_w[w]), D.R_w_DNS[w])) +
            sum(a[S, r] for r in D.R_accept[(w, S)]) <= 1)
        end
        for S in setdiff(outswappable, inswappable)
            JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[S], D.R_w[w]), D.R_w_DNS[w])) +
            sum(v[r] for r in intersect(D.S_list[S], R_swap)) <= 1)
        end
        for S in intersect(inswappable, outswappable)
            JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[S], D.R_w[w]), D.R_w_DNS[w])) +
            sum(v[r] for r in intersect(D.S_list[S], R_swap)) +
            sum(a[S, r] for r in D.R_accept[(w, S)]) <= 1)
        end
        for S in notswappable
            JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[S], D.R_w[w]), D.R_w_DNS[w])) <= 1)
        end

        if swappable
            # flow balance
            JuMP.@constraints m begin
                swap_flowbal[r in R_swap], v[r] == sum(a[S, r] for S in inswappable if r in D.R_accept[(w, S)])
            end

            # are there DNS requests?
            if isempty(D.R_w_DNS[w])
                JuMP.@objective(m, Min,
                - sum(
                    # remove insourced requests
                    D.remove[D.S_effective[w, S_xr[r]], r] * v[r] for r in R_swap
                ) + sum(
                    # accept insourced requests
                    D.accept[D.S_effective[w, S], r] * a[S, r] for S in inswappable for r in D.R_accept[(w, S)]
                ) +
                    # remove outsourced requests
                    sum((D.f + D.v[r] - D.remove[D.S_effective[w, S_xr[r]], r]) * y[r] for r in D.R_w[w]) +
                    # construct ad hoc itineraries
                    sum((D.c[e] - D.f) * z[e] for e in D.A_w[w])
                )
            else
                # ensure DNS requests aren't cancelled
                R_DNS_accept = intersect([r for S in inswappable for r in D.R_accept[w, S]], D.R_w_DNS[w])
                JuMP.@constraints m begin
                    [r in setdiff(D.R_w_DNS[w], R_DNS_accept)], y[r] == 1
                    [r in R_DNS_accept], y[r] + sum(a[S, r] for S in inswappable if r in D.R_accept[w, S]) == 1
                end
                JuMP.@objective(m, Min,
                    # same as before
                    - sum(D.remove[D.S_effective[w, S_xr[r]], r] * v[r] for r in R_swap) +
                    sum(D.accept[D.S_effective[w, S], r] * a[S, r] for S in inswappable for r in D.R_accept[w, S]) +
                    sum((D.f + D.v[r]) * y[r] for r in D.R_w[w]) -
                    # only decide to remove non-DNS requests
                    sum(D.remove[D.S_effective[w, S_xr[r]], r] * y[r] for r in setdiff(D.R_w[w], D.R_w_DNS[w])) +
                    sum((D.c[e] - D.f) * z[e] for e in D.A_w[w])
                )
            end
        else
            # it's just outsourcing even though we wanted to swap
            if isempty(D.R_w_DNS[w])
                JuMP.@objective(m, Min,
                    sum((D.f + D.v[r] - D.remove[D.S_effective[w, S_xr[r]], r]) * y[r] for r in D.R_w[w]) +
                    sum((D.c[e] - D.f) * z[e] for e in D.A_w[w])
                )
            else
                # ensure DNS requests aren't cancelled
                JuMP.@constraints m begin
                    [r in D.R_w_DNS[w]], y[r] == 1
                end
               JuMP.@objective(m, Min,
                    sum((D.f + D.v[r]) * y[r] for r in D.R_w[w]) -
                    sum(D.remove[D.S_effective[w, S_xr[r]], r] * y[r] for r in setdiff(D.R_w[w], D.R_w_DNS[w])) +
                    sum((D.c[e] - D.f) * z[e] for e in D.A_w[w])
                )
            end
        end
    else
        # one modification
        JuMP.@constraint(m, [S in S_x], sum(y[r] for r in setdiff(intersect(D.S_list[S], D.R_w[w]), D.R_w_DNS[w])) <= 1)

        if isempty(D.R_w_DNS[w])
            JuMP.@objective(m, Min,
                sum((D.f + D.v[r] - D.remove[D.S_effective[w, S_xr[r]], r]) * y[r] for r in D.R_w[w]) +
                sum((D.c[e] - D.f) * z[e] for e in D.A_w[w])
            )
        else
            # ensure DNS requests aren't cancelled
            JuMP.@constraints m begin
                [r in D.R_w_DNS[w]], y[r] == 1
            end
            JuMP.@objective(m, Min,
                sum((D.f + D.v[r]) * y[r] for r in D.R_w[w]) -
                sum(D.remove[D.S_effective[w, S_xr[r]], r] * y[r] for r in setdiff(D.R_w[w], D.R_w_DNS[w])) +
                sum((D.c[e] - D.f) * z[e] for e in D.A_w[w])
            )
        end
    end

    # obtain optimal decisions
    JuMP.optimize!(m)
    is_opt = JuMP.termination_status(m) == MOI.OPTIMAL
    y_val = is_opt ? JuMP.value.(y) : nothing
    z_val = is_opt ? JuMP.value.(z) : nothing
    obj = is_opt ? JuMP.objective_value(m) : nothing
    v_val = is_opt & (:v in keys(m.obj_dict)) ? JuMP.value.(v) : nothing
    a_val = is_opt & (:a in keys(m.obj_dict)) ? JuMP.value.(a) : nothing
    vals = (y_val, z_val, v_val, a_val)

    return vals, obj
end

"""
Given an incumbent first-stage solution, retrieve primal second-stage solutions

### Keywords
* `D` - DayBefore
* `S_xr` - dictionary of trips to their selected itineraries
* `swapping` - whether to implement swapping recourse formulation
* `time_limit_sec` - time limit in seconds to solve
* `mip_gap` - optimality gap
* `relaxed` - whether to use binary or cts variables
### Returns
* `vals` - dictionary (scenario_id => tuple of decision variable values)
* objective value
"""
function SecondStageSolutions(
        D::DayBefore,
        S_xr::Dict{String,String},
        swapping::Bool=false,
        time_limit::Int64=SP_TIME_LIMIT * 10,
        mip_gap::Float64=MIP_GAP,
        relaxed::Bool=false
    )
    obj_scen = Dict()
    vals = Dict()
    for w in D.O
        vals[w], obj_scen[w] = SecondStageScenario(w, D, S_xr, swapping, time_limit, mip_gap, relaxed)
    end
    obj = sum(D.p[w] * obj_scen[w] for w in D.O)

    return vals, obj
end

"""
CONSTRUCTOR
Full Benders subproblem

### Keywords
* `x_val` - first-stage decision
* `w` - scenario
* `D` - DayBefore
* `time_limit` - in seconds
* `swapping` - whether to implement swapping formulation
### Returns
* Dual model
"""
function FullSubproblemDual(
        x_val,
        w::String,
        D::DayBefore,
        time_limit::Int64=SP_TIME_LIMIT,
        swapping::Bool=false
    )
    # build
    SP_dual = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB), "OutputFlag" => 0, "DualReductions" => 0, "TimeLimit" => time_limit))

    JuMP.@variables SP_dual begin
        pi[D.R_w[w]] >= 0
        rho[D.R_w[w]] >= 0
        sigma[D.S_w[w]] >= 0
    end

    JuMP.@constraints SP_dual begin
        [(r, s) in D.A_w[w]], pi[r] + rho[s] >= D.f - D.c[(r, s)]
        [i in D.S_w[w], r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])], pi[r] + rho[r] - sigma[i] <= D.f + D.v[r] - D.remove[D.S_effective[w, i], r]
    end

    # DNS outsourcing
    if !isempty(D.R_w_DNS[w])
        JuMP.@variable(SP_dual, gamma[D.R_w_DNS[w]])
        JuMP.@constraint(SP_dual, [r in D.R_w_DNS[w]], pi[r] + rho[r] + gamma[r] <= D.f + D.v[r])
        JuMP.@objective(SP_dual, Max, sum(-x_val[i] * sigma[i] for i in D.S_w[w]) + sum(gamma))
    else
        JuMP.@objective(SP_dual, Max, sum(-x_val[i] * sigma[i] for i in D.S_w[w]))
    end

    if swapping & D.insourcing[w]
        JuMP.@variable(SP_dual, tau[D.R_swap[w]])
        JuMP.@constraints SP_dual begin
            # v variables
            [i in D.outswappable[w], r in intersect(D.S_list[i], D.R_swap[w])], -sigma[i] + tau[r] <= -D.remove[D.S_effective[w, i], r]
            # a variables
            [i in D.inswappable[w], r in setdiff(D.R_accept[w, i], D.R_w_DNS[w])], - sigma[i] - tau[r] <= D.accept[D.S_effective[w, i], r]
        end
        if !isempty(D.R_w_DNS[w])
            JuMP.@constraint(SP_dual, [i in D.inswappable[w], r in intersect(D.R_accept[w, i], D.R_w_DNS[w])], - sigma[i] + gamma[r] <= D.accept[D.S_effective[w, i], r])
        end
    end

    return SP_dual
end

"""
CONSTRUCTOR
Activated Benders subproblem

### Keywords
* `S_x` - Vector - itineraries chosen
* `S_xr` - Dict - request => itinerary assignment
* `w` - Scenario ID
* `D` - DayBefore
* `time_limit` - upper bound on solve time in seconds
* `swapping` - whether to build swapping version
### Returns
* activated subproblem dual
"""
function ActivatedSubproblemDual(
        S_x::Vector{String},
        S_xr::Dict,
        w::String,
        D::DayBefore,
        time_limit::Int64=SP_TIME_LIMIT,
        swapping::Bool=false
    )
    # swapping index HELL
    if swapping
        # swappable requests depend on activated itin. DNS aren't swapped: can be insourced, though
        R_swap = setdiff(unique([r for ((v, i), R_a) in D.R_accept for r in R_a if (v == w) & (i in S_x)]), D.R_w_DNS[w])

        # which itineraries can accept requests (DNS or not)a?
        inswappable = intersect(D.inswappable[w], S_x)

        # scenario amenable to insourcing?
        insourcing = length(inswappable) > 0

        # which itineraries have requests that are possible to insource?
        outswappable = [i for i in S_x if length(intersect(R_swap, D.S_list[i])) > 0]

        # itineraries irrelevant to insourcing
        notswappable = setdiff(S_x, union(inswappable, outswappable))
    end

    # build subproblem
    SP_dual = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB), "OutputFlag" => 0,
     "DualReductions" => 0, "TimeLimit" => time_limit))

    JuMP.@variables SP_dual begin
        pi[D.R_w[w]] >= 0
        rho[D.R_w[w]] >= 0
        sigma[intersect(S_x, D.S_w[w])] >= 0
    end

    # constraints corresponding to primal variables y and z (outsourcing)
    JuMP.@constraints SP_dual begin
        [(r, s) in D.A_w[w]], pi[r] + rho[s] >= D.f - D.c[(r, s)]
        [r in setdiff(D.R_w[w], D.R_w_DNS[w])], pi[r] + rho[r] - sigma[S_xr[r]] <= D.f + D.v[r] - D.remove[D.S_effective[w, S_xr[r]], r]
    end

    # DNS outsourcing
    if !isempty(D.R_w_DNS[w])
        JuMP.@variable(SP_dual, gamma[D.R_w_DNS[w]])
        JuMP.@constraint(SP_dual, [r in D.R_w_DNS[w]], pi[r] + rho[r] + gamma[r] <= D.f + D.v[r])
        JuMP.@objective(SP_dual, Max, - sum(sigma) + sum(gamma))
    else
        JuMP.@objective(SP_dual, Max, - sum(sigma))
    end

    if swapping
        if insourcing
            # swapping flow balance
            JuMP.@variable(SP_dual, tau[R_swap])
            JuMP.@constraints SP_dual begin
                [r in R_swap], -sigma[S_xr[r]] + tau[r] <= -D.remove[D.S_effective[w, S_xr[r]], r]
                [i in inswappable, r in setdiff(D.R_accept[w, i], D.R_w_DNS[w])], -sigma[i] - tau[r] <= D.accept[D.S_effective[w, i], r]
            end
            if !isempty(D.R_w_DNS[w])
                JuMP.@constraint(SP_dual, [i in inswappable, r in intersect(D.R_accept[w, i], D.R_w_DNS[w])], - sigma[i] + gamma[r] <= D.accept[D.S_effective[w, i], r])
            end
        end
    end

    return SP_dual
end

"""
Solve either the full or activated dual subproblem.

### Keywords
* `SP_dual` - subproblem dual JuMP Model
### Returns
* solve time in seconds
* objective value
* whether termination status is optimal
* (pi, rho, sigma, tau, gamma) values, or (pi, rho, sigma, tau) if swapping + insourcing was possible
"""
function solve_SP(
        SP_dual::JuMP.Model
    )
    time = @elapsed JuMP.optimize!(SP_dual)
    sp_optimal = JuMP.termination_status(SP_dual) == MOI.OPTIMAL
    if sp_optimal
        # get values
        obj = JuMP.objective_value(SP_dual)
        pi_val = JuMP.value.(SP_dual[:pi])
        rho_val = JuMP.value.(SP_dual[:rho])
        sigma_val = JuMP.value.(SP_dual[:sigma])
        tau_val = :tau in keys(SP_dual.obj_dict) ? JuMP.value.(SP_dual[:tau]) : nothing
        gamma_val = :gamma in keys(SP_dual.obj_dict) ? JuMP.value.(SP_dual[:gamma]) : nothing
        vals = (pi_val, rho_val, sigma_val, tau_val, gamma_val)
    else
        @info "Subproblem did not solve correctly. Termination status: " * string(JuMP.termination_status(SP_dual)) * ", primal status: " * string(JuMP.primal_status(SP_dual)) * ", solve time:" * string(time)
        obj, vals = nothing, nothing
    end
    return time, obj, sp_optimal, vals
end

################################################
###################### ACCELERATION
################################################

struct ParetoArguments
    pareto::Bool
    localpareto::Bool
    activated::Bool
    n::Int64
    dynamic_n::Int64
end


"""
Compute a core point for a master problem corresponding to a case study.
Return -1 as solve time if there is no core point.

### Keywords
* `D` - DayBefore
* `S_all` - set of itineraries over which to retrieve core point
* `time_limit` - solve time limit
* `mip_gap` - integrality gap
### Returns
* objective value if model solves, else nothing
* core point if one exists, else nothing
* solve time in seconds
"""
function core_point(
        D::DayBefore,
        S_all::Vector{String},
        time_limit::Int64=MP_TIME_LIMIT,
        mip_gap::Float64=MIP_GAP
    )
    @assert issubset(S_all, D.S)
    @assert all([length(intersect(D.S_r[r], S_all)) > 0 for r in D.R])

    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB),
        "TimeLimit" => time_limit,
        "MIPGap" => mip_gap,
        "OutputFlag" => 0))

    JuMP.@variables m begin
        x[S_all] >= 0
        zR[D.R]
        zW
        z
    end

    JuMP.@constraints m begin
        [r in D.R], sum(x[i] for i in intersect(D.S_r[r], S_all)) + zR[r] <= 1
        sum(x) + zW <= D.W
        z .<= x
        z .<= zR
        z <= zW
    end

    JuMP.@objective(m, Max, z)
    st = @elapsed JuMP.optimize!(m)
    if JuMP.termination_status(m) != MOI.OPTIMAL
        return nothing, nothing, -1
    elseif JuMP.objective_value(m) > 0
        return JuMP.objective_value(m), JuMP.value.(x), st
    end
    return 0, nothing, -1
end

"""
Compute an activated point.

If solve time returns as -1, there is no core point.

## Keywords
* `D` - DayBefore
* `n` - Number of feasible solutions to average
* `time_limit` - Solve time limit (seconds)
* `mip_gap` - Epsilon tolerance for integrality gap
### Returns
* optimal solution to core point model
* list of itineraries used in solution
* solve time in seconds
"""
function act_core_point(
        D::DayBefore,
        n::Int64,
        time_limit::Int64=MP_TIME_LIMIT,
        mip_gap::Float64=MIP_GAP
    )
    S_all = String[]

    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB),
        "TimeLimit" => time_limit,
        "MIPGap" => mip_gap,
        "OutputFlag" => 0))

    JuMP.@variable(m, x[D.S], Bin)

    JuMP.@constraints m begin
        [r in D.R], sum(x[i] for i in D.S_r[r]) == 1
        sum(x) <= D.W
    end

    #--- stage 1: obtain itinerary set
    for j in 1:n
        JuMP.optimize!(m)
        optimal = JuMP.termination_status(m) == MOI.OPTIMAL
        if !optimal & (j == 1)
            @error "No feasible MP solution..." n j
            return nothing, nothing, -1
        elseif !optimal & (j == 2)
            @info "No core point. Switching to regular ABD." n j
            return nothing, nothing, -1
        elseif !optimal
            @info "Terminating local branching prematurely." n j
            break
        end
        x_val = JuMP.value.(x)

        Ix = [i for i in D.S if x_val[i] > 0.5]
        S_all = vcat(S_all, Ix)

        # rule out previous optimal solution, obtain second optimal solution
        JuMP.@constraint(m, [i in Ix], x[i] == 0)
    end

    #--- stage 2: compute core point corresponding to the set
    obj, x0, st = core_point(D, S_all)
    if obj < 1e-8
        @info "No core point given itinerary set:" S_all
        st = -1
    end

    return x0, S_all, st
end

"""
Dynamically compute an activated core point.

## Keywords
* `D` - DayBefore
* `S_act` - list of itinerary lists used in previous `n` first-stage solutions
* `n` - Rolling window length
### Returns
* optimal obj value of core point model
* optimal solution to core point model
* solve time in seconds
* list of itineraries used in solution
"""
function act_core_point_dynamic(
        D::DayBefore,
        S_act::Vector{Vector{String}},
        n::Int64
    )
    S_cur = unique([S for S_list in S_act[(end - min(n, length(S_act)) + 1):end] for S in S_list])
    obj, x0, cp_st = core_point(D, S_cur)
    return obj, x0, cp_st, S_cur
end

"""
Build Pareto dual. First solve regular dual subproblem.
Then build Pareto dual subproblem using core point and optimal value.

### Keywords
* `x0` - Core point
* `SP_dual` - Regular subproblem dual
* `z_star` - optimal objective value of regular dual
* `S_core` - Set of itineraries represented in core point
* `x_val` - First-stage incumbent decision (JuMP container)
* `w` - Scenario ID
* `D` - DayBefore
* `time_limit` - Computation time limit in seconds
* `swapping` - Whether to implement swapping formulation
### Returns
* JuMP Model - Full Pareto dual
"""
function FullSubproblemParetoDual(
        x0,
        SP_dual,
        z_star::Float64,
        S_core::Vector{String},
        x_val,
        w::String,
        D::DayBefore
    )
    DNS = !isempty(D.R_w_DNS[w])

    # add Pareto changes
    sigma = SP_dual[:sigma]
    if DNS
        gamma = SP_dual[:gamma]
        JuMP.@constraint(SP_dual, pareto, - sum(x_val[i] * sigma[i] for i in D.S_w[w]) + sum(gamma) == z_star)
    else
        JuMP.@constraint(SP_dual, pareto, - sum(x_val[i] * sigma[i] for i in D.S_w[w]) == z_star)
    end

    if !isempty(intersect(D.S_w[w], S_core))
        if DNS
            JuMP.@objective(SP_dual, Max, - sum(x0[i] * sigma[i] for i in intersect(D.S_w[w], S_core)) + sum(gamma))
        else
            JuMP.@objective(SP_dual, Max, - sum(x0[i] * sigma[i] for i in intersect(D.S_w[w], S_core)))
        end
    end
    return SP_dual
end

"""
Build full subproblem Pareto dual when regular full subproblem has not yet been built (for use in ABD)

### Keywords
* `x0` - Core point
* `z_star` - Optimal objective value of regular subproblem dual
* `SP_dual` - Regular subproblem dual
* `z_star` - optimal objective value of regular dual
* `S_core` - Set of itineraries represented in core point
* `x_val` - First-stage incumbent decision (JuMP container)
* `w` - Scenario ID
* `D` - DayBefore
* `time_limit` - Computation time limit in seconds
* `swapping` - Whether to implement swapping formulation
### Returns
* JuMP Model - Full Pareto dual
"""
function FullSubproblemParetoDual(
        x0,
        z_star::Float64,
        S_core::Vector{String},
        x_val,
        w::String,
        D::DayBefore,
        time_limit::Int64=SP_TIME_LIMIT,
        swapping::Bool=false
    )
    SP_dual = FullSubproblemDual(x_val, w, D, time_limit, swapping)
    S_x, S_xr = get_roster(D, x_val)
    return FullSubproblemParetoDual(x0, SP_dual, z_star, S_core, x_val, w, D)
end

"""
Build activated Pareto dual subproblem using (activated) core point x0 and optimal value.

### Keywords
* `x0` - (Activated) core point
* `z_star` - optimal objective value of regular dual
* `S_core` - Set of itineraries represented in (act) core point
* `S_x` - Vector - itineraries chosen
* `S_xr` - Dict - request => itinerary assignment
* `w` - Scenario ID
* `D` - DayBefore
* `time_limit` - Computation time limit in seconds
* `swapping` - Whether to implement swapping formulation
### Returns
* JuMP Model - Activated Pareto dual
"""
function ActivatedSubproblemParetoDual(
        x0,
        z_star::Float64,
        S_core::Vector{String},
        S_x::Vector{String},
        S_xr::Dict,
        w::String,
        D::DayBefore,
        time_limit::Int64=SP_TIME_LIMIT,
        swapping::Bool=false
    )
    #------- Hold core point objects
    S_all = unique(vcat(S_core, S_x))
    S_all_w = intersect(D.S_w[w], S_all)

    #------- Build activated Pareto dual subproblem
    SP_dual = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB),
        "TimeLimit" => time_limit, "OutputFlag" => 0)
        )

    JuMP.@variables SP_dual begin
        pi[D.R_w[w]] >= 0
        rho[D.R_w[w]] >= 0
        sigma[S_all_w] >= 0
    end

    if isempty(D.R_w_DNS[w])
        # no gamma values if no DNS requests
        JuMP.@constraint(SP_dual, pareto, - sum(sigma[i] for i in intersect(S_x, D.S_w[w])) == z_star)
        if !isempty(intersect(S_core, D.S_w[w]))
            JuMP.@objective(SP_dual, Max, sum(-x0[i] * sigma[i] for i in intersect(S_core, D.S_w[w])))
        end
    else
        JuMP.@variable(SP_dual, gamma[D.R_w_DNS[w]])
        JuMP.@constraints SP_dual begin
            pareto, - sum(sigma[i] for i in intersect(S_x, D.S_w[w])) + sum(gamma) == z_star
            [r in D.R_w_DNS[w]], pi[r] + rho[r] + gamma[r] <= D.f + D.v[r]
        end
        if !isempty(intersect(S_core, D.S_w[w]))
            JuMP.@objective(SP_dual, Max, sum(-x0[i] * sigma[i] for i in intersect(S_core, D.S_w[w])) + sum(gamma))
        else
            JuMP.@objective(SP_dual, Max, sum(gamma))
        end
    end

    JuMP.@constraints SP_dual begin
        [(r, s) in D.A_w[w]], pi[r] + rho[s] >= D.f - D.c[r, s]
        [i in S_all_w, r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])], pi[r] + rho[r] - sigma[i] <= D.f + D.v[r] - D.remove[D.S_effective[w, i], r]
    end

    if swapping
        inswappable = intersect(D.inswappable[w], S_all_w)
        R_swap = setdiff(unique([r for ((v, i), R_a) in D.R_accept for r in R_a if (v == w) & (i in S_all_w)]), D.R_w_DNS[w])
        outswappable = [i for i in S_all_w if length(intersect(R_swap, D.S_list[i])) > 0]
        notswappable = setdiff(S_all_w, union(inswappable, outswappable))

        if length(inswappable) > 0
            JuMP.@variable(SP_dual, tau[R_swap])
            JuMP.@constraints SP_dual begin
                [i in outswappable, r in intersect(D.S_list[i], R_swap)], -sigma[i] + tau[r] <= -D.remove[D.S_effective[w, i], r]
                [i in inswappable, r in setdiff(D.R_accept[w, i], D.R_w_DNS[w])], -sigma[i] - tau[r] <= D.accept[D.S_effective[w, i], r]
            end
            if !isempty(D.R_w_DNS[w])
                JuMP.@constraint(SP_dual, [i in inswappable, r in intersect(D.R_accept[w, i], D.R_w_DNS[w])], - sigma[i] + gamma[r] <= D.accept[D.S_effective[w, i], r])
            end
        end
    end
    return SP_dual
end


################################################
###################### BENDERS
################################################


"""
Benders decomposition

### Keywords
* `activated` - whether to implement ABD (else BD)
* `D` - DayBefore
* `swapping` - whether to implement swapping formulation (else non-swapping)
* `pargs` - Pareto args
* `max_iterations` - max Benders iterations before timeout
* `mp_time_limit` - master problem time limit per iteration in seconds
* `sp_time_limit` - subproblem time limit per iteration/scenario in seconds
* `total_time_limit` - total time limit for algorithm in seconds
* `mip_gap` - master problem and subproblem optimality gap tolerance
* `benders_gap` - benders convergence criterion
### Returns
* lower bounds (sequential)
* upper bounds (sequential)
* mip bounds (sequential)
* MP times (sequential)
* max SP times (seqential)
* total SP times (sequential)
* first-stage decisions
* second-stage costs
* second-stage decisions
* whether each iteration was a speedy iteration (valid LB?)
"""
function benders(
        activated::Bool,
        D::DayBefore,
        swapping::Bool=false,
        pargs::ParetoArguments=NULLPARETO,
        speedy_mip_gap::Float64=SPEEDY_MIP_GAP,
        mip_gap::Float64=MIP_GAP,
        benders_gap::Float64=MIP_GAP,
        debug_label::String=DEBUG_LABEL,
        total_time_limit::Int64=TIME_LIMIT,
        max_iterations::Int64=MAX_ITER,
        max_speedy_iterations::Int64=MAX_SPEEDY_ITER,
        mp_time_limit::Int64=MP_TIME_LIMIT,
        sp_time_limit::Int64=SP_TIME_LIMIT
    )
    # first solve is to full optimality
    MP = MasterProblem(D, mp_time_limit, mip_gap)
    JuMP.@objective(MP, Min, sum(D.l[i] * MP[:x][i] for i in D.S))
    local speedy = (speedy_mip_gap > mip_gap)

    # obtain core point
    local x0
    local S_core
    if pargs.pareto
        S_core = D.S
        _, x0, cp_solvetime = core_point(D, S_core)
        if cp_solvetime < 0 # no core point
            pargs = NULLPARETO
            println("No activated core point, switching to regular ABD")
        end
    elseif pargs.localpareto
        x0, S_core, cp_solvetime = act_core_point(D, pargs.n)
        if cp_solvetime < 0 # no core point
            pargs = NULLPARETO
            println("No activated core point, switching to no accel")
        end
    end

    lower_bounds = []; upper_bounds = []; second_stage = []
    MP_times = []; SP_max_times = []; SP_times = []; mip_bounds = []
    is_speedy = []
    if pargs.dynamic_n > 0
        S_act = [S_core]
    end
    local MP_time # keep most recent MP time
    local x_incumbent = nothing
    local S_xr
    local iteration = 0
    local total_time = 0
    while true
        println("Iteration ", iteration += 1, "...")

        #--- main problem
        if (length(MP_times) > 0) & speedy
            # if this isn't the first solve and we are being speedy
            JuMP.set_optimizer_attribute(MP, "MIPGap", speedy_mip_gap)
        end
        MP_time, lb, x_val = solve_MP(MP)
        S_x, S_xr = get_roster(D, x_val)
        push!(is_speedy, speedy)

        if (pargs.pareto | pargs.localpareto) & (length(MP_times) == 0)
            # initial core point solve time
            MP_time += cp_solvetime
        end
        if iteration == 1
            # very low lb
            lb = - LARGE_NUMBER
        elseif iteration == 2
            # two-stage objective
            JuMP.@objective(MP, Min, sum(D.l[i] * MP[:x][i] for i in D.S) + sum(D.p[w] * MP[:theta][w] for w in D.O))
        end
        total_time += MP_time
        new_lb = speedy ? maximum(vcat(lower_bounds, [lb])) : lb
        push!(lower_bounds, new_lb)
        push!(MP_times, MP_time)

        #--- subproblems
        obj_SP = Dict()
        SP_time_all = []
        for w in D.O
            println("Scenario ", w, "...")
            SP_solvetime = 0

            # build and solve regular dual subproblem
            if activated
                SP_dual = ActivatedSubproblemDual(S_x, S_xr, w, D, sp_time_limit, swapping)
            else
                SP_dual = FullSubproblemDual(x_val, w, D, sp_time_limit, swapping)
            end
            SP_time, SP_obj, SP_optimal, SP_vals = solve_SP(SP_dual)
            SP_solvetime += SP_time

            if SP_optimal
                obj_SP[w] = SP_obj
                pi_val, rho_val, sigma_val, tau_val, gamma_val = SP_vals
            else
                bad_x = DataFrame([[i] for i in S_x])
                CSV.write(debug_label * "scenario_" * w * ".csv", bad_x)
                @info("Infeasible dual in scenario " * w * " at iteration " * string(iteration) * "\n")
            end

            # build and solve pareto dual
            if (pargs.localpareto | pargs.pareto) & SP_optimal
                if pargs.localpareto & pargs.activated
                    # A-LPO
                    SP_dual = ActivatedSubproblemParetoDual(x0, SP_obj, S_core, S_x, S_xr, w, D, sp_time_limit, swapping)
                elseif activated
                    # LPO or PO with ABD - need to build full SP dual
                    SP_dual = FullSubproblemParetoDual(x0, SP_obj, S_core, x_val, w, D, sp_time_limit, swapping)
                else
                    # LPO or PO with BD - add constraints to full SP dual
                    SP_dual = FullSubproblemParetoDual(x0, SP_dual, SP_obj, S_core, x_val, w, D)
                end
                SP_time, _, SP0_optimal, SP0_vals = solve_SP(SP_dual)
                SP_solvetime += SP_time

                if SP0_optimal
                    # replace second-stage values
                    pi_val, rho_val, sigma_val, tau_val, gamma_val = SP0_vals
                else
                    # we can still add the non-Pareto-optimized cut
                    bad_x = DataFrame([[i] for i in S_x])
                    bad_x0 = DataFrame([[i] for i in S_core])
                    CSV.write(debug_label * "scenario_" * w * ".csv", bad_x)
                    CSV.write(debug_label * "scenario_" * w * "_x0.csv", bad_x0)
                    @info("Feasible regular dual and infeasible Pareto dual in scenario " * w * " at iteration " * string(iteration) * "\n")
                end
            end
            total_time += SP_solvetime
            push!(SP_time_all, SP_solvetime)

            # add cut
            theta, x = MP[:theta][w], MP[:x]
            temp_deactivate = false
            if isnothing(SP_vals)
                @info "subproblem termination status is not optimal, no cuts added in scenario " * w * " at iteration " * string(iteration) * "\n"
            elseif (activated | pargs.activated)
                #--- construct full subproblem solution from activated solution
                if pargs.pareto | pargs.localpareto
                    if SP0_optimal
                        S_omit = vcat(S_x, S_core)
                    elseif activated
                        S_omit = S_x
                    else
                        temp_deactivate = true
                    end
                else
                    S_omit = S_x
                end
                if !temp_deactivate
                    # values that contribute to sigma_hat no matter what
                    sigma_hat = Dict(i => vcat([0], [pi_val[r] + rho_val[r] - D.f - D.v[r] + D.remove[D.S_effective[w, i], r] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])]) for i in setdiff(D.S_w[w], S_omit))
                    if !isnothing(tau_val)
                        # swaps considered as fcn of 1st stage decision
                        R_swap_x = [r[1] for r in keys(tau_val)]
                        for i in setdiff(D.outswappable[w], S_omit)
                            sigma_hat[i] = vcat(sigma_hat[i],
                                # insourcing removal
                                [tau_val[r] + D.remove[D.S_effective[w, i], r] for r in intersect(D.S_list[i], R_swap_x)],
                                [D.remove[D.S_effective[w, i], r] for r in setdiff(intersect(D.R_swap[w], D.S_list[i]), R_swap_x)]
                            )
                        end
                        for i in setdiff(D.inswappable[w], S_omit)
                            sigma_hat[i] = vcat(sigma_hat[i],
                                # insourcing acceptance
                                [- D.accept[D.S_effective[w, i], r] - tau_val[r] for r in intersect(D.R_accept[w, i], R_swap_x)],
                                [- D.accept[D.S_effective[w, i], r] for r in setdiff(D.R_accept[w, i], R_swap_x)]
                            )
                        end
                    end
                    if !isnothing(gamma_val) & swapping
                        for i in setdiff(D.inswappable[w], S_omit)
                            sigma_hat[i] = vcat(sigma_hat[i],
                                [- D.accept[D.S_effective[w, i], r] + gamma_val[r] for r in intersect(D.R_accept[w, i], D.R_w_DNS[w])]
                            )
                        end
                    end
                    sigma_hat = Dict(i => maximum(thelist) for (i, thelist) in sigma_hat)

                    # construct cut
                    if isnothing(gamma_val)
                        JuMP.@constraint(MP, theta >= - sum(sigma_val[i] * x[i] for i in intersect(S_omit, D.S_w[w])) - sum(sigma_hat[i] * x[i] for i in setdiff(D.S_w[w], S_omit)))
                    else
                        JuMP.@constraint(MP, theta >= sum(gamma_val[r] for r in D.R_w_DNS[w]) - sum(sigma_val[i] * x[i] for i in intersect(S_omit, D.S_w[w])) - sum(sigma_hat[i] * x[i] for i in setdiff(D.S_w[w], S_omit)))
                    end
                end
            end
            if !(activated | pargs.activated) | temp_deactivate
                # construct cut from full subproblem solution
                if length(D.R_w_DNS[w]) > 0
                    JuMP.@constraint(MP, theta >= sum(-sigma_val[i] * x[i] for i in D.S_w[w]) + sum(gamma_val[r] for r in D.R_w_DNS[w]))
                else
                    JuMP.@constraint(MP, theta >= sum(-sigma_val[i] * x[i] for i in D.S_w[w]))
                end
            end

            # check time limit
            if ((total_time >= total_time_limit - MP_time) & speedy) | ((total_time >= total_time_limit) & !speedy)
                @info "Out of time to solve subproblems, stopped at scenario" * w * " at iteration " * string(iteration) * "\n"
                break
            end
        end

        #--- subproblem solve times
        max_time = isempty(SP_time_all) ? 0 : maximum(SP_time_all)
        sum_time = isempty(SP_time_all) ? 0 : sum(SP_time_all)
        push!(SP_max_times, max_time)
        push!(SP_times, sum_time)

        #--- upper bounds and new incumbent
        valid_ub = isempty(setdiff(D.O, keys(obj_SP))) # solved all SPs
        ss_temp = isempty(obj_SP) ? 0 : sum(D.p[w] * obj for (w, obj) in obj_SP)
        ub = sum(D.l[i] * x_val[i] for i in D.S) + ss_temp
        prev_ub = isempty(upper_bounds) ? LARGE_NUMBER : upper_bounds[end]
        if valid_ub
            ub = isempty(upper_bounds) ? ub : min(upper_bounds[end], ub)
            # update incumbent solution
            if (ub < prev_ub)
                x_incumbent = deepcopy(x_val)
                S_x, S_xr = get_roster(D, x_incumbent)
                _, mip_obj = SecondStageSolutions(D, S_xr, swapping)
                mb = sum(D.l[i] * x_incumbent[i] for i in D.S) + mip_obj
            else
                mb = mip_bounds[end]
            end
        else
            # invalid upper bound - did not solve all subproblems, either by timeout or due to incorrect solve
            ub = isempty(upper_bounds) ? LARGE_NUMBER : upper_bounds[end]
            mb = isempty(mip_bounds) ? LARGE_NUMBER : mip_bounds[end]
        end
        push!(second_stage, ss_temp)
        push!(upper_bounds, ub)
        push!(mip_bounds, mb)

        @printf("MILP: %.2f - MIP: %.2f - Bound: %.2f\n", upper_bounds[end], mip_bounds[end], lower_bounds[end])

        #--- check for convergence or timeout
        if (speedy & (total_time >= total_time_limit - MP_time)) | (!speedy & (total_time >= total_time_limit)) | (iteration >= max_iterations)
            break
        elseif (abs((upper_bounds[end] - lower_bounds[end]) / lower_bounds[end]) < benders_gap) | (lower_bounds[end] > upper_bounds[end])
            if speedy
                # converged speedily, lower the tolerance
                speedy = false
                JuMP.set_optimizer_attribute(MP, "MIPGap", mip_gap)
            else
                # regular convergence
                break
            end
        elseif (iteration >= max_speedy_iterations) & speedy
            speedy = false
            JuMP.set_optimizer_attribute(MP, "MIPGap", mip_gap)
        end

        #--- LPO - compute new activated core point
        if pargs.dynamic_n > 0
            push!(S_act, S_x)
            _, x0, cp_solvetime, S_core = act_core_point_dynamic(D, S_act, pargs.dynamic_n)
            if cp_solvetime < 0
                # switch to no accel if no act core point
                pargs = NULLPARETO
            else
                MP_times[end] += cp_solvetime
                total_time += cp_solvetime
            end
        end
    end

    # obtain true lower bound
    if speedy
        JuMP.set_optimizer_attribute(MP, "MIPGap", mip_gap)
        MP_time, lb, _ = solve_MP(MP)
        push!(MP_times, MP_time)
        push!(lower_bounds, lb)

        # dummy column values
        push!(upper_bounds, upper_bounds[end])
        push!(mip_bounds, mip_bounds[end])
        push!(SP_max_times, SP_max_times[end])
        push!(SP_times, SP_times[end])
        push!(second_stage, second_stage[end])
        push!(is_speedy, false)
    end

    # obtain the final subproblem solutions
    if isnothing(x_incumbent)
        vals = nothing
    else
        _, S_xr = get_roster(D, x_incumbent)
        vals, _ = SecondStageSolutions(D, S_xr, swapping)
    end

    return lower_bounds, upper_bounds, mip_bounds, MP_times, SP_max_times, SP_times, x_incumbent, second_stage[end], vals, is_speedy
end
