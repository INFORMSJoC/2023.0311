#########################################################
################ FULL 2-STAGE MODEL
#########################################################


"""
Build full two-stage model.

### Keywords
* `D` - DayBefore
* `swapping` - whether recourse model is swapping model
* `relaxed` - whether second stage variables are relaxed
* `time_limit_sec` - time limit for Gurobi optimizer in seconds
* `mip_gap` - optimality tolerance for solver
### Returns
* JuMP Model
"""
function direct_model(
        D::DayBefore,
        swapping::Bool=false,
        relaxed::Bool=true,
        time_limit::Int64=TIME_LIMIT,
        mip_gap::Float64=GAP
    )
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB), "TimeLimit" => time_limit, "MIPGap" => mip_gap))

    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB)))

    JuMP.@variable(m, x[D.S], Bin)
    if relaxed
        JuMP.@variable(m, y[w in D.O, i in D.S_w[w], r in intersect(D.S_list[i], D.R_w[w])] >= 0)
        JuMP.@variable(m, z[w in D.O, D.A_w[w]] >= 0)
    else
        JuMP.@variable(m, y[w in D.O, i in D.S_w[w], r in intersect(D.S_list[i], D.R_w[w])], Bin)
        JuMP.@variable(m, z[w in D.O, D.A_w[w]], Bin)
    end

    JuMP.@constraints m begin
        # first stage
        [r in D.R], sum(x[i] for i in D.S_r[r]) == 1
        sum(x) <= D.W
        # second stage matching
        [w in D.O, r in D.R_w[w]], sum(z[w, e] for e in intersect(D.d_f[r], D.A_w[w])) <= sum(y[w, i, r] for i in intersect(D.S_r[r], D.S_w[w]))
        [w in D.O, r in D.R_w[w]], sum(z[w, e] for e in intersect(D.A_w[w], D.d_b[r])) <= sum(y[w, i, r] for i in intersect(D.S_r[r], D.S_w[w]))
    end

    if swapping
        O_swap = [w for w in D.O if length(D.inswappable[w]) > 0]
        O_inswap_only = [w for w in O_swap if length(setdiff(D.inswappable[w], D.outswappable[w])) > 0]
        O_outswap_only = [w for w in O_swap if length(setdiff(D.outswappable[w], D.inswappable[w])) > 0]
        O_both = [w for w in O_swap if length(intersect(D.outswappable[w], D.inswappable[w])) > 0]
        O_neither = [w for w in O_swap if length(D.notswappable[w]) > 0] # itineraries not participating in insourcing within insourcing scenarios

        # introduce new variables accordingly
        if relaxed & (length(O_swap) > 0)
            JuMP.@variables m begin
                v[w in O_swap, i in D.outswappable[w], r in intersect(D.R_swap[w], D.S_list[i])] >= 0
                a[w in O_swap, i in D.inswappable[w], r in D.R_accept[w, i]] >= 0
            end
        elseif length(O_swap) > 0
            JuMP.@variables m begin
                v[w in O_swap, i in D.outswappable[w], r in intersect(D.R_swap[w], D.S_list[i])], Bin
                a[w in O_swap, i in D.inswappable[w], r in D.R_accept[w, i]], Bin
            end
        end

        if length(O_inswap_only) > 0
            # one modification
            JuMP.@constraint(m, inswap_limit[w in O_inswap_only, i in setdiff(D.inswappable[w], D.outswappable[w])],
                sum(y[w, i, r] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) +
                sum(a[w, i, r] for r in D.R_accept[w, i]) <= x[i]
            )
        end
        if length(O_outswap_only) > 0
            # one modification
            JuMP.@constraint(m, outswap_limit[w in O_outswap_only, i in setdiff(D.outswappable[w], D.inswappable[w])],
                sum(y[w, i, r] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) +
                sum(v[w, i, r] for r in intersect(D.S_list[i], D.R_swap[w])) <= x[i]
            )
        end

        if length(O_both) > 0
            # one modification
            JuMP.@constraint(m, both_limit[w in O_both, i in intersect(D.outswappable[w], D.inswappable[w])],
                sum(y[w, i, r] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) +
                sum(v[w, i, r] for r in intersect(D.S_list[i], D.R_swap[w])) +
                sum(a[w, i, r] for r in D.R_accept[w, i]) <= x[i]
            )
        end
        if length(O_neither) > 0
            # one modification
            JuMP.@constraint(m, neither_limit[w in O_neither, i in D.notswappable[w]],
                sum(y[w, i, r] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) <= x[i]
            )
        end

        # scenarios without insourcing
        if length(setdiff(D.O, O_swap)) > 0
            JuMP.@constraint(m, noinsource_limit[w in setdiff(D.O, O_swap), i in D.S_w[w]],
                sum(y[w, i, r] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) <= x[i]
            )
        end


        # handle dns without acceptance option
        JuMP.@constraint(m, dns_noaccept[w in [v for v in D.O if !isempty(D.R_w_DNS[v])], r in setdiff(D.R_w_DNS[w], [r for i in D.inswappable[w] for r in D.R_accept[w, i]])],
            sum(y[w, i, r] for i in intersect(D.S_w[w], D.S_r[r])) == 1
        )

        # swapping flow balance
        if length(O_swap) > 0
            JuMP.@constraints m begin
                swap_flowbal[w in O_swap, r in D.R_swap[w]], sum(v[w, i, r] for i in intersect(D.S_r[r], D.S_w[w])) == sum(a[w, i, r] for i in D.inswappable[w] if r in D.R_accept[w, i])
                dns[w in [v for v in D.O if !isempty(D.R_w_DNS[v])], r in intersect([r for i in D.inswappable[w] for r in D.R_accept[w, i]], D.R_w_DNS[w])], sum(y[w, i, r] for i in intersect(D.S_w[w], D.S_r[r])) + sum(a[w, i, r] for i in D.inswappable[w] if r in D.R_accept[w, i]) == 1
            end
            JuMP.@expression(m, second_stage, sum(
                # scenarios with insourcing
                D.p[w] * (
                    # casting out requests
                    - sum(D.remove[D.S_effective[w, i], r] * v[w, i, r] for r in D.R_swap[w] for i in intersect(D.S_r[r], D.S_w[w]))
                )
                for w in O_swap
            ) + sum(
                # accepting insourced requests
                D.p[w] * D.accept[D.S_effective[w, i], r] * a[w, i, r] for w in O_swap for i in D.inswappable[w] for r in D.R_accept[w, i]
            ) + sum(
                # outsourcing to ad hocs - all scenarios
                D.p[w] * (
                    sum((D.f + D.v[r]) * y[w, i, r] for i in D.S_w[w] for r in intersect(D.S_list[i], D.R_w[w])) -
                    sum(D.remove[D.S_effective[w, i], r] * y[w, i, r] for i in D.S_w[w] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) +
                    sum((D.c[e] - D.f) * z[w, e] for e in D.A_w[w])
                )
                for w in D.O
            ))
        else
            # it's just outsourcing even though we wanted to swap
            JuMP.@expression(m, second_stage,
                sum(D.p[w] * (
                    sum((D.f + D.v[r]) * y[w, i, r] for i in D.S_w[w] for r in intersect(D.S_list[i], D.R_w[w])) -
                    sum(D.remove[D.S_effective[w, i], r] * y[w, i, r] for i in D.S_w[w] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) +
                    sum((D.c[e] - D.f) * z[w, e] for e in D.A_w[w])
                )
                for w in D.O
            ))
        end
    else
        JuMP.@constraints m begin
            # second stage limit
            [w in D.O, i in D.S_w[w]], sum(y[w, i, r] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) <= x[i]
        end

        JuMP.@expression(m, second_stage, sum(D.p[w] * (
                sum((D.f + D.v[r]) * y[w, i, r] for i in D.S_w[w] for r in intersect(D.S_list[i], D.R_w[w])) -
                sum(D.remove[D.S_effective[w, i], r] * y[w, i, r] for i in D.S_w[w] for r in setdiff(intersect(D.S_list[i], D.R_w[w]), D.R_w_DNS[w])) +
                sum((D.c[e] - D.f) * z[w, e] for e in D.A_w[w])
            )
            for w in D.O
        ))
    end
    JuMP.@objective(m, Min, sum(D.l[i] * x[i] for i in D.S) + second_stage)

    return m
end


"""
Append objective bounds callback function to model

### Keywords
* `D` - DayBefore
* `m` - JuMP MODEL
* `log_stem` - stem to location to output log
### Returns
* JuMP Model with callback function
"""
function append_callback(
        D::DayBefore,
        m::JuMP.Model,
        log_stem::String
    )
    log_filepath = log_stem * "_bounds.csv"

    #--- write initial line of log
    df = DataFrame(runtime=[], objbnd=[], objbst=[], call_type=[])
    CSV.write(log_filepath, df)

    ###--- FIND BOUNDS
    first_cb = true
    function cb_bounds(cb_data, cb_where::Cint)
        if cb_where != Gurobi.GRB_CB_MIPSOL
            return
        end

        # get attributes
        runtime, objbnd, objbst = Ref{Cdouble}(), Ref{Cdouble}(), Ref{Cdouble}()
        Gurobi.GRBcbget(cb_data, cb_where, Gurobi.GRB_CB_RUNTIME, runtime)
        Gurobi.GRBcbget(cb_data, cb_where, Gurobi.GRB_CB_MIPSOL_OBJBND, objbnd)
        Gurobi.GRBcbget(cb_data, cb_where, Gurobi.GRB_CB_MIPSOL_OBJBST, objbst)
        objbnd = first_cb ? missing : objbnd

        # get incumbent solution
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        x_val = callback_value.(cb_data, m[:x])
        incumbent_stem = log_stem * string(runtime)
        record_xval(x_val, D, incumbent_stem)

        # write
        df = DataFrame(runtime=runtime, objbnd=objbnd, objbst=objbst, call_type=cb_where)
        CSV.write(log_filepath, df, append=true)
        first_cb = false
        return
    end
    MOI.set(m, Gurobi.CallbackFunction(), cb_bounds)

    return m
end

"""
Build direct two-stage model with bounds callback fcn

### Keywords
* `log_stem` - front stem to logging output filepath for callback fcn
* `D` - DayBefore
* `swapping` - whether to implement swapping recourse
* `relaxed` - whether to relax second stage decisions
* `time_limit_sec` - time limit for Gurobi optimizer in seconds
* `mip_gap` - optimality tolerance for solver
### Returns
* JuMP Model with callback loaded in
"""
function direct_model_callback(
        log_stem::String,
        D::DayBefore,
        swapping::Bool,
        relaxed::Bool,
        time_limit_sec::Int64=TIME_LIMIT,
        mip_gap::Float64=GAP
    )
    m = direct_model(D, swapping, relaxed, time_limit_sec, mip_gap)
    return append_callback(D, m, log_stem)
end

"""
Build and solve two-stage model.
``
### Keywords
* `D` - DayBefore
* `swapping` - whether to implement swapping recourse
* `relaxed` - whether to relax second stage decisions
* `callback` - whether to implement callback
* `log_stem` - if callback, write callback outputs to file with this stem
### Returns
* time to solve in seconds
* objective value
* optimality gap
* x values
* second stage expected costs
"""
function direct_model_solve(
        D::DayBefore,
        swapping::Bool,
        relaxed::Bool,
        callback::Bool=false,
        log_stem::String="cb"
    )
    m = callback ? direct_model_callback(log_stem, D, swapping, relaxed) : direct_model(D, swapping, relaxed)
    solvetime = @elapsed optimize!(m)
    feas = JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    obj = feas ? JuMP.objective_value(m) : missing
    bnd = feas ? JuMP.objective_bound(m) : missing
    gap = feas ? (obj - bnd) / bnd : missing
    ss = feas ? JuMP.value(m[:second_stage]) : missing
    x_val = feas ? JuMP.value.(m[:x]) : missing
    if callback
        df = DataFrame(runtime=solvetime, objbnd=bnd, objbst=obj, call_type="")
        CSV.write(log_stem * "_bounds.csv", df, append=true)
    end
    return solvetime, obj, gap, x_val, ss
end

#########################################################
################ OOS MODEL
#########################################################

"""
Given first-stage decision and cancelled trips, build OOS second-stage model (activated)

### Keywords
* `D` - DayBefore
* `I_x` - List of 1st-stage itineraries
* `I_xr` - Dict (request => home itinerary of request)
* `cancelled` - list of cancelled trip IDs
* `R_w_dns` - list of DNS trip IDs
* `time_limit` - time limit in seconds
* `mip_gap` - integrality gap tolerance
### Returns
* JuMP Model of OOS recourse
"""
function oos_scenario(
        D::DayBefore,
        I_x::Vector{String},
        I_xr::Dict{String,String},
        cancelled::Vector{String},
        R_w_dns::Vector{String},
        time_limit::Int64=TIME_LIMIT,
        mip_gap::Float64=MIP_GAP
    )
    # oos parameters
    R_w = setdiff(D.R, cancelled)
    A_w = [(a, b) for (a, b) in D.A if (a in R_w) & (b in R_w)]
    I_w = [i for i in I_x if length(intersect(D.S_list[i], R_w)) > 0]
    I_effective = Dict()
    for path_id in I_w
        new_path = [path[i] in R_w_DNS ? DNS_label : path[i] for i in 1:length(path) if path[i] in R_w[w]]
        if isempty(setdiff(new_path, [DNS_label]))
            continue
        end

        # get rid of head/tail dns
        while new_path[1] == DNS_label
            new_path = new_path[2:end]
        end
        while new_path[end] == DNS_label
            new_path = new_path[1:(end-1)]
        end

        # remove DNS label duplicates
        new_new_path = []
        for i in 1:length(new_path)
            if i == 1
                push!(new_new_path, new_path[i])
            elseif !((new_path[i] == DNS_label) & (new_new_path[end] == DNS_label))
                push!(new_new_path, new_path[i])
            end
        end
        new_path = new_new_path

        # contribute to cost
        if length(new_path) > 0
            I_effective[path_id] = join(vcat("src", new_path, "sink"), '-')
        end
    end
    R_accept = Dict(
        path_id =>
            intersect(unique([r for w in D.O if (w, path_id) in keys(D.R_accept) for r in D.R_accept[w, path_id]]), R_w)
        for path_id in keys(I_effective)
    )

    # swapping sets
    R_swap = setdiff([r for i in I_w for r in R_accept[i]], R_w_DNS)
    inswappable = [i for i in I_w if length(R_accept[i]) > 0]
    swappable = !isempty([r for i in I_w for s in R_accept[i]])
    outswappable = [i for i in I_w if length(intersect(R_swap, D.S_list[I])) > 0]
    notswappable = setdiff(I_w, union(inswappable, outswappable))

    # build model
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GRB)))
    JuMP.@variables m begin
        y[R_w], Bin
        z[A_w], Bin
    end

    # outsourcing to ad hocs
    JuMP.@constraints m begin
        adhoc_out[r in R_w], sum(z[(r, s)] for (q, s) in A_w if q == r) <= y[r]
        adhoc_in[r in R_w], sum(z[(q, r)] for (q, s) in A_w if s == r) <= y[r]
    end

    if swappable
        JuMP.@variables m begin
            v[r in R_swap], Bin
            a[I in inswappable, r in R_accept[I]], Bin
        end
    end

    # single modifications
    for I in setdiff(inswappable, outswappable)
        JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[I], R_w), R_w_DNS)) +
        sum(a[I, r] for r in R_accept[I]) <= 1)
    end

    for I in setdiff(outswappable, inswappable)
        JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[I], R_w), R_w_DNS)) +
        sum(v[r] for r in intersect(D.S_list[I], R_swap)) <= 1)
    end

    for I in intersect(inswappable, outswappable)
        JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[I], R_w), R_w_DNS)) +
        sum(v[r] for r in intersect(D.S_list[I], R_swap)) +
        sum(a[I, r] for r in R_accept[I]) <= 1)
    end

    for I in notswappable
        JuMP.@constraint(m, sum(y[r] for r in setdiff(intersect(D.S_list[I], R_w), R_w_DNS)) <= 1)
    end

    if swappable
        # flow balance
        JuMP.@constraints m begin
            swap_flowbal[r in R_swap], v[r] == sum(a[I, r] for I in inswappable if r in R_accept[I])
        end

        JuMP.@objective(m, Min,
        - sum(
            # remove insourced requests
            D.remove[I_effective[I_xr[r]], r] * v[r] for r in R_swap
        ) + sum(
            # accept insourced requests
            D.accept[I_effective[I], r] * a[I, r] for I in inswappable for r in R_accept[I]
        ) +
            # remove outsourced requests
            sum((D.f + D.v[r]) * y[r] for r in R_w) -
            sum(D.remove[I_effective[I_xr[r]], r] * y[r] for r in setdiff(R_w, R_w_DNS)) +
            # construct ad hoc itineraries
            sum((D.c[e] - D.f) * z[e] for e in A_w)
        )
    else
        # it's just outsourcing even though we wanted to swap
        JuMP.@objective(m, Min,
            sum((D.f + D.v[r]) * y[r] for r in R_w) -
            sum(D.remove[I_effective[I_xr[r]], r] * y[r] for r in setdiff(R_w, R_w_DNS)) +
            sum((D.c[e] - D.f) * z[e] for e in A_w)
        )
    end

    return m
end

"""
Solve out-of-sample model

### Keywords
* `m` - OOS JuMP Model
### Returns
* solvetime in seconds
* objective value
* second stage values
"""
function solve_oos_model(
        m::JuMP.Model
    )
    st = @elapsed JuMP.optimize!(m)
    local obj, ss_vals
    if JuMP.termination_status(m) == MOI.OPTIMAL
        obj = JuMP.objective_value(m)
        y_val = JuMP.value.(m[:y])
        z_val = JuMP.value.(m[:z])
        a_val = :a in keys(m.obj_dict) ? JuMP.value.(m[:a]) : nothing
        v_val = :v in keys(m.obj_dict) ? JuMP.value.(m[:v]) : nothing
        ss_vals = y_val, z_val, a_val, v_val
    else
        ss_vals = nothing
    end
    return st, obj, ss_vals
end
