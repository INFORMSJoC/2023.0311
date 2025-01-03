MIN_PER_HOUR = 60

"""
Retrieve data path
"""
function get_data_path()
    filepath = joinpath(dirname(@__DIR__), "..", "data-path.txt")
    # Read first line.
    if isfile(filepath)
        open(filepath) do f
            return strip(readline(f))
        end
    else
        @warn "Please write path to data directory as the only line in a txt file called `data-path.txt` and store in the root directory."
    end
end

function record_xval(
        x_val,
        D::DayBefore,
        fp_stem::String
    )
    Ix = [i for i in D.S if x_val[i] > 0.5]

    # generate x
    trips = [r for i in Ix for r in D.S_list[i]]
    seq = [s for i in Ix for s in 1:length(D.S_list[i])]
    drivers = [i for i in Ix for _ in 1:length(D.S_list[i])]
    driver_id = [i for i in 1:length(Ix) for _ in 1:length(D.S_list[Ix[i]])]

    # build and write x
    xdat = DataFrame(trip_id = trips,
                    seq = seq,
                    driver_id = drivers,
                    driver_num = driver_id)
    CSV.write(fp_stem * "_x.csv", xdat)
    xval = unique(xdat[!,:driver_id])
    return xval
end

"""
Record itinerary solution as a CSV

### Keywords
* `x_val` - optimal first-stage solution
* `vals` - SP vals - y, z, v, a
* `D` - DayBefore
* `fp` - stem of each filepath to output CSVs
### Writes
* CSV with features (driver ID, sequence ID, trip ID)
"""
function record_sol(
        x_val,
        vals,
        D::DayBefore,
        fp_stem::String
    )
    Ix = [i for i in D.S if x_val[i] > 0.5]
    return record_sol(Ix, vals, D, fp_stem)
end

"""
Record itinerary solution as a CSV

### Keywords
* `Ix` - set of optimal itineraries
* `vals` - SP vals - y, z, v, a
* `D` - DayBefore
* `fp` - stem of each filepath to output CSVs
### Writes
* CSV with features (driver ID, sequence ID, trip ID)
"""
function record_sol(
        Ix::Vector{String},
        vals,
        D::DayBefore,
        fp_stem::String
    )
    y_val = Dict(); z_val = Dict(); a_val = Dict()
    for (w, val) in vals
        y, z, v, a = val
        y_val[w] = y
        z_val[w] = z
        if !isnothing(a)
            a_val[w] = a
        end
    end

    # build and write y
    chosen_outsourced = [(w, r) for w in D.O for r in D.R_w[w] if y_val[w][r] > 0.5]
    ydat = DataFrame(chosen_outsourced)
    if nrow(ydat) > 0
        rename!(ydat, ["scenario", "trip_id"])
        CSV.write(fp_stem * "_y.csv", ydat)
    end

    # build and write z
    chosen_connections = [(w, r, s) for w in D.O for (r, s) in D.A_w[w] if z_val[w][(r, s)] > 0.5]
    zdat = DataFrame(chosen_connections)
    if nrow(zdat) > 0
        rename!(zdat, ["scenario", "trip_origin", "trip_dest"])
        CSV.write(fp_stem * "_z.csv", zdat)
    end

    if !isempty(a_val)
        # write insourcing decisions
        swappable = [w for w in D.O if w in keys(a_val)]
        swapped = [(w, i, r) for w in swappable for i in intersect(Ix, D.inswappable[w]) for r in D.R_accept[w, i] if a_val[w][i, r] > 0.5]

        adat = DataFrame(swapped)
        if nrow(adat) > 0
            rename!(adat, ["scenario", "path_id", "trip_id"])
            CSV.write(fp_stem * "_a.csv", adat)
        end
    end

    return nothing
end

"""
Format first-stage solution into roster

### Keywords
* `D` - DayBefore
* `x_val` - first-stage solution {0,1}
### Returns
* list of itineraries
* dict of each request to its itinerary
"""
function get_roster(
        D::DayBefore,
        x_val
    )
    I_x = [i for i in D.S if x_val[i] > 0.5]
    I_xr = Dict(r => first([i for i in I_x if r in D.S_list[i]]) for r in D.R)
    return I_x, I_xr
end

"""
Compute:
* minutes online
* minutes deadheading
* minutes in service of a request
* minutes idling
* shift lengths
* start times and end times of each shift
* start times of each request

### Keywords
* `params` - outputs from `times` function in `params.jl`
* `Ix` - roster (list of path IDs)
"""
function service_metrics(
        params,
        Ix::Vector{String}
    )
    tp, td, tt, tt_d = params

    m = Dict()
    min_online, min_deadhead, min_service = 0, 0, 0
    shift_starttime, shift_endtime, trip_starttime = Dict(), Dict(), Dict()
    for path_id in Ix
        I = Array{String}(setdiff(split(path_id, "-"), ["src", "sink"]))
        if length(I) == 1
            shift_starttime[path_id] = tp[I[1]] - tt_d[("src", I[1])]
            shift_endtime[path_id] = tp[I[1]] + tt[I[1]] + tt_d[(I[1], "sink")]
            num_delays = 0

            min_online += shift_endtime[path_id] - shift_starttime[path_id]
            min_deadhead += tt_d[("src", I[1])] + tt_d[(I[1], "sink")]
            min_service += tt[I[1]]

            trip_starttime[I[1]] = tp[I[1]]
        else
            m[path_id] = itin_cost_model(I, tt, tt_d, tp, PENALTY, LABOR_MIN, WINDOW, BIG_M)
            JuMP.optimize!(m[path_id])

            shift_starttime[path_id] = value(m[path_id][:start_time])
            shift_endtime[path_id] = value(m[path_id][:end_time])
            min_online += shift_endtime[path_id] - shift_starttime[path_id]
            min_service += sum(tt[r] for r in I)
            min_deadhead += sum(tt_d[(I[i], I[i+1])] for i in 1:(length(I)-1))
            for r in I
                trip_starttime[r] = value(m[path_id][:pickup_time][r])
            end
        end
    end
    min_idle = min_online - min_service - min_deadhead
    return min_online, min_deadhead, min_service, min_idle, shift_starttime, shift_endtime, trip_starttime
end

"""
Compute global, shift, and trip metrics given a roster.

### Keywords
* `params` - outputs from 'times' function in `params.jl`
* `Ix` - list of itineraries
* `trial_id` - trial identifier
### Returns
* global metrics DF (minutes online, minutes deadheading, minutes in-service, minutes idle)
* shift metrics DF (path ID, shift start, shift end)
* trip metrics DF (trip ID, trip start, trip end, trip request time)
"""
function all_metrics(
        params,
        Ix::Vector{String},
        trial_id::String
    )
    # unravel parameters
    tp, td, tt, tt_d = params

    # instantiate dataframes
    global_metrics = DataFrame(trial=[], min_online=[], min_deadhead=[], min_service=[], min_idle=[])
    shift_metrics = DataFrame(trial=[], path_id=[], shift_start=[], shift_end=[])
    trip_metrics = DataFrame(trial=[], trip_id=[], trip_start=[], trip_request=[], trip_end = [])

    # get metrics
     sm = service_metrics(params, Ix)
     min_online, min_deadhead, min_service, min_idle, shift_starttime, shift_endtime, trip_starttime = sm

     # record the metrics
     push!(global_metrics, [trial_id, min_online, min_deadhead, min_service, min_idle])
     for path_id in Ix
         push!(shift_metrics, [trial_id, path_id, shift_starttime[path_id], shift_endtime[path_id]])
     end
     for r in keys(tp)
         push!(trip_metrics, [trial_id, r, trip_starttime[r], tp[r], trip_starttime[r] + tt[r]])
     end
    return global_metrics, shift_metrics, trip_metrics
end

"""
Obtain total adhoc hours and adhoc itineraries given recourse decision.

### Keywords
* `D` - DayBefore
* `w` - scenario ID
* `ss_vals` - second-stage outputs
* `params` - (tp, td, tt, tt_d)
### Returns
* total adhoc hours
* adhoc itineraries (DataFrame)
"""
function adhoc_hours(
        D::DayBefore,
        w::String,
        ss_vals,
        params
    )
    # instantiate
    y_val, z_val, _, _ = ss_vals
    tp, td, tt, _ = params
    chours = 0
    adhoc_itin = Vector{String}[]
    outsourced = [r for r in D.R_w[w] if y_val[r] > 0.5]
    connections = [(r, s) for (r, s) in D.A_w[w] if z_val[(r, s)] > 0.5]

    # build adhoc itineraries
    beginnings = [r for r in outsourced if isempty([(q, s) for (q, s) in connections if s == r])]
    @assert length(beginnings) == sum(y_val) - sum(z_val)
    for r in beginnings
        itin = String[]
        next = [r]
        while !isempty(next)
            cur = first(next) # there is only one
            push!(itin, cur)
            next = [s for (q, s) in connections if q == cur]
        end
        push!(adhoc_itin, itin)
    end

    # compute total time online and format itin into DataFrame
    trip_ids, seq, drivers = [], [], []
    driver_id = 0
    for itin in adhoc_itin
        # hours
        end_shift = maximum([td[r] for r in itin])
        beginning_shift = minimum([tp[r] for r in itin])
        chours += (end_shift - beginning_shift) / MIN_PER_HOUR

        # build dataframe
        driver_id += 1
        for (sequence, trip) in enumerate(itin)
            push!(trip_ids, trip)
            push!(seq, sequence)
            push!(drivers, driver_id)
        end
    end
    adhoc_dat = DataFrame(
        trip_id = trip_ids,
        seq = seq,
        driver_id = drivers,
        scenario_id = [w for _ in 1:length(trip_ids)]
    )

    return chours, adhoc_dat
end

"""
Obtain modified itineraries given recourse decisions, cancellations, and DNS. 
Compute number of total delays and total planned labor.

### Keywords
* `D` - DayBefore
* `w` - scenario ID
* `I_x` - 1st stage itin
* `I_xr` - Dict (r => itin to which r belongs)
* `ss_vals` - second-stage decisions
* `itin_attr_dict` - dict of itinerary attributes
* `add_trip` - added trips DataFrame
* `remove_trip` - removed trips DataFrame
### Returns
* number of delays
* planned labor total in minutes
"""
function recourse_itin(
        D::DayBefore,
        w::String,
        I_x::Vector{String},
        I_xr::Dict{String,String},
        ss_vals,
        itin_attr_dict::Dict{String,Tuple{Float64,Int64}},
        add_trip::DataFrame,
        remove_trip::DataFrame
    )
    y_val, z_val, v_val, a_val = ss_vals

    #--- count original delays and shift length
    labor_min = sum(itin_attr_dict[D.S_effective[w, i]][1] for i in I_x if (w, i) in keys(D.S_effective))
    delays = sum(itin_attr_dict[D.S_effective[w, i]][2] for i in I_x if (w, i) in keys(D.S_effective))

    #--- obtain modifications to each itinerary
    outsourced = Dict(I_xr[r] => r for r in setdiff(D.R_w[w], D.R_w_DNS[w]) if y_val[r] > 0.5)
    accepted = Dict(); removed = Dict()
    if !isnothing(a_val)
        R_swap = setdiff(unique([r for ((v, i), R_a) in D.R_accept for r in R_a if (v == w) & (i in I_x)]), D.R_w_DNS[w])
        accepted = Dict(i => r for i in intersect(D.inswappable[w], I_x) for r in D.R_accept[w, i] if a_val[i,r] > 0.5)
        removed = Dict(I_xr[r] => r for r in R_swap if v_val[r] > 0.5)
    end

    #--- modify each itinerary
    modified_itin = Dict()
    for i in vcat([a for a in keys(accepted)], [o for o in keys(outsourced)], [r for r in keys(removed)])
        path = D.S_list[i]

        # itin becomes empty
        if isempty(setdiff(path, vcat(setdiff(D.R, D.R_w[w]), D.R_w_DNS[w])))
            modified_itin[i] = []
            continue
        end

        # remove cancellations
        path = [path[j] for j in 1:length(path) if path[j] in D.R_w[w]]

        # replace dns with dns label
        new_path = [path[j] in D.R_w_DNS[w] ? DNS_label : path[j] for j in 1:length(path)]
        while new_path[1] == DNS_label
            # get rid of head dns
            new_path = (new_path[1] == DNS_label) ? new_path[2:end] : new_path
        end
        while new_path[end] == DNS_label
            # get rid of tail dns
            new_path = (new_path[end] == DNS_label) ? new_path[1:(end-1)] : new_path
        end
        new_new_path = [new_path[1]]
        for j in 2:length(new_path)
            # get rid of duplicate middle dns
            if !((new_path[j] == DNS_label) & (new_new_path[end] == DNS_label))
                push!(new_new_path, new_path[j])
            end
        end
        path = new_new_path

        modified_itin[i] = path
    end

    #--- count changes to delay due to recourse
    for (i, r) in accepted
        path = modified_itin[i]
        modified_path_id = join(vcat(["src"], modified_itin[i], ["sink"]), "-")
        add_info = eachrow(filter(row -> (row[:effective_path_id] == modified_path_id) & (row[:trip_id] == r), add_trip))[1]
        delays += add_info[:delta_delays]
        labor_min += add_info[:delta_shiftlength]
        #print(add_info[:delta_delays])
    end

    for (i, r) in merge(outsourced, removed)
        path = modified_itin[i]
        if !isempty(path)
            modified_path_id = join(vcat(["src"], modified_itin[i], ["sink"]), "-")
            remove_info = eachrow(filter(row -> (row[:effective_path_id] == modified_path_id) & (row[:trip_id] == r), remove_trip))[1]
            delays -= remove_info[:delta_delays]
            labor_min -= remove_info[:delta_shiftlength]
        end
    end

    return delays, labor_min
end
