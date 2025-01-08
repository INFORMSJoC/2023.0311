#########################################################
################ DATA STRUCTURE FOR MODEL INPUTS
#########################################################

#--- constants
REQUIRED_FILES = ["trips", "scenarios", "edges", "itin", "itin_attributes", "remove_trip", "add_trip"] .* ".csv"
const SL_PENALTY = 1000
const NOCANCEL_FLAG = "all"

"""
Data Structure for two-stage day-before model

### Sets
* `R` - requests
* `A` - VSN arcs (origin request, dest request)
* `S` - itinerary SDs
* `O` - scenario IDs

### Index sets
* `S_list` - itinerary ID -> ordered list of requests
* `S_effective` - (scenario ID, itinerary ID) -> effective itinerary ID (post-cancellation)
* `d_f` - request -> set of outgoing edge tuples
* `d_b` - request -> set of incoming edge tuples
* `S_r` - request -> set of itinerary IDs containing request
* `R_w` - scenario -> requests not cancelled
* `A_w` - scenario -> edge tuples among requests not cancelled
* `S_w` - scenario -> set of itinerary IDs with at least one request not cancelled
* `R_accept` - (scenario, itinerary) -> set of requests that can be accepted into itinerary in scenario

### Parameters
* `c` - edge tuple -> edge cost
* `p` - scenario -> probability
* `l` - itinerary -> expected shift cost over scenarios
* `l_nr` - itinerary -> shift cost
* `l_scen` - (itinerary ID, scenario ID) -> itinerary cost in scenario
* `remove` - (effective itinerary ID, request) -> marginal cost savings to remove request from itinerary
* `accept` - (effective itinerary ID, request) -> marginal cost to add request to itinerary
* `v` - request -> variable cost for a adhoc to serve request r
* `W` - number of workers
* `f` - fixed cost to hire adhoc

### Swapping-related index sets
* `insourcing`
* `inswappable`
* `outswappable`
* `notswappable`
* `R_swap`
"""
mutable struct DayBefore
    R::Vector{String}
    A::Vector{Tuple{String,String}}
    S::Vector{String}
    O::Vector{String}

    S_list::Dict{String,Vector{String}}
    S_effective::Dict{Tuple{String,String},String}
    d_f::Dict{String,Vector{Tuple{String,String}}}
    d_b::Dict{String,Vector{Tuple{String,String}}}
    S_r::Dict{String,Vector{String}}
    R_w::Dict{String,Vector{String}}
    R_w_DNS::Dict{String,Vector{String}}
    A_w::Dict{String,Vector{Tuple{String,String}}}
    S_w::Dict{String,Vector{String}}
    R_accept::Dict{Tuple{String,String},Vector{String}}

    c::Dict{Tuple{String,String},Float64}
    p::Dict{String,Float64}
    l::Dict{String,Float64}
    l_nr::Dict{String,Float64}
    l_scen::Dict{Tuple{String,String},Float64}
    remove::Dict{Tuple{String,String},Float64}
    accept::Dict{Tuple{String,String},Float64}
    v::Dict{String,Float64}

    W::Int64
    f::Float64

    insourcing::Dict{String,Bool}
    inswappable::Dict{String,Vector{String}}
    outswappable::Dict{String,Vector{String}}
    notswappable::Dict{String,Vector{String}}
    R_swap::Dict{String,Vector{String}}
end

"""
CONSTRUCTOR
Build DayBefore struct from data

### Keywords
* `data_dir` - location of CSV files
* `num_scenarios` - number of scenarios to sample from scenario dataset
* `num_itin` - number of itineraries to sample from itinerary set
* `W` - maximum number of workers
* `f` - fixed cost of hiring a adhoc
* `labor_min_planned` - minutely wage for planned drivers
* `labor_min_adhoc` - minutely wage for adhocs
* `delay_penalty` - cost of a delay
* `min_shift_length` - minimum shift length in minutes
### Returns
* DayBefore
"""
function DayBefore(
        data_dir::String,
        num_scenarios::Int64,
        num_itin::Int64,
        W::Int64,
        f::Float64,
        labor_min_planned::Float64,
        labor_min_adhoc::Float64,
        delay_penalty::Float64,
        min_shift_length::Float64=0.0,
        required_files::Vector{String}=REQUIRED_FILES
    )
    @assert issubset(joinpath.(data_dir, required_files), Glob.glob("*.csv", data_dir))
    data = Dict(replace(filename, ".csv" => "") => CSV.read(joinpath(data_dir, filename), DataFrame) for filename in required_files)

    #--- get first num_itin shortest itineraries
    itin_keep = data["itin_attributes"][.!occursin.(DNS_label, data["itin_attributes"][:,:effective_path_id]), :]
    itin_keep = leftjoin(data["itin"], itin_keep, on = (:path_id => :effective_path_id))
    itin_keep = sort(tuple.(itin_keep[:, :labor_min], itin_keep[:, :path_id]), by=first)[1:min(num_itin, nrow(data["itin_attributes"]))]

    # effective itinerary IDs
    effective_shift_attr = Dict{String,Tuple{Float64,Int64}}(data["itin_attributes"][:, :effective_path_id] .=> tuple.(data["itin_attributes"][:, :labor_min], round.(data["itin_attributes"][:, :delays], digits=0)))

    #--- sets
    R = string.(data["trips"].trip_id)
    A = [(row["origin_trip_id"], row["dest_trip_id"]) for row in eachrow(data["edges"]) if row["adhoc_edge"]]
    S = [path_id for (_, path_id) in itin_keep if path_id in keys(effective_shift_attr)]

    #--- index sets
    S_list = Dict{String,Vector{String}}(
        path_id => setdiff(split(path_id, "-"), ["src", "sink"]) for path_id in S
    )
    d_f = Dict{String,Vector{Tuple{String,String}}}(
        r => [(a, b) for (a, b) in A if a == r] for r in R
    )
    d_b = Dict{String,Vector{Tuple{String,String}}}(
        r => [(a, b) for (a, b) in A if b == r] for r in R
    )
    S_r = Dict{String,Vector{String}}(
        r => [path_id for (path_id, path) in S_list if r in path] for r in R
    )

    #--- sample scenarios
    scen_str = "scenarios.csv" in required_files ? "scenarios" : "nodns_scenarios"
    n = min(num_scenarios, length(unique(data[scen_str].scenario_id)))
    O = string.(unique(data[scen_str].scenario_id)[1:n])
    scen_temp = filter(row -> string(row[:scenario_id]) in O, data[scen_str])

    #--- build scenario index-sets
    R_w = Dict{String,Vector{String}}(w => String[] for w in O)
    R_w_DNS = Dict{String,Vector{String}}(w => String[] for w in O)
    for row in eachrow(scen_temp)
        scen_id = string(row[:scenario_id])
        trip_id = string(row[:trip_id])
        if trip_id == NOCANCEL_FLAG # TODO check this
            continue
        end
        if row[:dns]
            # DNS only
            push!(R_w_DNS[scen_id], trip_id)
        else
            # cancellations only
            push!(R_w[scen_id], trip_id)
        end
    end
    # all requests that are served (not cancelled but maybe DNS)
    R_w = Dict{String,Vector{String}}(w => setdiff(R, Rw) for (w, Rw) in R_w)
    A_w = Dict{String,Vector{Tuple{String,String}}}(
        w => [(a, b) for (a, b) in A if (a in R_w[w]) & (b in R_w[w])] for w in O
    )
    S_w = Dict{String,Vector{String}}(
        w => [i for i in S if length(intersect(S_list[i], R_w[w])) > 0] for w in O
    )


    #--- parameters
    p = Dict{String,Float64}(w => 1 / n for w in O)

    # straightforward shift costs
    l_nr = Dict{String,Float64}(
        path_id => labor_min_planned * effective_shift_attr[path_id][1] + delay_penalty * effective_shift_attr[path_id][2] for path_id in S
    )

    # expected shift costs
    l = Dict{String,Float64}(path_id => 0 for path_id in S)
    l_scen = Dict{Tuple{String,String},Float64}()
    S_effective = Dict{Tuple{String,String},String}()
    for path_id in S, w in O
        path = S_list[path_id]
        new_path = [path[i] in R_w_DNS[w] ? DNS_label : path[i] for i in 1:length(path) if path[i] in R_w[w]]
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
            S_effective[(w, path_id)] = join(vcat("src", new_path, "sink"), '-')
            if !(S_effective[w, path_id] in keys(effective_shift_attr))
                print(w, path_id)
            end
            labor, delays = effective_shift_attr[S_effective[w, path_id]]
            l_scen[path_id, w] = labor_min_planned * labor + delay_penalty * delays
            l[path_id] += p[w] * l_scen[path_id, w]
        else
            # if shift disappears post-cancellation, no cost incurred: p_w * 0
            l_scen[path_id, w] = 0
        end
    end
    # penalize the shifts that aren't long enough
    not_long_enough = [row[:effective_path_id] for row in eachrow(data["itin_attributes"]) if row[:end_min] - row[:start_min] <= min_shift_length]
    for path_id in intersect(not_long_enough, S)
        l[path_id] = SL_PENALTY
    end

    # marginal costs
    add = Dict{Tuple{String,String},Float64}(
        tuple.(data["add_trip"][:, "effective_path_id"], string.(data["add_trip"][:, "trip_id"])) .=> labor_min_planned * data["add_trip"][:, "delta_shiftlength"] .+ delay_penalty * data["add_trip"][:, "delta_delays"]
    )

    # marginal cost savings
    remove = Dict{Tuple{String,String},Float64}(
        tuple.(data["remove_trip"][:, "effective_path_id"], string.(data["remove_trip"][:, "trip_id"])) .=> labor_min_planned * data["remove_trip"][:, "delta_shiftlength"] .+ delay_penalty * data["remove_trip"][:, "delta_delays"]
    )

    # which requests can be accepted into an itinerary in each scenario?
    effective_paths = unique(values(S_effective))
    R_accept_effective = Dict{String,Vector{String}}(path_id => [] for path_id in unique(data["add_trip"][:, "effective_path_id"]))
    for (path_id, r) in keys(add)
        push!(R_accept_effective[path_id], r)
    end
    # for each itinerary without DNS
    R_accept = Dict{Tuple{String,String},Vector{String}}()
    for w in O, path_id in intersect(S_w[w], keys(R_accept_effective))
        if (w, path_id) in keys(S_effective)
            if S_effective[(w, path_id)] in keys(R_accept_effective)
                if isempty(intersect(S_list[path_id], R_w_DNS[w]))
                    R_accept[(w, path_id)] = intersect(R_accept_effective[S_effective[(w, path_id)]], R_w[w])
                end
            end
        end
    end

    # edge costs
    c = Dict{Tuple{String,String},Float64}(
        (row["origin_trip_id"], row["dest_trip_id"]) => labor_min_adhoc * (row["deadhead_time_min"] + row["idle_min"]) + delay_penalty * row["delay_proportion"] for row in eachrow(filter(row -> row["adhoc_edge"], data["edges"]))
    )

    # variable adhoc costs
    v = Dict{String,Float64}(
        string.(data["trips"][:, :trip_id]) .=> labor_min_adhoc * data["trips"][:, :travel_time_min]
    )

    #--- swapping model indices

    # all of the requests that can be swapped out in a given scenario
    R_swap = Dict{String,Vector{String}}(w => setdiff([r for ((v, _), R_a) in R_accept for r in R_a if v == w], R_w_DNS[w]) for w in O)

    # which itineraries can accept requests?
    inswappable = Dict{String,Vector{String}}(w => [] for w in O)
    for ((w, i), R_a) in R_accept
        if length(R_a) > 0
            push!(inswappable[w], i)
        end
    end

    # which itineraries have requests that are possible to insource? can't outswap DNS because that's not a removal choice
    outswappable = Dict{String,Vector{String}}(w => [i for i in S_w[w] if length(intersect(R_swap[w], S_list[i])) > 0] for w in O)

    # is this scenario amenable to insourcing?
    insourcing = Dict{String,Bool}(w => length(inswappable[w]) > 0 for w in O)

    # itineraries that are irrelevant to insourcing
    notswappable = Dict{String,Vector{String}}(w => setdiff(S_w[w], union(inswappable[w], outswappable[w])) for w in O)

    return DayBefore(R, A, S, O, S_list, S_effective, d_f, d_b, S_r, R_w, R_w_DNS, A_w, S_w, R_accept, c, p, l, l_nr, l_scen, remove, add, v, W, f, insourcing, inswappable, outswappable, notswappable, R_swap)
end
