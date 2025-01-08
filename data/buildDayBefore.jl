using ArgParse
using CSV
using DataFrames
using Dates
using JSON
using JuMP
using Gurobi
using Printf


#---filepaths
TRIP_FN = "trips.csv"
ITIN_FN = "itin.csv"
ITIN_ATTR_FN = "itin_attributes.csv"
EDGE_FN = "edges.csv"
REMOVE_FN = "remove_trip.csv"
SCENARIO_FN = "scenarios.csv"
ADD_FN = "add_trip.csv"

SRC_ID = "src"
SINK_ID = "sink"


#---constants
ADD = true
GUROBI_ENV = Gurobi.Env()
WINDOW = 15.0
MIN_PER_HOUR = 60
BIG_M = 1.0e6
DNS_label = "dns"
LATE_LIMIT_MIN = 60


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config", "-c"
            help = "filepath to config file"
            required = true
    end

    return parse_args(s)
end


function get_fp(
        output_dir::String
    )
    itin_fp = joinpath(output_dir, ITIN_FN)
    itin_attr_fp = joinpath(output_dir, ITIN_ATTR_FN)
    trip_fp = joinpath(output_dir, TRIP_FN)
    edge_fp = joinpath(output_dir, EDGE_FN)
    scen_fp = joinpath(output_dir, SCENARIO_FN)
    remove_fp = joinpath(output_dir, REMOVE_FN)
    add_fp = joinpath(output_dir, ADD_FN)
    return itin_fp, itin_attr_fp, trip_fp, edge_fp, scen_fp, remove_fp, add_fp
end

"""
Read in CSV files. If `scenario_fp` or `itin_fp` are empty strings, nothing to read in.

### Keywords
* `itin_fp` - filepath to itinerary CSV
* `trip_fp` - filepath to trip CSV
* `edge_fp` - filepath to edge CSV
* `scenario_fp` - filepath to scenario CSV
### Returns
* DataFrame of itineraries
* Dict of itin ID to itinerary
* Dict of effective itin ID to itinerary (post-cancellation and post-DNS)
* DataFrame of trips
* DataFrame of edges
* DataFrame of scenarios
"""
function read_data(
        itin_fp::String,
        trip_fp::String,
        edge_fp::String,
        scenario_fp::String
    )
    # optional data
    if length(itin_fp) > 0
        itin = CSV.read(itin_fp, DataFrame)

        # format itineraries
        I_list = Dict{String,Vector{String}}(path_id => split(path_id, "-")[2:(end-1)] for path_id in itin[:, "path_id"])
    else
        itin, I_list = nothing, nothing
    end
    scenario = length(scenario_fp) > 0 ? CSV.read(scenario_fp, DataFrame) : nothing

    # rest of data
    trip = CSV.read(trip_fp, DataFrame)
    edge = CSV.read(edge_fp, DataFrame)

    #--- get effective itineraries post-cancellation
    I_list_effective = copy(I_list)
    scenario_ids = unique(scenario[:, :scenario_id])
    for scen_id in scenario_ids
        #--- compute effective path IDs across scenarios
        # an effective path ID is the set of all trips that haven't cancelled together with any intermediate DNS
        scen_temp = filter(row -> row[:scenario_id] == scen_id, scenario)
        cancelled = [string(row[:trip_id]) for row in eachrow(scen_temp) if !row[:dns]]
        dns = setdiff(string.(scen_temp[:, :trip_id]), cancelled)
        for (path_id, path) in I_list
            # empty itinerary check
            if isempty(setdiff(path, vcat(cancelled, dns)))
                continue
            end
            new_path = [path[i] in dns ? DNS_label : path[i] for i in 1:length(path) if !(path[i] in vcat(cancelled))]

            # get rid of head/tail dns
            while new_path[1] == DNS_label
                new_path = (new_path[1] == DNS_label) ? new_path[2:end] : new_path
            end
            while new_path[end] == DNS_label
                new_path = (new_path[end] == DNS_label) ? new_path[1:(end-1)] : new_path
            end

            # get rid of duplicate dns in middle
            new_new_path = [new_path[1]]
            for i in 2:length(new_path)
                if !((new_path[i] == DNS_label) & (new_new_path[end] == DNS_label))
                    push!(new_new_path, new_path[i])
                end
            end
            new_path = new_new_path

            # won't be empty due to initial empty itinerary check
            effective_path_id = join(vcat(SRC_ID, new_path, SINK_ID), '-')
            I_list_effective[effective_path_id] = new_path
        end
    end

    return itin, I_list, I_list_effective, trip, edge, scenario
end

"""
Get parameters needed for delay model
Deadhead times (minutes)
Trip travel times (minutes)
Requested pickups and dropoffs (normalized)
### Keywords
* `trip` - trip DataFrame
* `edge` - edge DataFrame
### Returns
* Dict of requested pickup times (trip => pickup minute)
* Dict of requested dropoff times (trip => dropoff minute)
* Dict of travel times (trip => time in mins)
* Dict of deadhead times (edge => time in mins)
"""
function times(
        trip::DataFrame,
        edge::DataFrame
    )
    #--- convert to datetime objects
    trip.req_pickup_dt = replace.(trip.req_pickup_dt, "Z" => "")
    trip.req_pickup_dt = DateTime.(trip.req_pickup_dt, "y-m-dTH:M:S")
    trip.req_dropoff_dt = replace.(trip.req_dropoff_dt, "Z" => "")
    trip.req_dropoff_dt = DateTime.(trip.req_dropoff_dt, "y-m-dTH:M:S")

    #--- format requested pickup times and dropoff times into dictionaries
    trip.trip_id = string.(trip.trip_id)
    tp = Dict(trip.trip_id .=> MIN_PER_HOUR * hour.(trip.req_pickup_dt) .+ minute.(trip.req_pickup_dt))
    td = Dict(trip.trip_id .=> MIN_PER_HOUR * hour.(trip.req_dropoff_dt) .+ minute.(trip.req_dropoff_dt))
    tt = Dict(trip.trip_id .=> trip.travel_time_min)
    tt_d = Dict(tuple.(edge.origin_trip_id, edge.dest_trip_id) .=> edge.deadhead_time_min)

    return tp, td, tt, tt_d
end


########################################
################ ITINERARY PARAMS
########################################

"""
Itinerary cost model
### Keywords
* `I` - ordered list of trips
* `tt` - travel time dict
* `tt_d` - deadhead time dict
* `tp` - requested pickup time dict
* `P` - binary delay penalty
* `C` - minute-wise cost of labor
* `window` - time window radius
* `M` - a very large number
### Returns
* delay model
"""
function itin_cost_model(
        I::Vector{String},
        tt::Dict,
        tt_d::Dict,
        tp::Dict,
        P::Float64=PENALTY,
        C::Float64=LABOR_MIN,
        window::Float64=WINDOW,
        M::Float64=BIG_M
    )
    #--- get DNS
    R_dns = []
    for (i, r) in enumerate(I)
        if r == DNS_label
            push!(R_dns, I[i-1])
        end
    end
    R_dns_trip = setdiff(R_dns, [DNS_label])

    I = setdiff(I, [DNS_label])
    n = length(I)
    R_dns = [i for i in 1:n if I[i] in R_dns_trip]

    #--- build model
    m = JuMP.Model(JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0))

    JuMP.@variable(m, pickup_time[I] >= 0)
    JuMP.@variable(m, z[I], Bin)

    JuMP.@expressions m begin
        # start and end times include depot travel
        start_time, pickup_time[I[1]] - tt_d[("src", I[1])]
        end_time, pickup_time[I[n]] + tt[I[n]] + tt_d[(I[n], "sink")]
        delays, sum(z)
    end

    # don't count labor during DNS periods
    if isempty(R_dns)
        JuMP.@expression(m, labor, end_time - start_time)
    else
        JuMP.@expression(m, labor, end_time - start_time - sum(pickup_time[I[r+1]] - (pickup_time[I[r]] + tt[I[r]]) for r in R_dns))
    end

    JuMP.@objective(m, Min,
        labor * C + delays * P
    )

    JuMP.@constraints m begin
        # trip 1 on time
        pickup_time[I[1]] <= tp[I[1]] + window
        # all trips picked up on time or late
        [r in I], pickup_time[r] >= tp[r] - window
        # no trips picked up too late
        [r in I], pickup_time[r] - (tp[r] + window) <= LATE_LIMIT_MIN
        # pickup after arrival
        [j in 2:n], pickup_time[I[j-1]] + tt[I[j-1]] + tt_d[(I[j-1], I[j])] <= pickup_time[I[j]]
        # whether to count penalty
        [r in I], pickup_time[r] - (tp[r] + window) <= M * z[r]
    end

    return m
end

"""
For an itinerary, compute:
* shift length
* delays
* start time (min)
* end time (min)

### Keyword Args
* `I` - ordered list of trips
* `tt` - travel time dict
* `tt_d` - deadhead time dict
* `tp` - requested pickup time dict
* `P` - binary delay penalty
* `C` - minute-wise cost of labor
* `window` - time window radius
* `M` - a very large number
### Returns
* optimal number of planned delays
* itinerary start time in mins
* itinerary end time in mins
* labor in minutes
"""
function itin_cost(
        I::Vector{String},
        tt::Dict,
        tt_d::Dict,
        tp::Dict,
        P::Float64=PENALTY,
        C::Float64=LABOR_MIN,
        window::Float64=WINDOW,
        M::Float64=BIG_M
    )
    if length(I) == 0
        delays = 0
        shift_start = missing
        shift_end = missing
        labor = missing
    elseif length(I) == 1
        delays = 0
        shift_start = tp[I[1]] - tt_d[("src", I[1])]
        shift_end = tp[I[1]] + tt[I[1]] + tt_d[(I[1], "sink")]
        labor = shift_end - shift_start
    else
        m = itin_cost_model(I, tt, tt_d, tp, P, C, window, M)
        JuMP.optimize!(m)
        if JuMP.termination_status(m) == MOI.OPTIMAL
            delays = JuMP.value(m[:delays])
            shift_start = JuMP.value(m[:start_time])
            shift_end = JuMP.value(m[:end_time])
            labor = JuMP.value(m[:labor])
        else
            delays, shift_start, shift_end, labor = nothing, nothing, nothing, nothing
        end
    end
    return delays, shift_start, shift_end, labor
end



"""
For all post-cancellation and -DNS itineraries (dependent on scenarios), compute:
- number of delays
- start time
- end time
- labor in minutes

### Returns 
- Dict: path ID => (minutes of labor, number of delays)
- Dict: path ID => (shift start, shift end)
"""
function all_shift_costs(
        output_fp::String,
        I_list_effective::Dict{String,Vector{String}},
        tt::Dict,
        tt_d::Dict,
        tp::Dict,
        P::Float64=PENALTY,
        C::Float64=LABOR_MIN,
        window::Float64=WINDOW,
        M::Float64=BIG_M
    )
    num_itineraries = length(I_list_effective)
    count = 0
    shift_attr = Dict{String,Tuple{Float64,Int64}}()
    shift_timestamps = Dict{String,Tuple{Float64,Float64}}()
    for (path_id, I) in I_list_effective
        count += 1
        if count % 1000 == 0
            @printf("Calculated costs for %i out of %i itineraries\n", count, num_itineraries)
        end
        thecosts = itin_cost(I, tt, tt_d, tp, P, C, window, M)
        # don't include planned itineraries that require cancellations
        if any(isnothing.(thecosts))
            continue 
        end
        delays, shift_start, shift_end, labor_min = thecosts
        shift_attr[path_id] = (labor_min, delays)
        shift_timestamps[path_id] = (shift_start, shift_end)
    end
    @printf("Calculated all %i itineraries\n", num_itineraries)
    return shift_attr, shift_timestamps
end


"""
Generate remove_trip.csv or add_trip.csv

### Keywords
* `output_fp` -  filepath to write outputs
* `add` - whether we are computing marginal cost savings to add trips (as opposed to marginal costs to remove trips)
* `I_list_effective` - Dictionary of post-cancellation itineraries (path_id => ordered list of trip IDs served)
* `shift_attr` - Dict of shift lengths and delay quantities (path_id => (shift length in minutes, num delays))
* `trips` - list of trip IDs
* `tt` - travel times
* `tt_d` - deadhead times
* `tp` - requested pickup times
### Writes
* DataFrame: add_trip.csv if `add` else remove_trip.csv
"""
function all_marginal_itin(
        output_fp::String,
        add::Bool,
        I_list_effective::Dict{String,Vector{String}},
        shift_attr::Dict{String,Tuple{Float64,Int64}},
        shift_timestamps::Dict{String,Tuple{Float64,Float64}},
        trips::Vector{String},
        tt::Dict,
        tt_d::Dict,
        tp::Dict
    )
    f = open(output_fp, "w")
    write(f, "effective_path_id,trip_id,delta_shiftlength,delta_delays\n")

    # test every itin/trip combination
    for (path_id, I) in I_list_effective
        # model assumption: don't add any trips to DNS itineraries
        if (DNS_label in I) & add
            continue
        end

        # edge case: remove a trip from a single-trip itinerary
        old_shiftlength, old_delays = shift_attr[path_id]
        if (length(I) == 1) & !add
            output_str = join([path_id, I[1], old_shiftlength, old_delays], ',') * "\n"
            write(f, output_str)
            continue
        end

        if add 
            # consider adding all trips not in the itinerary 
            thetrips = setdiff(trips, I) 
        else 
            # consider removing all trips from the itinerary
            thetrips = setdiff(I, [DNS_label])
        end

        for trip_id in thetrips
            if add
                # all simple insertions
                paths_to_test = [join(vcat(["src"], I[1:i], [trip_id], I[(i+1):end], ["sink"]), "-") for i in 0:length(I)]
            else
                candidate_path = [I[i] for i in 1:length(I) if I[i] != trip_id]

                # remove start/end/duplicate DNS labels
                candidate_path = candidate_path[findfirst(candidate_path .!= DNS_label):findlast(candidate_path .!= DNS_label)]
                thepath = [candidate_path[1]]
                for j in 2:length(candidate_path)
                    t = candidate_path[j]
                    if !((t == DNS_label) & (thepath[end] == DNS_label))
                        # if not a duplicate DNS label
                        push!(thepath, t)
                    end
                end
                paths_to_test = [join(vcat(["src"], thepath, ["sink"]), "-")]
            end

            # calculate the leftover itinerary for each removed trip if unseen
            if !add & isempty(intersect(paths_to_test, keys(shift_attr))) 
                remove_path = paths_to_test[1]
                I = setdiff(String[s for s in split(remove_path, "-")], ["src", "sink"])
                delays, shift_start, shift_end, labor = itin_cost(I, tt, tt_d, tp)
                shift_attr[remove_path] = (labor, delays)
                shift_timestamps[remove_path] = (shift_start, shift_end)
            end
            paths_to_test = intersect(paths_to_test, keys(shift_attr))

            # if there are any candidates, find the best and compute marginal cost
            if length(paths_to_test) > 0
                thecosts = [LABOR_MIN * shift_attr[path][1] + PENALTY * shift_attr[path][2] for path in paths_to_test]
                new_path_id = paths_to_test[thecosts .== minimum(thecosts)][1]

                # params of interest
                new_shiftlength, new_delays = shift_attr[new_path_id]
                if add
                    delta_shiftlength = new_shiftlength - old_shiftlength
                    delta_delays = new_delays - old_delays
                else
                    delta_shiftlength = old_shiftlength - new_shiftlength
                    delta_delays = old_delays - new_delays
                end

                # record the outputs
                output_str = join([path_id, trip_id, delta_shiftlength, delta_delays], ',') * "\n"
                write(f, output_str)
            end
        end
    end
    close(f)
    return shift_attr, shift_timestamps
end

function write_itinerary_attributes(
        shift_attr, 
        shift_timestamps,
        output_fp:: String
    )
    f = open(output_fp, "w")
    write(f, "effective_path_id,delays,start_min,end_min,labor_min\n")
    for path_id in keys(shift_attr)
        labor_min, delays = shift_attr[path_id]
        shift_start, shift_end= shift_timestamps[path_id]
        output_str = join([path_id, delays, shift_start, shift_end, labor_min], ',') * "\n"
        write(f, output_str)
    end
    close(f)
end

########################################
################ MAIN
########################################


function main()
    itin_fp, itin_attr_fp, trip_fp, edge_fp, scen_fp, remove_fp, add_fp = get_fp(OUTPUT_DIR)

    #--- set up data
    itin, I_list, I_list_effective, trip, edge, scenario = read_data(itin_fp, trip_fp, edge_fp, scen_fp)
    tp, td, tt, tt_d = times(trip, edge)
    trips = string.(trip.trip_id)

    print("Calculating itinerary costs.\n")

    # record original (first-stage) itinerary attributes 
    shift_attr, shift_timestamps = all_shift_costs(itin_attr_fp, I_list_effective, tt, tt_d, tp)
    # edit any infeasible itineraries out of I_list_effective
    I_list_effective = Dict{String,Vector{String}}(
        path_id => path 
        for (path_id, path) in I_list_effective 
            if path_id in keys(shift_attr)
    )

    println("Calculating marginal costs to remove trips. ")
    shift_attr, shift_timestamps = all_marginal_itin(remove_fp, !ADD, I_list_effective, shift_attr, shift_timestamps, trips, tt, tt_d, tp)
    
    println("Calculating marginal cost savings to insert trips.")
    shift_attr, shift_timestamps = all_marginal_itin(add_fp, ADD, I_list_effective, shift_attr, shift_timestamps, trips, tt, tt_d, tp)
    
    write_itinerary_attributes(shift_attr, shift_timestamps, itin_attr_fp)

    println("Done with itinerary additions.")
    println("Done with all.")
end

# Get the config filepath from command line arguments 
parsed_args = parse_commandline()
config_filepath = parsed_args["config"]
config_file = open(config_filepath, "r") 
config_data = JSON.parse(read(config_file, String)) 
close(config_file)

# build SIPPAR parameters
LABOR_MIN = config_data["hourly_wage"] / MIN_PER_HOUR
PENALTY = convert(Float64, config_data["penalty"])
OUTPUT_DIR = config_data["output_dir"]
main()