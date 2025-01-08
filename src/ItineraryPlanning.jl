module ItineraryPlanning
    using ArgParse
    using CSV
    using DataFrames
    using Glob
    using Gurobi
    using JuMP
    using Printf
    using StatsBase

    # data structure 
    export DayBefore

    # local pareto optimality cuts
    export ParetoArgugments
    export core_point, act_core_point, act_core_point_dynamic

    # regular and activated benders decomposition
    export MasterProblem, FullSubproblemDual, ActivatedSubproblemDual, FullSubproblemParetoDual, ActivatedSubproblemParetoDual, SecondStageScenario, SecondStageSolutions
    export solve_MP, solve_SP, benders

    # direct and OOS implementation
    export direct_model, oos_scenario
    export append_callback, direct_model_callback, direct_model_solve, solve_oos_model

    # utilities
    export get_data_path, record_xval, record_sol, get_roster, service_metrics, all_metrics, adhoc_hours, recourse_itin

    include("DayBefore.jl")
    include("benders.jl")
    include("direct.jl")
    include("utils.jl")

    const GRB = Gurobi.Env()
    const MAX_ITER = 1000          # benders
    const MAX_SPEEDY_ITER = 15
    const LARGE_NUMBER = 1e6
    const SP_TIME_LIMIT = 60 * 10       # time limits all in seconds
    const MP_TIME_LIMIT = 60 * 30
    const MP_TIME_LIMIT_MAX = 60 * 60
    const TIME_LIMIT = 60 * 60
    const SPEEDY_MIP_GAP = 5e-3
    const MIP_GAP = 1e-4
    const DEBUG_LABEL = "dual_infeas_"
    const RAND_CUT = 0.5
    const MIN_WORKERS = true
    const MIN_COST = true
    const GAP = 1e-4
    const DNS_label = "dns"
    const NULLPARETO = ParetoArguments(false, false, false, 0, 0)
end