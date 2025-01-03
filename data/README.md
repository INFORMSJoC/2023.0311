# Synthetic data and case study creation

The data in "Activated Benders Decomposition for Day-ahead Paratransit Planning" is proprietary. We have created a module to simulate synthetic trips and planning attributes so that others can test our methodology. 

## Simulate a SIPPAR instance

The simulator inputs are in ./data_sim. Trip origins and destinations are treated as nodes in an undirected graph with deterministic travel times. 

**Data schema for trip simulator**
* data_sim/trip_pattern.csv 
  - origin_id : origin node of trip 
  - dest_id : destination node of trip 
  - trip_pattern_id : ID for OD pair 
  - cdf : point in empirical CDF
* data_sim/travel_time.csv
  - travel_time_min : travel time in minutes between `origin_id` and `dest_id`
  - origin_id : origin node of travel 
  - dest_id : destination node of travel
* pickup_time_density.csv
  - pickup_t_min : requested pickup time - units are minutes elapsed since midnight
  - density : probability of that pickup time 
* data_sim/prob_cancel_density.csv
  - prob_cancel : probability of cancellation
  - density : probability of that cancellation probability

Empirical densities for requested pickup times and cancellation probabilities were created using kernel density estimation (`density` function in *R* with bw = "nrd0" and n=1000). Cancellation probability density was post-processed to exclude negative values. Travel times were calculated using OSRM R-package (`osrmTable` function).

**Run the simulation**

To simulate SIPPAR model inputs, complete the following steps. 

1. Configure parameters in `configure_parameters.py`:
    * num_trips : number of requested trips 
    * output_dir : output directory for case study
    * hourly_wage : employee driver hourly wage
    * penalty : penalty for operator to delay a trip 
    * num_scenarios : number of uncertainty scenarios (both DNS and cancellation events)
    * num_itin : max number of itineraries generated

2. Run data/build_sippar.sh, which executes the following:
    * Run `configure_parameters.py` (config.json)
    * Run `simulate_rides.R` (trips)
    * Run `build_VSNs_scenarios.R` (edges, scenarios)
    * Run `all_paths.py` (itineraries)
    * Run `buildDayBefore.jl` (itinerary attributes, marginal costs)


The authors set a seed of "4747" when simulating the provided case study (in data/case_study). Be sure to reset the seed in `simulate_rides.R` and `build_VSNs_scenarios.R` to generate a different case study.