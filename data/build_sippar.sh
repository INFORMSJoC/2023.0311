#!/bin/bash

# user specifications
CONFIG_LOC=$(python3 configure_parameters.py)
echo $CONFIG_LOC

# trips
Rscript simulate_rides.R "$CONFIG_LOC"

# edges, scenarios
Rscript build_VSNs_scenarios.R "$CONFIG_LOC"

# itineraries
python3 all_paths.py "$CONFIG_LOC"

# marginal costs
julia buildDayBefore.jl -c "$CONFIG_LOC"