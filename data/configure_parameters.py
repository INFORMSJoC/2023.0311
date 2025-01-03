import json 
import os

# USER CONFIGURES SIMULATION HERE
config = {
    # number of requested trips
    "num_trips" : 30,
    # output directory for case study
    "output_dir" : "case_study",
    # employee driver hourly wage in $
    "hourly_wage" : 20,
    # $ penalty for operator to pickup one passenger with delay
    "penalty" : 100,
    # number of uncertainty scenarios
    "num_scenarios" : 100,
    # max number of generated employee driver itineraries
    "num_itin" : 10000
}

if not os.path.exists(config["output_dir"]):
    os.makedirs(config["output_dir"])
config_fn = os.path.join(config["output_dir"], "config.json")
with open(config_fn, "w") as f:
    json.dump(config, f, indent=4) 

# for bash 
print(config_fn)