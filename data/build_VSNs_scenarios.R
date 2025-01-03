"
Build VSNs and uncertain scenarios.

### INPUTS 
* trips.csv
* travel times 
### OUTPUTS 
* edges.csv
* scenarios.csv
"
library(tidyverse)
library(readr)
library(lubridate)
library(igraph)
library(caret)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(file=args[1])

# input parameters
PLANNED_HOURLY_WAGE <- config$hourly_wage # worker hourly wage
PENALTY <- config$penalty                 # penalty for operator to delay a trip (same as P)
NUM_SCENARIOS <- config$num_scenarios     # number of in-sample scenarios

# filepaths
setwd(".")
output_dir <- config$output_dir
print(output_dir)
TRIP_FN <- file.path(output_dir, "trips.csv")
EDGE_FN <- file.path(output_dir, "edges.csv")
SCEN_FN <- file.path(output_dir, "scenarios.csv")


# constants
MIN_PER_HOUR <- 60
SRC_ID <- "src"
SINK_ID <- "sink"
DEPOT_ID <- 0
DNS_PROB <- 0.07
MAX_DNS <- 0.20
HALF_WINDOW <- 15         # time HALF_WINDOW in minutes (same as Delta/2)
MAX_CANCELLATIONS <- 4    # maximum trip cancellations in a scenario


#--- data
trip <- read_csv(file.path("case_study", "trips.csv"))
tt <- read_csv(file.path("data_sim", "travel_time.csv"))

set.seed(4747)



########################## BUILD EDGES


# look at every possible OD combo for drivers
edge <- expand.grid(origin_trip_id = trip$trip_id,
                    dest_trip_id = trip$trip_id) %>%
  filter(origin_trip_id != dest_trip_id) %>%
  # get pickup and drop-off times
  left_join(select(trip,
                   origin_trip_id=trip_id,
                   origin_node_id=dest_id,
                   origin_pickup_dt = req_pickup_dt,
                   origin_dropoff_dt = req_dropoff_dt),
            by="origin_trip_id") %>%
  left_join(select(trip,
                   dest_trip_id=trip_id,
                   dest_node_id=origin_id,
                   dest_pickup_dt = req_pickup_dt,
                   dest_dropoff_dt = req_dropoff_dt),
            by="dest_trip_id") %>%
  # get deadhead times
  left_join(tt %>%
              rename(origin_node_id=origin_id,
                     dest_node_id=dest_id),
            by=c("origin_node_id",
                 "dest_node_id")) %>%
  filter() %>%
  mutate(deadhead_time_min = ifelse(origin_node_id==dest_node_id,
                                0,
                                travel_time_min),
         # t_i^d - w + T_ij <= t_j^p + w ==> sequence is physically possible w/o delay
         planned_driver_edge = origin_dropoff_dt - HALF_WINDOW + deadhead_time_min < 
           dest_pickup_dt + HALF_WINDOW,
         # these conditions guarantee the adhoc VSN is acyclic
         adhoc_edge = (origin_pickup_dt < dest_pickup_dt) | 
           ((origin_pickup_dt == dest_pickup_dt) & 
              (origin_dropoff_dt < dest_dropoff_dt))) %>%
  # cost components
  mutate(dropoff_min = MIN_PER_HOUR * hour(origin_dropoff_dt) + 
           minute(origin_dropoff_dt),
         pickup_min = MIN_PER_HOUR * hour(dest_pickup_dt) + 
           minute(dest_pickup_dt),
         idle_min = pmax(0, pmin(pickup_min - dropoff_min - deadhead_time_min, 
                             2*HALF_WINDOW)) / 2 + 
           pmax(0, pickup_min - dropoff_min - 2*HALF_WINDOW - deadhead_time_min),
         # eq (27) Appendix A.2.1
         delay_proportion = pmax(0, pmin(1, (deadhead_time_min - pickup_min + dropoff_min) / 
                               (2 * HALF_WINDOW))),
         # approximate cost for k-shortest paths
         edge_cost = PLANNED_HOURLY_WAGE / MIN_PER_HOUR * (deadhead_time_min + idle_min) +
           PENALTY * delay_proportion) %>%
  select(-c(dropoff_min, 
            pickup_min,
            origin_dropoff_dt, 
            dest_pickup_dt,
            travel_time_min,
            origin_node_id,
            dest_node_id,
            dest_dropoff_dt,
            origin_pickup_dt)) %>%
  mutate(depot = F)

#--- build source/sink edges to/from depot
edge <- data.frame(origin_trip_id = trip$trip_id,
                      dest_trip_id = SINK_ID,
                      origin_node_id = trip$dest_id,
                      dest_node_id = DEPOT_ID) %>%
  rbind(data.frame(origin_trip_id = SRC_ID,
                   dest_trip_id = trip$trip_id,
                   origin_node_id = DEPOT_ID, 
                   dest_node_id = trip$origin_id)) %>%
  left_join(tt %>%
              rename(origin_node_id=origin_id,
                     dest_node_id=dest_id),
            by=c("origin_node_id",
                 "dest_node_id")) %>%
  rename(deadhead_time_min = travel_time_min) %>%
  select(-c(origin_node_id,
            dest_node_id)) %>%
  mutate(depot = T,
         adhoc_edge = F,
         planned_driver_edge = T,
         idle_min = 0,
         delay_proportion = 0,
         edge_cost = PLANNED_HOURLY_WAGE / MIN_PER_HOUR * deadhead_time_min) %>%
  rbind(edge)

edge %>% write_csv(EDGE_FN)

########################## CANCELLATION SCENARIOS

scenarios <- data.frame(scenario_id = c(),
                        trip_id = c(),
                        dns = c())
while (length(unique(scenarios$scenario_id)) < NUM_SCENARIOS) {
  # cancel some trips 
  try_again = T 
  cancellations <- c()
  while (try_again) {
    cancellations <- trip[runif(nrow(trip)) <= trip$cancel_prob,]$trip_id %>% sort()
    if ((length(cancellations) <= MAX_CANCELLATIONS) & (length(cancellations) > 0)) {
      try_again <- F
    }
  }

  # DNS some trips
  dns_candidates <- setdiff(trip$trip_id, cancellations)
  try_again = T
  dns_trips <- c()
  while (try_again) {
    dns_trips <- dns_candidates[runif(length(dns_candidates)) <= DNS_PROB] %>% sort()
    if (length(dns_trips) / nrow(trip) <= MAX_DNS) {
      try_again <- F
    }
  }

  # identify scenario
  id <- paste("cancel", paste(cancellations, collapse=",")) %>%
    paste(., paste("dns", paste(dns_trips, collapse=",")))
  if (!(id %in% scenarios$scenario_id)) {
    scenarios <- data.frame(scenario_id = id,
               trip_id = cancellations,
               dns = F) %>% 
      rbind(scenarios)
    if (length(dns_trips) > 0) {
      scenarios <- data.frame(scenario_id = id,
                 trip_id = dns_trips,
                 dns = T) %>%
        rbind(scenarios)
    }
    # If you want: check for the scenario with no cancellations or DNS here
  }
}

scenarios %>% write_csv(SCEN_FN)