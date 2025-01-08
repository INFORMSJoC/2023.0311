library(tidyverse) 
library(readr)
library(lubridate)
library(caret)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
config <- fromJSON(file=args[1])


# USER CONFIGURATION: number of rides in the SIPPAR instance
NUM_TRIPS <- config$num_trips
output_dir <- config$output_dir
SEED <- 4747

set.seed(SEED)
setwd(".")

# simulation tools
tt <- read_csv(file.path("data_sim", "travel_time.csv"))
prob_density <- read_csv(file.path("data_sim", "prob_cancel_density.csv"))
pickup_density <- read_csv(file.path("data_sim", "pickup_time_density.csv"))
trip_pattern <- read_csv(file.path("data_sim", "trip_pattern.csv")) %>%
  mutate(cdf = round(cdf, digits=6),
         lb = lag(cdf),
         lb = ifelse(is.na(lb), 0, lb - 1e-6))


### Sample cancellation probabilities 
# NOTE: A continuous empirical distribution can be leveraged with this input 
# using inverse transform style sampling.
cancel_prob <- sample(prob_density$prob_cancel, 
                      size = NUM_TRIPS, 
                      prob = prob_density$density, 
                      replace = TRUE) 

### Sample pickup times 
# NOTE: same inverse transform comment as above
req_pickup_time <- sample(pickup_density$pickup_t_min, 
                          size = NUM_TRIPS, 
                          prob = pickup_density$density, 
                          replace = TRUE) %>% 
  data.frame(pickup_t_min = .) %>%
  # convert to datetime
  mutate(pickup_t_hour = floor(pickup_t_min / 60),
         pickup_t_minute = floor(pickup_t_min - pickup_t_hour * 60),
         pickup_t = paste(pickup_t_hour, pickup_t_minute, sep=":"),
         pickup_t = as.POSIXct(pickup_t, format="%H:%M", tz="utc")) %>%
  select(-pickup_t_min, pickup_t_hour) %>%
  .$pickup_t

### Sample trip patterns
pattern <- c()
for (t in 1:NUM_TRIPS) {
  u <- runif(1) %>% round(digits=6)
  tp <- trip_pattern %>%
    filter(lb <= u,
           cdf > u) %>%
    .$trip_pattern_id
  pattern <- c(pattern, tp)
}

### Compile trip characteristics

rides <- data.frame(
  trip_pattern_id = pattern,
  cancel_prob = cancel_prob,
  req_pickup_dt = req_pickup_time) %>%
  left_join(trip_pattern %>%
              select(origin_id, 
                     dest_id,
                     trip_pattern_id), 
            by="trip_pattern_id") %>%
  select(-trip_pattern_id) %>%
  mutate(trip_id = 1:n()) %>%
  left_join(tt,
            by=c("origin_id", 
                 "dest_id")) %>% 
  mutate(req_dropoff_dt = req_pickup_dt + seconds(travel_time_min * 60))
rides %>% write_csv(file.path(output_dir, "trips.csv"))
