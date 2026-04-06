# size of simulation
states <- states(cb = TRUE)
# Convert to simple feature (sf) object
states_sf <- st_as_sf(states)
states <- unique(states$NAME)
p <- length(states)
n <- 1000