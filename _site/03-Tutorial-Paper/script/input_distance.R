# dist_matrix <- matrix(0, nrow = length(states),
#                             ncol = length(states),
#                             dimnames = list(states, states))
# 
# for (i in 1:length(states)) {
#   print(i)
#   for (j in 1:length(states)) {
#     if (i != j) {
#       state1 <- states_sf[states_sf$NAME == states[i], ]
#       state2 <- states_sf[states_sf$NAME == states[j], ]
#       dist_matrix[i, j] <- min(st_distance(state1, state2))
#     }
#   }
# }
# save(dist_matrix, file="dist-USA.RData")
load(file="input/dist-USA.RData")

dist_matrix <- log(dist_matrix/100000+1)
range(dist_matrix)

y_log <- as.vector(sin(-dist_matrix/1.5))
dist <- as.vector(dist_matrix)
fit <- gam(y_log ~ s(dist) )
coefficients <- coef(fit)
spline_function <- function(dist_matrix, fit) {
  return(predict(fit, newdata = data.frame(dist=dist_matrix)))
}
# Visualization
data <- data.frame(x = as.vector(dist_matrix),
                   y = spline_function(as.vector(dist_matrix),
                                       fit))
coefficients_plot <- ggplot(data=data, aes(x, y)) +
  geom_line(data=data,
            linetype = 1, col=2,
            size=1)

coefficients_plot