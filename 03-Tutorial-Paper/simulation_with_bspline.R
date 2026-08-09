source("script/packages.R")
source("script/states.R")
source("script/input_distance.R")
library(mgcv)

# m + 1 is the order of the spline
# e.g. m = 2 for cubic splines
# m + 1 is the degree of the polynomial
# m + 2 is the order of the polynomial
# k: number of parameters / control points
# k + m + 2 knots
# x1, x2, ..., x_{k+m+2}

# Evaluate the ith B-spline basis function of order m at the values in x
# given knot locations in k - default is the cubic B-spline
bspline <- function(x, k, i, m = 2) {
  if (m == -1) {
    # Base case of the recursion:
    # when m is -1, we are at the lowest order basis function
    res <- as.numeric(x < k[i + 1] & x >= k[i])
    # Returns 1 if x is within the interval [k[i], k[i+1])
  } else {
    # Recursive case: construct the basis function
    # from lower order basis functions
    
    # Calculate the first term's coefficient
    z0 <- (x - k[i]) / (k[i + m + 1] - k[i])
    # Calculate the second term's coefficient
    z1 <- (k[i + m + 2] - x) / (k[i + m + 2] - k[i + 1])
    
    # Recursive calls to the lower order basis functions
    res <- z0 * bspline(x, k, i, m - 1) + z1 * bspline(x, k, i + 1, m - 1)
  }
  # Return the evaluated B-spline basis function
  return(res)
}

bspline_basis <- function(x, df, m = 2, xlim = c(0, 1)) {
  mg <- m + 1
  nk <- df - mg + 1
  if (nk <= 0) stop("basis dimension too small for b-spline order")
  
  xl <- xlim[1]
  xu <- xlim[2]
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl) / (nk - 1)
  knots <- seq(xl - dx * mg, xu + dx * mg, length = nk + 2 * mg)
  
  X <- sapply(1:df, function(i) bspline(x, knots, i, m))
  
  return(list(X = X, knots = knots))
}

k_dim <- 10
m <- 2
set.seed(1)
beta_true <- runif(k_dim, -1, 1)
dist_range <- range(dist_matrix, na.rm = TRUE)

cols <- c(
  "#A9A9A9",  # black
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#66A61E",  # green
  "#A6761D"   # brown
)

xg <- seq(dist_range[1], dist_range[2], length.out = 100)
b  <- bspline_basis(x = xg, df = k_dim, m = m, xlim = dist_range)
B  <- b$X
knots <- b$knots

sm <- smoothCon(s(xg, bs = "bs"), data.frame(xg))[[1]]
sm$knots

# check that this coincides with the Bspline construction in mgcv
knots == sm$knots
all(round(B, 8) == round(sm$X, 8))

XB <- sweep(B, 2, beta_true, "*")

common_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

plot(xg, rowSums(XB), type = "l", col = 1, lwd = 2)
fx <- function(x, xlim, beta, df = 10, m = 2) {
  B <- bspline_basis(x = x, df = k_dim, m = m, xlim = xlim)$X
  XB <- sweep(B, 2, beta, "*")
  return(rowSums(XB))
}
lines(xg, fx(xg, xlim = c(min(dist_matrix), max(dist_matrix)), beta = beta_true),
      col = 2, lty = 3, lwd = 2)

# ---- Multi-panel basis functions plot ----
pdf("pictures/modelling/bspline-basis-functions.pdf", width = 11, height = 8)
op <- par(
  mfrow = c(2, 5),
  mar = c(4, 4, 3, 1),
  mgp = c(2.2, 0.7, 0),
  cex.main = 2.0,   # <- larger panel titles
  cex.lab  = 1.6,
  cex.axis = 1.4
)
for (i in 1:k_dim) {
  plot(xg, B[, i],
       type = "l",
       lwd = 2,
       col = cols[i],
       ylim = c(0, max(B)),
       xlab = "x",
       ylab = bquote(b.(i)),
       main = bquote(b[.(i)]))
}
dev.off()
par(op)

op <- par(
  mfrow = c(1, 1),
  mar = c(4, 4, 3, 1),
  mgp = c(2.2, 0.7, 0),
  cex.main = 1.4,
  cex.lab = 1.4,
  cex.axis = 1.2
)
matplot(xg, B,
        type = "l",
        lwd = 2,
        lty = 1,
        col = cols,
        xlab = "x",
        ylab = expression(bj),
        main = "Cubic B-spline basis functions",
        ylim = c(0, 1.3))
legend("topright",
       legend = paste0("b", 1:k_dim),
       col = cols,
       lwd = 3,
       bty = "n",
       cex = 1.3,     # <- larger legend text
       ncol = 2)
par(op)

# ---- Linear basis expansion plot ----
pdf("pictures/modelling/linear-basis-expansion.pdf", width = 11, height = 8)
op <- par(
  mar = c(5, 5, 4, 1),
  mgp = c(3, 0.9, 0),
  cex.main = 2.2,   # <- larger main title
  cex.lab  = 1.7,
  cex.axis = 1.4
)
matplot(xg, XB,
        type = "l",
        lwd = 2,
        lty = 2,
        col = cols,
        xlab = "Distance",
        ylab = expression(theta[j] * bj),
        main = "Linear Basis Expansion",
        ylim = c(-0.6, 1.3))

lines(xg, rowSums(XB), lwd = 3, col = "black")

abline(h = 0, col = "grey85", lty = 3)

legend("topright",
       legend = c(sapply(1:k_dim, function(i) as.expression(bquote(theta[.(i)] * b[.(i)]))),
                  expression(f(x))),
       col = c(cols, "black"),
       lwd = c(rep(3, k_dim), 4),   # slightly thicker keys to match bigger text
       lty = c(rep(2, k_dim), 1),
       bty = "n",
       cex = 1.5,       # <- larger legend text
       ncol = 2)
dev.off()
par(op)

event_simulator <- function(states, p, n, seed, beta_true, xlim) {
  
  set.seed(seed)
  print(seed)

  current.tm <- 0
  
  tms <- numeric(n)
  dat <- matrix(0, nrow = n, ncol = 5)       
  nondat <- matrix(0, nrow = n, ncol = 4)
  
  last_contacted <- matrix(NA, nrow = n + 1, ncol = p)
  last_contacted[1, ] <- sample(states, p)
  colnames(last_contacted) <- states
  
  dist.m <- t(sapply(states, function(st)
    sapply(states, function(pot.st)
      dist_matrix[pot.st, last_contacted[1, st]])))
  distance <- matrix(fx(x = as.vector(dist.m), 
                        xlim = xlim, beta = beta_true), nrow = p, ncol = p)
  
  p2 <- p^2
  all_ids <- 1:p2
  
  for (i in 1:n) {
    if (i %% 100 == 0) print(i)
    
    hazard <- exp(distance) / p2
    haz_sum <- sum(hazard)
    dt <- rexp(1, haz_sum)
    current.tm <- current.tm + dt
    tms[i] <- current.tm
    
    event.id <- sample(all_ids, 1, prob = hazard / haz_sum)
    s <- (event.id - 1) %% p + 1
    r <- (event.id - 1) %/% p + 1
    dat[i, ] <- c(current.tm, s, r, event.id, dist.m[s, r])
    
    non.event.id <- sample(all_ids[-event.id], 1)
    non.s <- (non.event.id - 1) %% p + 1
    non.r <- (non.event.id - 1) %/% p + 1
    nondat[i, ] <- c(non.s, non.r, non.event.id, dist.m[non.s, non.r])
    
    last_contacted[i + 1, ] <- last_contacted[i, ]
    last_contacted[i + 1, states[s]] <- states[r]
    
    dist.m <- t(sapply(states, function(st)
      sapply(states, function(pot.st)
        dist_matrix[pot.st, last_contacted[i + 1, st]])))
    distance <- matrix(fx(x = as.vector(dist.m), 
                          xlim = xlim, beta = beta_true), nrow = p, ncol = p)
  }
  
  remdat <- cbind(1:n, dat, nondat, dat[, 5] - nondat[, 4], 1)
  colnames(remdat) <- c("index", "tms", "s", "r", "event.id",
                        "distance", "non.s", "non.r", "non.event.id",
                        "non.distance", "delta.dist", "one")
  remdat <- as.data.frame(remdat)
  return(list(remdat, last_contacted))
}

n_sim <- 1
set.seed(1234)
remdats <- sapply(1:n_sim, function(x)
  event_simulator(states = states, p = length(states),
                  n = n, seed = x, beta_true = beta_true, 
                  xlim = c(min(dist_range), max(dist_range))),
  simplify = FALSE)

remdats_data <- lapply(remdats, function(x) x[[1]])

gam.fits <- lapply(remdats_data, function(x) {
  dist_MAT <- cbind(x$distance, x$non.distance)
  ONE <- cbind(x$one, -x$one)
  fit <- gam(one ~ s(dist_MAT, by = ONE, bs = "ps", k = k_dim) - 1,
             family = "binomial", data = x)
  return(fit)
})

b.dist_shift.objects <- lapply(gam.fits, getViz)
xg <- seq(min(dist_range), 
          max(dist_range), length = 200)
data <- data.frame(x = xg,
                   y = fx(xg, xlim = dist_range, beta_true)-0.4)

coefficients_plot <- ggplot(data=data, aes(x, y)) +
  geom_line(data=data,
            linetype = 1, col="lightcoral",
            size=1)

for (iter in 1:n_sim){
  b.iter <- b.dist_shift.objects[[iter]]
  o.iter <- plot(sm(b.iter, 1))
  coefficients_plot <- coefficients_plot + 
    geom_line(data=o.iter$data$fit, 
              linetype = 3, col="lightgray", 
              size=0.5)
}

coefficients_plot <- coefficients_plot +
  labs(
    title = "Estimated vs True Non-linear Effect",
    x = "Covariate",
    y = "Log-hazard contribution"
  ) +
  common_theme +
  theme(legend.position = "none")
coefficients_plot
