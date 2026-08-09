set.seed(1234)
dist_matrix <- matrix(rnorm(p^2, mean = 0, sd = 1), nrow = p, ncol = p)

bspline <- function(x, k, i, m = 2) {
  
  if (m == -1) {
    res <- as.numeric(x < k[i + 1] & x >= k[i])  
  } else {
    z0 <- (x - k[i]) / (k[i + m + 1] - k[i])
    z1 <- (k[i + m + 2] - x) / (k[i + m + 2] - k[i + 1])
    res <- z0 * bspline(x, k, i, m - 1) + z1 * bspline(x, k, i + 1, m - 1)
  }
  
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

set.seed(1234)
beta_true <- runif(k_dim, -1, 1)
dist_range <- range(dist_matrix, na.rm = TRUE)
fx <- function(x, xlim, beta, df = 10, m = 2){
  B = bspline_basis(x = x, 
                    df = k_dim, 
                    m = m, 
                    xlim = xlim)$X
  XB <- sweep(B, 2, beta, "*")
  return(rowSums(XB))
}