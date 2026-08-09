## -----------------------------------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
.gof_extract <- function(rem_fit) {
 
  X    <- model.matrix(rem_fit)
  n    <- nrow(X)
  mu   <- rem_fit$fitted.values
  eta  <- rem_fit$linear.predictors
  y    <- rem_fit$y             
  sig2 <- rem_fit$sig2
 
  mu.eta <- rem_fit$family$mu.eta(eta)
  Vmu    <- rem_fit$family$variance(mu)
 
  w <- mu.eta * (y - mu) / (sig2 * Vmu)
 
  v <- (mu.eta^2) / (sig2 * Vmu)
  U <- X * w
 
  list(X = X, n = n, mu = mu, eta = eta, y = y, sig2 = sig2,
       mu.eta = mu.eta, Vmu = Vmu, w = w, v = v, U = U)
}

## -----------------------------------------------------------------------------
.gof_info <- function(U, rem_fit) {

  InfoPlusS.inv <- rem_fit$Vp
  if (is.null(InfoPlusS.inv)) {
    stop("rem_fit$Vp is NULL -- refit the gam() call so that Vp is ",
         "available (it is returned by default for essentially any gam fit).")
  }

  InfoPlusS <- tryCatch(solve(InfoPlusS.inv), 
                        error = function(e) MASS::ginv(InfoPlusS.inv))
  U_hat <- colSums(U)
  list(InfoPlusS = InfoPlusS, InfoPlusS.inv = InfoPlusS.inv, U_hat = U_hat)
}

## -----------------------------------------------------------------------------
.gof_term_cols <- function(rem_fit, term) {

  # 1) smooth terms: match `term` against rem_fit$smooth[[k]]$term 
  # (the covariate name the smooth is built on) 
  sm <- rem_fit$smooth
  if (!is.null(sm)) {
    for (s in sm) {
      if (term %in% s$term) return(s$first.para:s$last.para)
    }
  }

  # 2) column matching `term` by name
  # (covers linear terms, where h_{sr,j_l}(v) is just the raw covariate
  X   <- model.matrix(rem_fit)
  idx <- which(colnames(X) == term)
  if (length(idx) == 0) {
    stop("Could not find `term` = '", term, "' among rem_fit's smooth ",
         "terms (rem_fit$smooth[[k]]$term) or parametric columns ",
         "(colnames(model.matrix(rem_fit))).")
  }
  idx
}

## -----------------------------------------------------------------------------
.gof_time_grid <- function(rem_fit, X, term, times, max_grid = 500) {

  n   <- nrow(X)
  idx <- .gof_term_cols(rem_fit, term)
  delta_h <- X[, idx, drop = FALSE]        

  stopifnot(length(times) == n)
  pooled <- times
  Fhat   <- ecdf(times) 
  x.grid <- sort(unique(times))
  
  if (length(x.grid) > max_grid) {
    keep   <- round(seq(1, length(x.grid), length.out = max_grid))
    x.grid <- x.grid[unique(keep)]
  }
  u.grid <- Fhat(x.grid)
  ind_time <- outer(times, x.grid, "<=")

  phi_list <- lapply(seq_len(ncol(delta_h)), function(k) ind_time * delta_h[, k])

  list(x.grid = x.grid, u.grid = u.grid, G.n = length(x.grid),
       phi_list = phi_list, K = ncol(delta_h), xlab = "t  (event time)")
}

## -----------------------------------------------------------------------------
.gof_process_grid <- function(covariate_event, covariate_nonevent, max_grid = 500) {

  pooled <- c(covariate_event, covariate_nonevent)   # N = 2n pooled values
  Fhat   <- ecdf(pooled)                             # empirical CDF, hat F_N

  x.grid <- sort(unique(pooled))
  if (length(x.grid) > max_grid) {
    keep   <- round(seq(1, length(x.grid), length.out = max_grid))
    x.grid <- x.grid[unique(keep)]
  }

  u.grid <- Fhat(x.grid)
  ind_case <- outer(covariate_event,    x.grid, "<=")
  ind_ctrl <- outer(covariate_nonevent, x.grid, "<=")
  phi_diff <- ind_case - ind_ctrl

  list(x.grid = x.grid, u.grid = u.grid, G.n = length(x.grid),
       phi_list = list(phi_diff), K = 1,
       xlab = "u  (empirical quantile of the covariate)")
}

## -----------------------------------------------------------------------------
.gof_covariate_grid <- function(rem_fit, term,
                                 covariate_event, covariate_nonevent,
                                 X_event, X_nonevent,
                                 max_grid = 500) {

  # basis columns for `term` (smooth or parametric), same lookup as .gof_term_cols
  idx <- .gof_term_cols(rem_fit, term)

  h_event    <- X_event[,    idx, drop = FALSE]   # h_{s_i r_i, j_l}(t_i),     n x K
  h_nonevent <- X_nonevent[, idx, drop = FALSE]   # h_{s_i* r_i*, j_l}(t_i),  n x K

  stopifnot(length(covariate_event)    == nrow(h_event),
            length(covariate_nonevent) == nrow(h_nonevent),
            nrow(h_event) == nrow(h_nonevent))

  # pooled covariate values fix the quantile grid / order statistics x^l_(1..N)
  pooled <- c(covariate_event, covariate_nonevent)   # N = 2n
  Fhat   <- ecdf(pooled)

  x.grid <- sort(unique(pooled))
  if (length(x.grid) > max_grid) {
    keep   <- round(seq(1, length(x.grid), length.out = max_grid))
    x.grid <- x.grid[unique(keep)]
  }
  u.grid <- Fhat(x.grid)

  # indicators I(x_{sr}(t) <= x_(floor(Nu)) ), built separately for event/nonevent
  ind_event    <- outer(covariate_event,    x.grid, "<=")   # n x G
  ind_nonevent <- outer(covariate_nonevent, x.grid, "<=")   # n x G

  K <- ncol(h_event)
  # Delta phi^l_i(t_i,u) = I(x_event<=grid) * h_event_{,k}  -  I(x_nonevent<=grid) * h_nonevent_{,k}
  phi_list <- lapply(seq_len(K), function(k) {
    ind_event * h_event[, k] - ind_nonevent * h_nonevent[, k]
  })

  list(x.grid = x.grid, u.grid = u.grid, G.n = length(x.grid),
       phi_list = phi_list, K = K,
       xlab = "u  (empirical quantile of the covariate)")
}

## -----------------------------------------------------------------------------
.gof_bootstrap <- function(Bmat, corr_Gm, w, phi_diff, n, n.sim, Gm) {
 
  raw       <- t(phi_diff) %*% (w * Gm)       
  drift_sim <- t(Bmat) %*% corr_Gm          
  Gstar     <- (raw - drift_sim) / sqrt(n)
 
  list(Gstar = Gstar)
}

## -----------------------------------------------------------------------------
.gof_plot <- function(test, u.grid, norm.obs, norm.sim,
                       n.sim, pvalue, term, main, subtitle, xlab) {

  force(test); force(u.grid); force(norm.obs); force(norm.sim)
  force(n.sim); force(pvalue); force(term); force(main); force(subtitle); force(xlab)

  function() {
    old_par <- par(mar = c(6, 4, 5, 2))
    on.exit(par(old_par))

    n.show <- seq_len(min(20, n.sim))
    yr <- range(norm.obs, norm.sim, finite = TRUE)

    if (!is.null(subtitle)){
      main.title <- main
    } else {
      main.title <- switch(test,
                           functional_form       = "Test for correctly specified functional form",
                           functional_form_score = "Test for correctly specified functional form (score)",
                           time                  = "Test for correct time specification")
    }

    ylab <- switch(test,
      functional_form       = "|Cumulative residual process G(1, u)|",
      functional_form_score = "Cumulative residual process G(1, u) - Norm",
      time                  = "Norm of cumulative residual process G(t,1)")

    xlim <- if (test == "functional_form") c(0, 1) else range(u.grid)

    matplot(u.grid, norm.sim[, n.show, drop = FALSE], type = "l",
            col = "grey80", lty = 1, xlim = xlim, ylim = yr,
            xlab = xlab, ylab = ylab, main = main.title)
    lines(u.grid, norm.obs, col = "black", lwd = 2)

    if (!is.null(subtitle)) mtext(subtitle, side = 3, line = 0.5, cex = 0.9)
    mtext("Observed process vs. simulated realizations under true model",
          side = 1, line = 4, cex = 1, font = 3)

    legend("topright", legend = sprintf("p-value = %.3f", pvalue),
           bty = "n", cex = 1, text.font = 2)
  }
}

## -----------------------------------------------------------------------------
GOF_test <- function(rem_fit,
                      covariate_event = NULL, covariate_nonevent = NULL,
                      X_event = NULL, X_nonevent = NULL,
                      test = c("time",
                               "functional_form",
                               "functional_form_score"),
                      times = NULL, term = NULL,
                      n.sim = 2000, seed = NULL, plot = FALSE,
                      max_grid = 500,
                      main = NULL,
                      subtitle = NULL) {

  test <- match.arg(test)
  if (!is.null(seed)) set.seed(seed)

  q <- .gof_extract(rem_fit)
  X <- q$X; n <- q$n; mu <- q$mu; eta <- q$eta; y <- q$y; sig2 <- q$sig2
  mu.eta <- q$mu.eta; Vmu <- q$Vmu; w <- q$w; v <- q$v; U <- q$U

  info <- .gof_info(U = U, rem_fit = rem_fit)
  InfoPlusS <- info$InfoPlusS; InfoPlusS.inv <- info$InfoPlusS.inv; U_hat <- info$U_hat

  if (test == "functional_form") {
    
    # a(h) = 1 (constant map)
    if (is.null(covariate_event) || is.null(covariate_nonevent)) {
      stop("`covariate_event` and `covariate_nonevent` must be supplied ",
           "when test = 'functional form - identity'.")
    }
    gr <- .gof_process_grid(covariate_event = covariate_event,
                             covariate_nonevent = covariate_nonevent,
                             max_grid = max_grid)

  } else if (test == "functional_form_score") {
    # a(h) = h (identity map)
    if (is.null(covariate_event) || is.null(covariate_nonevent) ||
        is.null(X_event) || is.null(X_nonevent) || is.null(term)) {
      stop("`covariate_event`, `covariate_nonevent`, `X_event`, `X_nonevent`, ",
           "and `term` must all be supplied when test = 'functional form - score'.")
    }
    gr <- .gof_covariate_grid(rem_fit = rem_fit, term = term,
                              covariate_event = covariate_event,
                              covariate_nonevent = covariate_nonevent,
                              X_event = X_event, X_nonevent = X_nonevent,
                              max_grid = max_grid)

  } else {  # "time"
    if (is.null(times) || is.null(term)) {
      stop("`times` and `term` must both be supplied when test = 'time'.")
    }
    gr <- .gof_time_grid(rem_fit = rem_fit, X = X, term = term,
                         times = times, max_grid = max_grid)
  }

  x.grid <- gr$x.grid; u.grid <- gr$u.grid; G.n <- gr$G.n
  phi_list <- gr$phi_list; K <- gr$K; xlab <- gr$xlab

  Bmat_list  <- lapply(phi_list, function(phi) t(mu.eta * X) %*% phi)
  drift_list <- lapply(Bmat_list, function(B) as.vector(t(B) %*% InfoPlusS.inv %*% U_hat))
  Gobs_list  <- Map(function(phi, drift)
    as.vector(w %*% phi) / sqrt(n) - drift / sqrt(n),
    phi_list, drift_list)

  G.obs <- do.call(cbind, Gobs_list)

  Gm <- matrix(rnorm(n * n.sim), n, n.sim)
  correction.i <- U %*% InfoPlusS.inv
  corr_Gm      <- t(correction.i) %*% Gm
  boot_list <- Map(function(phi, Bmat)
    .gof_bootstrap(Bmat = Bmat, corr_Gm = corr_Gm,
                   w = w, phi_diff = phi, n = n,
                   n.sim = n.sim, Gm = Gm),
    phi_list, Bmat_list)
  Gstar_arr <- array(unlist(lapply(boot_list, `[[`, "Gstar")),
                      dim = c(G.n, n.sim, K))

  if (test == "functional_form") {
    norm.obs <- abs(G.obs[, 1])
    norm.sim <- abs(Gstar_arr[, , 1])
  } else {
    norm.obs <- sqrt(rowSums(G.obs^2))
    norm.sim <- apply(Gstar_arr, c(1, 2), function(v) sqrt(sum(v^2)))
  }

  T.obs <- max(norm.obs)
  T.sim <- apply(norm.sim, 2, max)
  pvalue <- mean(T.sim >= T.obs)

  plot_fn <- .gof_plot(test = test, u.grid = u.grid,
                        norm.obs = norm.obs, norm.sim = norm.sim,
                        n.sim = n.sim, pvalue = pvalue,
                        term = term, main = main,
                        subtitle = subtitle, xlab = xlab)

  if (plot) plot_fn()

  list(pvalue = pvalue, T.obs = T.obs, T.sim = T.sim,
       G.obs = G.obs, norm.obs = norm.obs, Gstar_arr = Gstar_arr, K = K,
       x.grid = x.grid, u.grid = u.grid, test = test, term = term,
       info_matrix = "Vp", plot_fn = plot_fn)
}

## -----------------------------------------------------------------------------
run_gof <- function(gam_fits, event_list = NULL, nonevent_list = NULL,
                     X_event_list = NULL, X_nonevent_list = NULL,
                     n_sim, n.sim.boot = 200,
                     test = c("time",
                              "functional_form",
                              "functional_form_score"),
                     times_list = NULL, term = NULL,
                     max_grid = 500, plot = FALSE, main = NULL, subtitle = NULL) {
  test <- match.arg(test)
  sapply(seq_len(n_sim), function(x) {
    GOF_test(gam_fits[[x]],
             covariate_event    = if (!is.null(event_list))     event_list[[x]]     else NULL,
             covariate_nonevent = if (!is.null(nonevent_list))  nonevent_list[[x]]  else NULL,
             X_event    = if (!is.null(X_event_list))    X_event_list[[x]]    else NULL,
             X_nonevent = if (!is.null(X_nonevent_list)) X_nonevent_list[[x]] else NULL,
             test      = test,
             times     = if (!is.null(times_list)) times_list[[x]] else NULL,
             term      = term,
             n.sim     = n.sim.boot,
             seed      = x,
             max_grid  = max_grid,
             plot      = plot,
             main = main,
             subtitle  = subtitle)$pvalue
  })
}