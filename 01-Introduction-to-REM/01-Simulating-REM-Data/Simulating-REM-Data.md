## Introduction

In the two practicals of the workshop *Introduction to Relational Event
Models*, we will show how to use `amore` to simulate relational events
and fit relational event models with increasing levels of complexity.

If you are running this tutorial locally, make sure to uncomment the
`remotes::install_github` command below to install the package.

``` r
# remotes::install_github("franciscorichter/amore")
library(amore)
```

While exploring the main functions of the package, we will cover how:

-   to simulate relational event data in `amore` \[Theoretical
    reference\]
-   to structure data for inference \[Theoretical reference\]
-   to fit models with linear effects \[Theoretical reference\]
-   different inference techniques compare \[Theoretical reference\]

## 1. How to simulate relational event data in `amore`

As seen in \[Theoretical reference\], *relational event data* consists
of a sequence of $n$ triplets:

$$E = \{e_i = (s_i, r_i, t_i), i=1,\dots,n\}$$

where $s_i\in V^S$ and $r_i\in V^R$. $V^S$ and $V^R$ are the *sender
set* and the *receiver set*, respectively.

The *number of events* (`n_events`), the *sender set* (`senders`), and
the *receiver set* (`receivers`) are the first three arguments that we
provide to simulate our relational event sequence by means of the
function `simulate_relational_events`.

We saw in \[Theoretical reference\] that relational events can be both
described and simulated using a *rate function*. This rate function can
be decomposed into three main blocks:

1.  *Risk indicator*: equal to 1 if the dyad $(s,r)$ is at risk at time
    $t$. In this first practical, we assume that all potential dyads -
    except *self-loops* (by default, `allow_loops = FALSE`) - can be
    involved in an event. Please refer to the second practical of this
    tutorial for circumstances in which this is not the case.

2.  *Baseline hazard*: a time-varying function expressing the *global
    determinants* of the event rate. In its current state, `amore`
    supports a **time-invariant** baseline function, which can be
    specified via the `baseline_rate` argument using a positive scalar.

3.  *Event log-rate contribution*: expresses how *dyadic covariates*
    impact the event log-rate. In this first practical we focus on a
    *linear contribution*. To simulate events following a linear REM,
    the covariates of interest must be specified in the
    `endogenous_stats` argument (chosen among the endogenous mechanisms
    listed in the corresponding documentation). Each statistic is then
    associated with an *effect* in the `endogenous_effects` argument.

Whenever the true parameters are unknown, inference is required. This
requires not only *events* - i.e. interactions that have been observed -
but also *non-events*, i.e. interactions involving dyads at risk at that
time (i.e. included in the risk set). Depending on the inference method,
a different number of non-events may be needed. The
`simulate_relational_events` function allows you to specify this via the
`n_controls` argument, i.e. the number of non-events to be stored for
each event.

``` r
set.seed(1)
raw_data_m1 <- simulate_relational_events(
  n_events           = 1000,
  senders            = paste0("a", 1:20),
  receivers          = paste0("a", 1:20),
  baseline_rate      = 1,
  endogenous_stats   = "reciprocity_count", 
  endogenous_effects = c(reciprocity_count = 0.6),
  n_controls         = 1
)
```

``` r
head(raw_data_m1,16)
```

<pre><code>##    stratum event sender receiver        time reciprocity_count
## 1        1     1     a3       a8 0.001987321                 0
## 2        1     0    a17       a9 0.001987321                 0
## 3        2     1    a20      a18 0.002354409                 0
## 4        2     0    a14      a16 0.002354409                 0
## 5        3     1    a11      a20 0.003031029                 0
## 6        3     0    a17      a10 0.003031029                 0
## 7        4     1     a2      a14 0.007734052                 0
## 8        4     0    a10       a5 0.007734052                 0
## 9        5     1    a14       a5 0.009142085                 0
## 10       5     0     a8      a18 0.009142085                 0
## 11       6     1    a16      a18 0.010586697                 0
## 12       6     0     a3       a5 0.010586697                 0
## 13       7     1     a5       a3 0.011374750                 0
## 14       7     0    a12       a6 0.011374750                 0
## 15       8     1    a20      a11 0.014583015                 1
## 16       8     0    a14       a9 0.014583015                 0
</code></pre>

## 2. How to structure data for inference

In \[Theoretical reference\], we covered a wide range of inference
methods. Since most of the REM literature you are likely to encounter
relies on *(sampled) partial likelihood estimation*, we start with a
practical example of this approach. We then provide an example using the
*case-1-control partial likelihood*, which opens the door to subsequent
practicals where we move beyond linearity.

Each approach requires a particular data structure. This tutorial
distinguishes between two main data formats:

1.  **Long format**: each event and its associated non-events are stored
    as separate rows. An `event` variable distinguishes events from
    non-events, and all rows referring to the same event share the same
    value of `stratum`.

2.  **Wide format**: each non-event is paired with its corresponding
    event, so the information about both can be stored on the same row.

### 2.1. Data for Sampled partial likelihood

We simulate the events again, this time storing $7$ non-events per event
in the long format, which is the default output format of
`simulate_relational_events`.

``` r
set.seed(1)
raw_data_m7 <- simulate_relational_events(
  n_events           = 1000,
  senders            = paste0("a", 1:20),
  receivers          = paste0("a", 1:20),
  baseline_rate      = 1,
  n_controls         = 7,
  endogenous_stats   = "reciprocity_count", 
  endogenous_effects = c(reciprocity_count = 0.6),
)
```

``` r
head(raw_data_m7,10)
```

<pre><code>##    stratum event sender receiver        time reciprocity_count
## 1        1     1     a3       a8 0.001987321                 0
## 2        1     0    a17       a9 0.001987321                 0
## 3        1     0    a16       a7 0.001987321                 0
## 4        1     0    a15      a16 0.001987321                 0
## 5        1     0     a5      a15 0.001987321                 0
## 6        1     0    a18      a10 0.001987321                 0
## 7        1     0     a4      a17 0.001987321                 0
## 8        1     0    a10       a5 0.001987321                 0
## 9        2     1    a15       a5 0.003404473                 0
## 10       2     0     a8      a18 0.003404473                 0
</code></pre>

### 2.2. Data for Case-1-control partial likelihood

As discussed in \[Theoretical reference\], when there is only one
non-event per event, a *logistic regression* can be fitted. In this
case, the statistical unit is the event/non-event pair. The logistic
regression is fitted with a fixed response equal to 1 (which must
therefore be included in the data) and covariates equal to the
differences between the covariates of the event and those of the
non-event. These differences must also be recorded in the data.

The function `widen_case_control` transforms data from long format to
wide format. To do so, we need to specify the `event` argument - the
dummy variable that distinguishes events from non-events - and the
`stratum` argument, namely the variable that identifies which non-event
is associated with each event.

``` r
wide_data_m1 <- widen_case_control(raw_data_m1, case = "event", stratum= "stratum")
head(wide_data_m1)
```

<pre><code>##   stratum sender_ev sender_nv receiver_ev receiver_nv     time_ev     time_nv
## 1       1        a3       a17          a8          a9 0.001987321 0.001987321
## 2      10       a15        a4         a16         a16 0.022153557 0.022153557
## 3     100        a3        a2         a15          a1 0.207398477 0.207398477
## 4    1000       a18        a8         a20          a2 0.565565424 0.565565424
## 5     101        a2       a12         a17         a19 0.210725053 0.210725053
## 6     102       a11       a19         a13         a11 0.212042188 0.212042188
##   d_time reciprocity_count_ev reciprocity_count_nv d_reciprocity_count
## 1      0                    0                    0                   0
## 2      0                    0                    0                   0
## 3      0                    1                    0                   1
## 4      0                  336                    2                 334
## 5      0                    0                    0                   0
## 6      0                    0                    1                  -1
</code></pre>

**Remark**: Time.

## 3. How to fit models with linear effects

Once the data has been structured, we are ready to fit the model. We do
that by using the `rem` function in `amore`. The desired inference
method can be selected via the `method` argument.

### 3.1. Sampled partial likelihood via Conditional logistic regression

As discussed in \[Theoretical reference\], when considering more than
one non-event per event, it can be shown that the partial likelihood
coincides with that of a *conditional logistic regression* - hence we
set the method to `clogit`.

``` r
fit_clogit <- rem(event ~ reciprocity_count, 
                  data = raw_data_m7, 
                  method = "clogit")
summary(fit_clogit)
```

<pre><code>## Call:
## coxph(formula = Surv(rep(1, 8000L), .case) ~ reciprocity_count + 
##     strata(.strat), data = cl, method = "exact")
## 
##   n= 8000, number of events= 1000 
## 
##                      coef exp(coef) se(coef)     z Pr(>|z|)    
## reciprocity_count 0.66848   1.95126  0.08316 8.039 9.07e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##                   exp(coef) exp(-coef) lower .95 upper .95
## reciprocity_count     1.951     0.5125     1.658     2.297
## 
## Concordance= 0.898  (se = 0.066 )
## Likelihood ratio test= 517.3  on 1 df,   p=<2e-16
## Wald test            = 64.62  on 1 df,   p=9e-16
## Score (logrank) test = 798.1  on 1 df,   p=<2e-16
</code></pre>

**Interpretation:** the estimated coefficient for `reciprocity_count` is
0.669 and is highly significant. This means that the rate is multiplied
by $\exp(0.669) \approx 1.95$; that is, a dyad with one prior
reciprocated interaction is roughly twice as likely to interact next
compared to a dyad with none. This was expected since the true value
employed in simulation was 0.6.

### 3.2. Case-1-control partial likelihood via Degenerate logistic regression

When considering only one non-event per event (\[Theoretical
reference\]), it can be shown that the partial likelihood coincides with
that of a *degenerate logistic regression*:

-   without an intercept term,
-   with a constant response equal to 1,
-   with covariates equal to the difference between the covariates of
    the event and those of the sampled non-event.

This formulation allows the framework to be extended to additive
effects, which is why the method here is `gam`.

``` r
fit_glm <- rem(~ reciprocity_count, 
               data = wide_data_m1, 
               method = "gam")
summary(fit_glm)
```

<pre><code>## 
## Family: binomial 
## Link function: logit 
## 
## Formula:
## one ~ -1 + reciprocity_count
## 
## Parametric coefficients:
##                   Estimate Std. Error z value Pr(>|z|)    
## reciprocity_count  0.54154    0.09815   5.518 3.44e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## 
## R-sq.(adj) =   -Inf   Deviance explained = -Inf%
## UBRE = -0.55893  Scale est. = 1         n = 1000
</code></pre>

**Interpretation:** the estimated coefficient for `reciprocity_count`
(0.542 in this case) is interpreted exactly as before.

Why might we need more than one non-event per event, then?

## 4. How different inference techniques compare

To compare different inference methods - and specifically different
numbers of non-events per event - we run a simulation study `N_SIM`
times, with 100 replications. All other simulation parameters remain as
before.

``` r
# parameters
N_SIM        <- 100
N_EVENTS     <- 1000
TRUE_BETA    <- 0.6
SENDERS      <- paste0("a", 1:20)
RECEIVERS    <- paste0("a", 1:20)
```

``` r
# storage
coefs <- data.frame(
  glm_1   = numeric(N_SIM),   # degenerate logistic, 1 control
  clogit_7  = numeric(N_SIM), # conditional logit,   7 controls
  clogit_20 = numeric(N_SIM)  # conditional logit,  20 controls
)
```

For each iteration:

-   we store the case-1-control dataset, compute the corresponding
    wide-format dataset, and fit a degenerate logistic regression;
-   we store case-7-control and case-20-control datasets and fit a
    conditional logistic regression to each;
-   we store the estimated coefficients from the three models.

``` r
for (i in seq_len(N_SIM)) {

  set.seed(i)

  # case-1-control dataset
  d1 <- simulate_relational_events(
    n_events           = N_EVENTS,
    senders            = SENDERS,
    receivers          = RECEIVERS,
    baseline_rate      = 1,
    n_controls         = 1,
    endogenous_stats   = "reciprocity_count",
    endogenous_effects = c(reciprocity_count = TRUE_BETA),
  )

  # fit degenerate logistic regression
  wide_d1 <- widen_case_control(d1, case = "event", stratum = "stratum")
  fit_glm <- rem(~ reciprocity_count, data = wide_d1, method = "gam")
  coefs$glm_1[i] <- coef(fit_glm)[["reciprocity_count"]]


  # case-7-control dataset
  d7 <- simulate_relational_events(
    n_events           = N_EVENTS,
    senders            = SENDERS,
    receivers          = RECEIVERS,
    baseline_rate      = 1,
    n_controls         = 7,
    endogenous_stats   = "reciprocity_count",
    endogenous_effects = c(reciprocity_count = TRUE_BETA)
  )

  # fit conditional logistic regression
  fit_clogit_7 <- rem(event ~ reciprocity_count, data = d7, method = "clogit")
  coefs$clogit_7[i] <- coef(fit_clogit_7)[["reciprocity_count"]]

  # case-20-control dataset
  d20 <- simulate_relational_events(
    n_events           = N_EVENTS,
    senders            = SENDERS,
    receivers          = RECEIVERS,
    baseline_rate      = 1,
    n_controls         = 20,
    endogenous_stats   = "reciprocity_count",
    endogenous_effects = c(reciprocity_count = TRUE_BETA),
  )

  # fit conditional logistic regression
  fit_clogit_20 <- rem(event ~ reciprocity_count, data = d20, method = "clogit")
  coefs$clogit_20[i] <- coef(fit_clogit_20)[["reciprocity_count"]]

  if (i %% 10 == 0) message(sprintf("  Completed %d / %d replications", i, N_SIM))


}
```

    ##   Completed 10 / 100 replications

    ##   Completed 20 / 100 replications

    ##   Completed 30 / 100 replications

    ##   Completed 40 / 100 replications

    ##   Completed 50 / 100 replications

    ##   Completed 60 / 100 replications

    ##   Completed 70 / 100 replications

    ##   Completed 80 / 100 replications

    ##   Completed 90 / 100 replications

    ##   Completed 100 / 100 replications

We then plot the coefficients across iterations for the three models.

``` r
coefs_long <- stack(coefs)
levels(coefs_long$ind) <- c(
  "Logistic\n(1 control)",
  "Cond. logit\n(7 controls)",
  "Cond. logit\n(20 controls)"
)
cols <- c("#AEC6CF", "#E8D57A", "#FF9999")
```

``` r
boxplot(
  values ~ ind,
  data       = coefs_long,
  col        = cols,
  border     = darken <- adjustcolor(cols, red.f = 0.6, green.f = 0.6, blue.f = 0.6),
  outcol     = cols,
  outpch     = 16,
  cex        = 0.6,
  las        = 1,
  xlab       = "",
  ylab       = expression(hat(beta)),
  main       = expression(
    "Estimated "*beta*" for reciprocity_count across 100 replications"
  ),
  cex.main   = 1.0,
  cex.lab    = 0.95,
  frame.plot = FALSE
)
abline(h = TRUE_BETA, lty = 2, lwd = 1.8, col = "black")
```

![](Simulating-REM-Data_files/figure-markdown/simulation_varying_m-1.png)

As expected, reducing the number of controls **increases the variance**
of the estimates, while the estimator remains **consistent** - all three
distributions are centered around the true value of 0.6, but the
case-1-control logistic regression shows noticeably more spread than the
conditional logistic regression with 7 or 20 controls.

However, this additional variance comes with a key advantage: the
possibility to go **beyond linear effects**, which we will explore in
the next practical.
