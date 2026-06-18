## Introduction and Set up

In the two practicals of the workshop *Introduction to Relational Event
Models*, we will show how to use `amore` to simulate relational events
and fit relational event models with increasing levels of complexity.

If you are running this tutorial locally, make sure to uncomment the
`remotes::install_github` command below to install the package.

``` r
# remotes::install_github("franciscorichter/amore")
library(amore)
library(dplyr)
library(survival)
```

This tutorial is inspired by the analysis conducted in (Lerner et al.,
2025), where, for illustrative purposes, RHEM is applied to the network
of actor co-appearances in *Les Misérables*, compiled by Donald Knuth
(Knuth, 1993). The data comprises 80 actors (co-)appearing in one or
several of 288 chapters. The
(binary`\footnote{only for mathematical purposes}`{=tex}) gender of the
actors is known.

## 1. Les Misérables: Data and Covariate Computation

a)  **Given the information above, interpret this data as a relational
    hyper-event network. Define the components of the hyperevent
    sequence, the vertex set, its cardinality, and the available
    attributes.**

-   **80 actors** $|V| = 80$ (co-)appearing in one or several of **288
    chapters** ($n = 288$);
-   $e_i = (S_i, t_i)$ is a **hyperevent** (a chapter) involving a set
    of co-appearing actors $S_i$. In this case, $t_i$ is simply a
    counter for the chapter;
-   $\text{female}(i)$ is a **binary gender attribute** for each actor.

b)  **Open the file `jean_events.csv`. How are the events reported? Does
    it contain only event information?**

The file `jean_events.csv` contains three types of entries (which can be
distinguished using the column `type`), namely `add.actor`, `is.female`,
and `chapter`. Each entry involves only one actor at a time. A chapter
is therefore represented by as many rows as the number of actors
involved.

c)  **Open the configuration file `jean.config.small.txt` in `eventnet`
    and inspect the different windows.**

See the tutorial notebook for more information on how the data is
treated in `eventnet`.

d)  **Load the input data `jean_events_EVENTS.csv` located in `01-Data`.
    Is this a case-1-control dataset? If not, how do you think it can be
    constructed?**

In`eventnet` we computed the necessary statistics to support the
inference techniques that will be explored. Upon inspecting the dataset,
you may notice the presence of missing values (NA) in some columns.

In particular, the `TARGET` column contains only missing values, while
the `SOURCE` column includes multiple entries. This indicates that we
are dealing with `\textbf{undirected relational hyper-events}`{=tex},
where the sender set is a subset of the vertices in a relational
hypergraph $V$, and the receiver set is empty.

Additionally, the `EVENT_INTERVAL_ID` column does not hold interpretable
values, as in this application we just know the order of events. The
`IS_OBSERVED` column is used to distinguish between actual observed
events and non-events (i.e., events that could have occurred but did
not).

``` r
data_original <- read.csv("01-Data/jean_events_EVENTS.csv")
```

``` r
head(data_original)
```

<pre><code>##   IS_OBSERVED     SOURCE TARGET    TYPE EVENT_INTERVAL_ID EVENT INTEGER_TIME
## 1           1 |MY|NP|MB|     NA chapter             1.1.1   109          108
## 2           0 |CL|GR|MT|     NA                             111          110
## 3           0 |MA|BB|CV|     NA                             111          110
## 4           0 |PG|GT|CM|     NA                             111          110
## 5           0 |TH|MM|HL|     NA                             111          110
## 6           0 |TH|LL|SS|     NA                             111          110
##   TIME_POINT TIME_UNIT EVENT_INTERVAL num.actors individual.activity
## 1        109       108              2          3                   0
## 2        111       110              2          3                   0
## 3        111       110              2          3                   0
## 4        111       110              2          3                   0
## 5        111       110              2          3                   0
## 6        111       110              2          3                   0
##   dyadic.activity closure female diff.female
## 1               0       0      1           2
## 2               0       0      2           2
## 3               0       0      0           0
## 4               0       0      0           0
## 5               0       0      1           2
## 6               0       0      2           2
</code></pre>

In this application, the `num.actors` column is also unnecessary. This
is because non-events have been sampled to match the cardinality of the
observed events, meaning that the effect of event size cannot be
meaningfully estimated.

To simplify visualization, we proceed by removing all unnecessary
columns from the dataset.

``` r
data <- data_original[,setdiff(colnames(data_original), 
                               c("TARGET", "TYPE", "EVENT_INTERVAL_ID", "EVENT", 
                                 "INTEGER_TIME", "TIME_POINT", 
                                 "TIME_UNIT", "num.actors"))]
rm(data_original)
head(data)
```

<pre><code>##   IS_OBSERVED     SOURCE EVENT_INTERVAL individual.activity dyadic.activity
## 1           1 |MY|NP|MB|              2                   0               0
## 2           0 |CL|GR|MT|              2                   0               0
## 3           0 |MA|BB|CV|              2                   0               0
## 4           0 |PG|GT|CM|              2                   0               0
## 5           0 |TH|MM|HL|              2                   0               0
## 6           0 |TH|LL|SS|              2                   0               0
##   closure female diff.female
## 1       0      1           2
## 2       0      2           2
## 3       0      0           0
## 4       0      0           0
## 5       0      1           2
## 6       0      2           2
</code></pre>

The object `data` is not a case-1-control dataset. To construct it, we
first create two separate objects for the datasets containing events and
non-events, respectively. We then merge these datasets using a set of
variables that identify which event each non-event refers to. Among the
non-events satisfying this condition, we sample only one.

``` r
data_ev <- data %>% filter(IS_OBSERVED == 1)
data_nv <- data %>% filter(IS_OBSERVED == 0)
```

``` r
# Add a new column to label each row in 'data_ev' and 'data_nv' 
# as event data and non-event data respectively
data_ev_tagged <- data_ev %>%
  mutate(.row_type = "ev")
data_nv_tagged <- data_nv %>%
  mutate(.row_type = "nv")
head(data_ev_tagged)
```

<pre><code>##   IS_OBSERVED              SOURCE EVENT_INTERVAL individual.activity
## 1           1          |MY|NP|MB|              2                   0
## 2           1          |MY|ME|MB|              3                   2
## 3           1                |MY|              4                   2
## 4           1 |MY|ME|CL|GE|MC|MB|              5                   6
## 5           1          |MY|MB|ME|              6                   9
## 6           1             |ME|MY|              7                   8
##   dyadic.activity closure female diff.female .row_type
## 1               0       0      1           2        ev
## 2               1       2      2           2        ev
## 3               0       0      0           0        ev
## 4               4       8      3           9        ev
## 5               7      32      2           2        ev
## 6               3      12      1           1        ev
</code></pre>

``` r
# head(data_nv_tagged)
rm(data_ev)
rm(data_nv)
```

``` r
set.seed(1234)
merge_id_cols <- c("EVENT_INTERVAL")

# For each group defined by the merge_id_cols, take one random non-event
data_nv_1_sampled <- data_nv_tagged %>%
  group_by(across(all_of(merge_id_cols))) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Perform a left join: 
# for each row in data_ev, attach the corresponding non-event row (if available)
# based on matching values in the merge_id_cols
# The suffixes "_ev" and "_nv" will be added to columns from data_ev and data_nv 
dat_gam_1 <- data_ev_tagged %>%
  left_join(data_nv_1_sampled, 
            by = merge_id_cols, suffix = c("_ev", "_nv"))
rm(data_nv_1_sampled)
```

## 2. Model fitting: conditional and degenerate logistic regression

a)  **Based on the data structures you were given and those you
    constructed, which inference techniques can you apply?**

**OPTION 1**

Since we do not have access to timing information for the events and we
have a case-control sampling with 20 non-events, we use the `clogit`
function in R to fit our model (conditional logistic regression). When
applied in this context, `clogit` internally calls the `coxph` routine.

``` r
clogit_fit <- clogit(IS_OBSERVED ~ 
                         + diff.female
                         + female
                         + individual.activity 
                         + dyadic.activity 
                         + strata(EVENT_INTERVAL)
                         , data = data)
summary(clogit_fit)
```

<pre><code>## Call:
## coxph(formula = Surv(rep(1, 6048L), IS_OBSERVED) ~ +diff.female + 
##     female + individual.activity + dyadic.activity + strata(EVENT_INTERVAL), 
##     data = data, method = "exact")
## 
##   n= 6048, number of events= 288 
## 
##                          coef exp(coef)  se(coef)      z Pr(>|z|)    
## diff.female         -0.290310  0.748032  0.117345 -2.474   0.0134 *  
## female              -0.293382  0.745737  0.144822 -2.026   0.0428 *  
## individual.activity  0.042973  1.043909  0.003918 10.967   <2e-16 ***
## dyadic.activity      0.460395  1.584701  0.049289  9.341   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##                     exp(coef) exp(-coef) lower .95 upper .95
## diff.female            0.7480     1.3368    0.5943    0.9415
## female                 0.7457     1.3410    0.5615    0.9905
## individual.activity    1.0439     0.9579    1.0359    1.0520
## dyadic.activity        1.5847     0.6310    1.4388    1.7454
## 
## Concordance= 0.896  (se = 0.013 )
## Likelihood ratio test= 819.9  on 4 df,   p=<2e-16
## Wald test            = 245.8  on 4 df,   p=<2e-16
## Score (logrank) test = 1418  on 4 df,   p=<2e-16
</code></pre>

We find a
`\textcolor{myred}{positive effect of individual, dyadic}`{=tex}. This
suggests that prior (co-)presence in previous chapters - whether alone
or in pairs - increases the rate of participating in events. A negative
effect for diff.female suggests a positive effect for gender homophily,
favouring co-appearence with other actors of the same gender.

**OPTION 2**

Using the case-1-control dataset, it is also possible to fit a
degenerate logistic regression without an intercept, where the
covariates are defined as the differences between the covariates for the
events and those for the non-events, and the response variable is
constant and equal to 1.

``` r
# covariates defined as difference between 
# covariate values for event and corresponding non-event
dat_gam_1$female <- 
  dat_gam_1$female_ev - dat_gam_1$female_nv
dat_gam_1$diff_female <- 
  dat_gam_1$diff.female_ev - dat_gam_1$diff.female_nv
dat_gam_1$individual_activity <- 
  dat_gam_1$individual.activity_ev - dat_gam_1$individual.activity_nv
dat_gam_1$dyadic_activity <- 
  dat_gam_1$dyadic.activity_ev - dat_gam_1$dyadic.activity_nv
dat_gam_1$closure <- 
  dat_gam_1$closure_ev - dat_gam_1$closure_nv

# constant response equal to 1
dat_gam_1$y <- 1
```

``` r
gam_fit <- glm(y ~ 
                + diff_female
                + female
                + individual_activity 
                + dyadic_activity 
                - 1 # no intercept
               , data = dat_gam_1, 
               family="binomial")
```

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

``` r
summary(gam_fit)
```

<pre><code>## 
## Call:
## glm(formula = y ~ +diff_female + female + individual_activity + 
##     dyadic_activity - 1, family = "binomial", data = dat_gam_1)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## 0.00000  0.00002  0.05696  0.69837  2.94123  
## 
## Coefficients:
##                     Estimate Std. Error z value Pr(>|z|)    
## diff_female         -0.67999    0.37963  -1.791 0.073260 .  
## female              -0.71924    0.31532  -2.281 0.022550 *  
## individual_activity  0.07649    0.01704   4.490 7.13e-06 ***
## dyadic_activity      1.81184    0.54181   3.344 0.000826 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 399.25  on 288  degrees of freedom
## Residual deviance: 121.75  on 284  degrees of freedom
## AIC: 129.75
## 
## Number of Fisher Scoring iterations: 11
</code></pre>

As before, we find a positive effect of individual, dyadic and a
negative main effect for female gender and a positive effect for gender
homophily.

## References

-   Knuth, D. E. (1993). The Stanford GraphBase: a platform for
    combinatorial computing (Vol. 1). New York: AcM Press.

-   Lerner, J., Hâncean, M. G., & Perc, M. (2025). Modeling temporal
    hypergraphs. Journal of Complex Networks, 13(6).

-   Lerner, J., & Lomi, A. (2020). Reliability of relational event model
    estimates under sampling: How to fit a relational event model to 360
    million dyadic events. Network Science, 8(1), 97-135.

-   Lerner, J., & Lomi, A. (2023). Relational hyperevent models for
    polyadic interaction networks. Journal of the Royal Statistical
    Society Series A: Statistics in Society, 186(3), 577-600.
