---
layout: default
title: "Practical 1: Les Misérables"
render_with_liquid: false
---

# Practical 1: Les Misérables

## Introduction

*Write your introduction here.*

---

## 1.

1)  **Given the information above, interpret this data as a relational
    hyper-event network. Define the components of the hyperevent
    sequence, the vertex set, its cardinality, and the available
    attributes.**

- **80 actors** $`|V| = 80`$ (co-)appearing in one or several of **288
  chapters** ($`n = 288`$);
- $`e_i = (S_i, t_i)`$ is a **hyperevent** (a chapter) involving a set
  of co-appearing actors $`S_i`$. In this case, $`t_i`$ is simply a
  counter for the chapter;
- $`\text{female}(i)`$ is a **binary gender attribute** for each actor.

2)  **Open the file `jean_events.csv`. How are the events reported? Does
    it contain only event information?**

The file `jean_events.csv` contains three types of entries (which can be
distinguished using the column `type`), namely `add.actor`, `is.female`,
and `chapter`. Each entry involves only one actor at a time. A chapter
is therefore represented by as many rows as the number of actors
involved.

3)  **Open the configuration file `jean.config.small.txt` in `eventnet`
    and inspect the different windows.**

See the tutorial notebook for more information on how the data is
treated in `eventnet`.

4)  **Load the input data `jean_events_EVENTS.csv` located in `01-Data`.
    Is this a case-1-control dataset? If not, how do you think it can be
    constructed?**

In`eventnet` we computed the necessary statistics to support the
inference techniques that will be explored. Upon inspecting the dataset,
you may notice the presence of missing values (NA) in some columns.

In particular, the `TARGET` column contains only missing values, while
the `SOURCE` column includes multiple entries. This indicates that we
are dealing with , where the sender set is a subset of the vertices in a
relational hypergraph $`V`$, and the receiver set is empty.

Additionally, the `EVENT_INTERVAL_ID` column does not hold interpretable
values, as in this application we just know the order of events. The
`IS_OBSERVED` column is used to distinguish between actual observed
events and non-events (i.e., events that could have occurred but did
not).

``` r
data_original <- read.csv("01-Data/jean_events_EVENTS.csv")
```

test.R

``` r
head(data_original)
```

    ##   IS_OBSERVED     SOURCE TARGET    TYPE EVENT_INTERVAL_ID EVENT INTEGER_TIME
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

test.R

---

## 2.

*Content here.*

---

## 3.

*Content here.*

---

## 4.

*Content here.*

---

## References