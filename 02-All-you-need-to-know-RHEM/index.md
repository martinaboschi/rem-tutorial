---
layout: default
title: All You Need to Know About RHEM
render_with_liquid: false
---

# All You Need to Know About Relational Hyper-Event Modeling

- [Practical 1: Application to Les Misérables](01-Les-Miserables/)

This tutorial shows how to apply Relational Hyper-Event Modeling (RHEM) to the network of actor co-appearances in Les Misérables, based on 80 actors across 288 chapters with binary gender as a node-level covariate. It covers the full pipeline step by step — from raw data import and time-varying statistic computation in eventnet (individual activity, dyadic activity, gender homophily), through case-control sampling, to model fitting in R using the `amorem` package.
Two inference approaches are presented and compared: conditional logistic regression using 20 non-events per event, and degenerate logistic regression in a case-1-control design, the latter extendable to additive/nonparametric effects via GAM. Both methods consistently find that prior individual and dyadic co-presence increases the rate of future joint appearance, and that co-appearances are more likely among actors of the same gender.



- [Practical 2: Application to Scientific Innovation](02-Scientific-Innovation/)

This tutorial shows how to employ Relational Hyper-Event Modeling (RHEM) to a large-scale dataset composed of the co-authorship and citation network derived from the DBLP/Aminer Citation Network, spanning publications from 1939 to 2023. Each relational hyperevent represents a journal article, defined by a set of authors (sender set) and a set of cited papers (receiver set), forming a two-mode directed hypergraph. A sampled dataset of 10,046 events in case-1-control wide format is used, with covariates capturing prior publication activity, heterogeneity in publication experience, paper outdegree popularity, and paper citation popularity — all constructed via subset-repetition statistics and exponential time-downweighting.
The tutorial first fits a linear RHEM, finding that author teams tend to involve members with heterogeneous publication histories (e.g., a PhD student and their supervisor), and that citation popularity exhibits a "rich gets richer" dynamic. It then progressively extends the model to capture time-varying non-linear effects using GAMs, revealing a non-monotonic relationship between experience heterogeneity and co-authorship rates, as well as a time-varying trend. 
