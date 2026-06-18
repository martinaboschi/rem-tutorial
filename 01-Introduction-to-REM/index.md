---
layout: default
title: Introduction to REMs
render_with_liquid: false
---

# Introduction to Relational Event Models

- [Practical 1: Simulating Relational Events](01-Simulating-REM-Data/)

This tutorial introduces Relational Event Models (REMs) through simulation and inference using the `amorem` R package. Starting from the core representation of relational event data as a sequence of sender-receiver-time triplets, we show how to simulate event sequences governed by endogenous network statistics — specifically reciprocity — under a model with a fixed linear effect. We then walk through two inference frameworks: sampled partial likelihood estimation via conditional logistic regression, and case-1-control partial likelihood via degenerate logistic regression, explaining the data structures each requires. A simulation study with 100 replications compares these approaches across varying numbers of sampled non-events (1, 7, and 20 controls per event), illustrating the bias-variance tradeoff: reducing the control set increases estimation variance while preserving consistency.


- [Practical 2: Application to Alien Species Invasions](02-Alien-Species-Invasions/)

This tutorial demonstrates the application REMs beyond linear effect specifications, using data on the first recorded invasions of insect species into world regions between 1880 and 2005. Framing alien species invasions as relational events — where species act as senders and regions as receivers — we show how to prepare ecological data for REM inference using the `amorem` R package, including risk set construction that accounts for native ranges and the non-recurrent nature of first records. Three endogenous, time-dependent covariates are derived by combining event history with exogenous data sources: climatic dissimilarity (minimum temperature difference between a target region and previously invaded regions), international trade flows, and geographic distance. Using case-1-control sampling and the degenerate logistic regression framework, we fit and interpret four model specifications — linear, time-varying, nonlinear, and random-intercept effects — illustrating how each captures a distinct aspect of invasion dynamics. 