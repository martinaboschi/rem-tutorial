## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
# install.packages("mgcv")
# install.packages("RColorBrewer")
# install.packages("survival")
# install.packages("dplyr")
library(mgcv)
library(RColorBrewer)
library(survival)
library(dplyr)

## ---- echo=FALSE--------------------------------------------------------------
pal.blue <- brewer.pal(9, "Blues")
pal.rose <- brewer.pal(9, "RdPu")
colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#999999", "#66C2A5", "#FC8D62")

## -----------------------------------------------------------------------------
data_original <- read.csv("01-Data/jean_events_EVENTS.csv")

## -----------------------------------------------------------------------------
head(data_original)

## -----------------------------------------------------------------------------
data <- data_original[,setdiff(colnames(data_original), 
                               c("TARGET", "TYPE", "EVENT_INTERVAL_ID", "EVENT", 
                                 "INTEGER_TIME", "TIME_POINT", 
                                 "TIME_UNIT", "num.actors"))]
rm(data_original)
head(data)

## -----------------------------------------------------------------------------
data_ev <- data %>% filter(IS_OBSERVED == 1)
data_nv <- data %>% filter(IS_OBSERVED == 0)

## -----------------------------------------------------------------------------
# Add a new column to label each row in 'data_ev' and 'data_nv' 
# as event data and non-event data respectively
data_ev_tagged <- data_ev %>%
  mutate(.row_type = "ev")
data_nv_tagged <- data_nv %>%
  mutate(.row_type = "nv")
head(data_ev_tagged)
# head(data_nv_tagged)
rm(data_ev)
rm(data_nv)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
clogit_fit <- clogit(IS_OBSERVED ~ 
                         + diff.female
                         + female
                         + individual.activity 
                         + dyadic.activity 
                         + strata(EVENT_INTERVAL)
                         , data = data)
summary(clogit_fit)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
gam_fit <- glm(y ~ 
                + diff_female
                + female
                + individual_activity 
                + dyadic_activity 
                - 1 # no intercept
               , data = dat_gam_1, 
               family="binomial")
summary(gam_fit)

