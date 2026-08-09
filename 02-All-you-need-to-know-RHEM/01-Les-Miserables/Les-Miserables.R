## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)
knitr::knit_hooks$set(output = function(x, options) {  
  paste0('<pre><code>', x, '</code></pre>')})

## ----message=FALSE, warning=FALSE---------------------------------------------
# remotes::install_github("franciscorichter/amorem")
library(amorem)

## ----message=FALSE, warning=FALSE---------------------------------------------
# install.packages("dplyr")
library(dplyr)

## -----------------------------------------------------------------------------
data_original <- read.csv("00-eventnet/jean_events.csv")
head(data_original)

## -----------------------------------------------------------------------------
table(data_original$type)

## -----------------------------------------------------------------------------
data_preproc <- read.csv("01-Data/jean_events_EVENTS.csv")
head(data_preproc)

## -----------------------------------------------------------------------------
raw_data_m20 <- data_preproc[,setdiff(colnames(data_preproc),
                              c("TARGET", "TYPE", "EVENT_INTERVAL_ID", "EVENT",
                                "INTEGER_TIME", "TIME_POINT",
                                "TIME_UNIT", "num.actors"))]
rm(data_preproc)
head(raw_data_m20)

## -----------------------------------------------------------------------------
data_ev <- raw_data_m20 %>% filter(IS_OBSERVED == 1)
data_nv <- raw_data_m20 %>% filter(IS_OBSERVED == 0)
set.seed(1234)
merge_id_cols <- c("EVENT_INTERVAL")
# For each group defined by the merge_id_cols, take one random non-event
data_nv_sampled <- data_nv %>% 
  group_by(across(all_of(merge_id_cols))) %>% 
  slice_sample(n = 1) %>%
  ungroup()
raw_data_m1 <- bind_rows(data_ev, data_nv_sampled) %>%
  arrange(EVENT_INTERVAL)

## -----------------------------------------------------------------------------
ncc_data <- widen_case_control(raw_data_m1, 
                               case = "IS_OBSERVED", 
                               stratum = "EVENT_INTERVAL")

## -----------------------------------------------------------------------------
fit_clogit <- rem(IS_OBSERVED ~ diff.female
                  + female + individual.activity 
                  + dyadic.activity, 
                  method = "clogit", 
                  data = raw_data_m20, 
                  stratum = "EVENT_INTERVAL")
summary(fit_clogit)

## -----------------------------------------------------------------------------
fit_glm <- rem(~ d_diff.female + d_female 
               + d_individual.activity 
               + d_dyadic.activity 
               - 1, data = ncc_data, 
               method="gam")
summary(fit_glm)

