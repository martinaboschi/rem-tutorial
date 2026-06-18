## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::knit_hooks$set(purl = knitr::hook_purl)

## ----warning = FALSE----------------------------------------------------------
# install.packages("remotes")
# install.packages("pbapply")
# remotes::install_github("franciscorichter/amore")

library(amore)
library(pbapply)

## -----------------------------------------------------------------------------
load(file="01-Data/01-Inputs/input00.RData")

head(FR[,c("LifeForm", "Taxon",
           "Region", "FirstRecord",
           "Source")])

## -----------------------------------------------------------------------------
native[sample(1:nrow(native), 10), c("species", "region")]

## -----------------------------------------------------------------------------
head(first_records[,c("year","lf",
           "species", 
           "region")])

## -----------------------------------------------------------------------------
event_log <- standardize_event_log(first_records[,c(1,3:4)],
                                   sender_col = "species",
                                   receiver_col = "region",
                                   time_col = "year")
  
head(event_log)

## -----------------------------------------------------------------------------
native_df<-native[, c("species", "region")]
colnames(native_df)<-c("sender", "receiver")
head(native_df)

## -----------------------------------------------------------------------------
set.seed(1234)
cc <- sample_non_events(event_log,
                        n_controls = 1,              # 1 control is default
                        scope         = "all",
                        mode          = "two",
                        risk          = "remove",          # drops past + concurrent
                        exclude_pairs = native_df[,1:2])   ### excludes native dyads


head(cc)

