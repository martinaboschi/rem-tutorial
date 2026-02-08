## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
# install.packages("mgcv")
# install.packages("mgcViz")
# install.packages("RColorBrewer")
# install.packages("dplyr")
library(mgcv)
library(mgcViz)
library(RColorBrewer)
library(dplyr)

## -----------------------------------------------------------------------------
load(file="01-Data/dat_gam.RData")
nrow(dat_gam)

## -----------------------------------------------------------------------------
range(dat_gam$TIME_ev)

## -----------------------------------------------------------------------------
head(dat_gam[, c("SOURCE_ev", "TARGET_ev", "TIME_ev")])

## -----------------------------------------------------------------------------
# set of authors
dat_gam[481, "SOURCE_ev"]

## -----------------------------------------------------------------------------
# set of cited papers
dat_gam[481, "TARGET_ev"]
dat_gam[481, "target.size_ev"]

## -----------------------------------------------------------------------------
# publication time
dat_gam[481, "TIME_ev"]

## -----------------------------------------------------------------------------
head(dat_gam[, c("SOURCE_nv", "TARGET_nv", "TIME_nv")])

## -----------------------------------------------------------------------------
# just a check :)
all(dat_gam[, "source.size_ev"] == dat_gam[, "source.size_nv"])
all(dat_gam[, "target.size_ev"] == dat_gam[, "target.size_nv"])

## -----------------------------------------------------------------------------
# just a check :)
all(dat_gam[, "TIME_nv"] == dat_gam[, "TIME_ev"])

## -----------------------------------------------------------------------------
colnames(dat_gam)[7:19]

## -----------------------------------------------------------------------------
colnames(dat_gam)[26:38]

## -----------------------------------------------------------------------------
dat_gam[,c("diff.author.publication.activity",
                         "avg.author.citation.popularity",
                         "diff.author.citation.popularity",  
                         "paper.outdegree.popularity",
                         "author.sub.rep.1", "author.sub.rep.2", 
                         "collaborate.with.citing.author","reference.sub.rep.1",
                         "reference.sub.rep.2", "author.ref.paper.sub.rep.1.1",
                         "cite.paper.and.its.references.1","author.self.citation")] <-
dat_gam[,c("diff.author.publication.activity_ev",
                         "avg.author.citation.popularity_ev",
                         "diff.author.citation.popularity_ev",  
                         "paper.outdegree.popularity_ev",
                         "author.sub.rep.1_ev", "author.sub.rep.2_ev", 
                         "collaborate.with.citing.author_ev","reference.sub.rep.1_ev",
                         "reference.sub.rep.2_ev", "author.ref.paper.sub.rep.1.1_ev",
                         "cite.paper.and.its.references.1_ev","author.self.citation_ev")] -
dat_gam[,c("diff.author.publication.activity_nv",
                         "avg.author.citation.popularity_nv",
                         "diff.author.citation.popularity_nv",  
                         "paper.outdegree.popularity_nv",
                         "author.sub.rep.1_nv", "author.sub.rep.2_nv", 
                         "collaborate.with.citing.author_nv","reference.sub.rep.1_nv",
                         "reference.sub.rep.2_nv", "author.ref.paper.sub.rep.1.1_nv",
                         "cite.paper.and.its.references.1_nv","author.self.citation_nv")]

## -----------------------------------------------------------------------------
gam_linear <- gam(y ~ -1 + diff.author.publication.activity +
                           paper.outdegree.popularity +
                           author.sub.rep.1 +
                           reference.sub.rep.1,
                           family="binomial",
                           data = dat_gam)

## -----------------------------------------------------------------------------
summary(gam_linear)

## -----------------------------------------------------------------------------
summary(dat_gam$diff.author.publication.activity_ev)

## -----------------------------------------------------------------------------
# just a check :) 
              
# diff.author.publication.activity
all(dat_gam[,42] == dat_gam[,40] - dat_gam[,41])
# avg.author.citation.popularity
all(dat_gam[,45] == dat_gam[,43] - dat_gam[,44])
# diff.author.citation.popularity
all(dat_gam[,48] == dat_gam[,46] - dat_gam[,47])
# paper.outdegree.popularity
all(dat_gam[,51] == dat_gam[,49] - dat_gam[,50])
# author.sub.rep.1
all(dat_gam[,54] == dat_gam[,52] - dat_gam[,53])
# author.sub.rep.2
all(dat_gam[,57] == dat_gam[,55] - dat_gam[,56])
# collaborate.with.citing.author
all(dat_gam[,60] == dat_gam[,58] - dat_gam[,59])
# reference.sub.rep.1
all(dat_gam[,63] == dat_gam[,61] - dat_gam[,62])
# reference.sub.rep.2
all(dat_gam[,66] == dat_gam[,64] - dat_gam[,65])
# author.ref.paper.sub.rep.1.1
all(dat_gam[,69] == dat_gam[,67] - dat_gam[,68])
# cite.paper.and.its.references
all(dat_gam[,72] == dat_gam[,70] - dat_gam[,71])
# author.self.citation
all(dat_gam[,75] == dat_gam[,73] - dat_gam[,74])

## -----------------------------------------------------------------------------
gam_tve <- gam(y ~ -1 + s(TIME_ev, by=diff.author.publication.activity) +
                          paper.outdegree.popularity +
                          author.sub.rep.1 +
                          reference.sub.rep.1,
                          family="binomial",
                          data = dat_gam)

## -----------------------------------------------------------------------------
summary(gam_tve)

## -----------------------------------------------------------------------------
plot(gam_tve)

## -----------------------------------------------------------------------------
diff.author.publication.activity_matrix <-
  cbind(dat_gam$transform_diff.author.publication.activity_ev,
        dat_gam$transform_diff.author.publication.activity_nv)
W <- diff.author.publication.activity_matrix
W[,1] <- 1
W[,2] <- -1
gam_nle <- gam(y ~ -1 + s(diff.author.publication.activity_matrix, by=W) +
                          paper.outdegree.popularity +
                          author.sub.rep.1 +
                          reference.sub.rep.1,
                          family="binomial",
                          data = dat_gam)

## -----------------------------------------------------------------------------
summary(gam_nle)

## -----------------------------------------------------------------------------
plot(gam_nle)

## -----------------------------------------------------------------------------
time_matrix <-
  cbind(dat_gam$transformed_time,
        dat_gam$transformed_time)
gam_tvnle <- gam(y ~ -1 + te(time_matrix, diff.author.publication.activity_matrix, by=W) +
                          paper.outdegree.popularity +
                          author.sub.rep.1 +
                          reference.sub.rep.1,
                          family="binomial",
                          data = dat_gam)
viz <- getViz(gam_tvnle)

## -----------------------------------------------------------------------------
summary(gam_tvnle)

## -----------------------------------------------------------------------------
plot_obj <- plot(viz)
plot_data <- plot_obj$plots[[1]]$ggObj$data
plot_data <- plot_data[!is.na(plot_data$z),]
  plot_data <- plot_data %>%
    group_by(x) %>%
    mutate(z_centered = z - mean(z)) %>%
    ungroup()
ggplot(plot_data, aes(x = x, y = y, fill = z_centered)) +
    geom_tile() +
    geom_contour(mapping = aes(x = x, y = y, z = z_centered, group = 1), 
                 color = "black", inherit.aes = FALSE) +
    scale_fill_viridis_c()

