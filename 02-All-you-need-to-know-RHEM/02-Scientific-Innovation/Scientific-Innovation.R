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

## ----message=FALSE, warning=FALSE---------------------------------------------
# install.packages("mgcViz")
library(mgcViz)

## -----------------------------------------------------------------------------
load(file="01-Data/dat_gam.RData")
nrow(dat_gam)

## -----------------------------------------------------------------------------
range(dat_gam$TIME_ev)

## -----------------------------------------------------------------------------
ncol(dat_gam)

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
gam_linear <- rem(~ diff.author.publication.activity +
                    paper.outdegree.popularity +
                    author.sub.rep.1 +
                    reference.sub.rep.1,
                  method = "gam",
                  data = dat_gam)
summary(gam_linear)

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
gam_tve <- gam_tve <- rem(~ tv(diff.author.publication.activity) +
                 paper.outdegree.popularity +
                 author.sub.rep.1 +
                 reference.sub.rep.1,
                 method ="gam", time = "TIME_ev",
                 data = dat_gam)

## -----------------------------------------------------------------------------
summary(gam_tve)

## ----gam_tve------------------------------------------------------------------
plot(gam_tve)

## -----------------------------------------------------------------------------
gam_tve_transformed <- rem(~ tv(diff.author.publication.activity) +
                           paper.outdegree.popularity +
                           author.sub.rep.1 +
                           reference.sub.rep.1,
                           method ="gam", time = "transformed_time",
                           data = dat_gam)

## ----gam_tve_transformed------------------------------------------------------
plot(gam_tve_transformed)

## -----------------------------------------------------------------------------
gam_nle <- rem(~ nl(diff.author.publication.activity) 
                 # automatically uses the transformed one when computing nl
               
                 + paper.outdegree.popularity
                 + author.sub.rep.1
                 + reference.sub.rep.1,
                 method ="gam",
                 data = dat_gam)

## -----------------------------------------------------------------------------
summary(gam_nle)

## ----gam_nle------------------------------------------------------------------
plot(gam_nle)

## -----------------------------------------------------------------------------
gam_tvnle <- rem(~ tvnl(diff.author.publication.activity) 
                   # automatically uses the transformed one when computing nl
                 
                   + paper.outdegree.popularity +
                   + author.sub.rep.1 +
                   + reference.sub.rep.1,
                   method="gam", time = "transformed_time",
                   data = dat_gam)
viz <- getViz(gam_tvnle$fit)

## -----------------------------------------------------------------------------
summary(gam_tvnle)

## ----gam_tvnle----------------------------------------------------------------
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

## -----------------------------------------------------------------------------
AIC(gam_tve)
AIC(gam_tve_transformed)
AIC(gam_nle)
AIC(gam_tvnle)

