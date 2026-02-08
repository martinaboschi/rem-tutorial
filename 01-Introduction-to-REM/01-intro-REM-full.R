## ----setup, include=FALSE--------------------------------------------------------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE------------------------------------------------------------------------------------------
# install.packages("mgcv")
# install.packages("RColorBrewer")
library(mgcv)
library(RColorBrewer)

## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
pal.blue <- brewer.pal(9, "Blues")
pal.rose <- brewer.pal(9, "RdPu")
colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#999999", "#66C2A5", "#FC8D62")

## --------------------------------------------------------------------------------------------------------------------------
load(file="01-Data/01-Inputs/input00.RData")

## --------------------------------------------------------------------------------------------------------------------------
head(FR[,c("LifeForm", "Taxon",
           "Region", "FirstRecord",
           "Source")])

## --------------------------------------------------------------------------------------------------------------------------
native[sample(1:nrow(native), 10), c("species", "region")]

## --------------------------------------------------------------------------------------------------------------------------
head(first_records[,c("year","lf",
           "species", 
           "region")])

## --------------------------------------------------------------------------------------------------------------------------
invaded.regions <- function(sp.n, r.n, y, native, first_records){
  
  # Convert input arguments to numeric type if not already
  sp.n <- as.numeric(sp.n)
  r.n <- as.numeric(r.n)
  y <- as.numeric(y)
  
  # Get unique combinations of species number and region number from native data
  t <- unique(as.vector(subset(native, sp.num == sp.n, r.num)))
  
  # Get region numbers from first records data where species number matches
  pr <- as.vector(subset(first_records, sp.num == sp.n, r.num))
  
  # If the invasion is both present in first records data and in native range
  # consider the former as actual piece of information
  t <- setdiff(t, pr)
  
  # Find indices of first records occurring before end date for the species
  set.sp <- which(first_records$sp.num == sp.n & first_records$year < y)
  
  # Combine regions in native range with regions in first records before date
  t <- na.omit(c(t, first_records$r.num[set.sp]))
  
  # Do not consider the involved region
  inv <- unlist(setdiff(t, r.n))
  
  # Return invaded regions
  return(inv)
}

## --------------------------------------------------------------------------------------------------------------------------
# Define sender set - species
spec <- unique(first_records$species)
(s <- length(spec))
# Define receiver set - regions
reg.lf <- unique(first_records$region)
(r <- length(reg.lf))

## --------------------------------------------------------------------------------------------------------------------------
reg <- colnames(data_distance)

## --------------------------------------------------------------------------------------------------------------------------
creating_case_control_dataset <- function(first_records,
                                 spec, reg.lf, reg,
                                 seed=1234
                                 ){
  
  reg.lf.num <- match(reg.lf, reg)
  s <- length(spec)
  r <- length(reg.lf)
  
  ## POSSIBLE (S,R) INTERACTIONS ####
  # Initialize a vector to keep track of the risk set size over time
  at.risk <- NULL
  # Create a matrix to indicate possible (species, region) interactions
  alien.occ <- matrix(0, nrow = s, ncol = r) 
  rownames(alien.occ) <- spec
  colnames(alien.occ) <- reg.lf
  for (n.sp in 1:s){
    # Identify native regions for each species
    nat.id <- unique(native$r.num[native$sp.num == n.sp]) 
    # Identify regions where the species is not native
    possible.to <- setdiff(reg.lf.num,nat.id)
    # Mark these regions as possible invasion sites for the species
    alien.occ[n.sp,reg[possible.to]] <- 1
  }
  
  ## COLLECTING INFORMATION ####
  dat.gam <- data.frame(matrix(NA, nrow=nrow(first_records),ncol=6))
  colnames(dat.gam) <- c("y", "year", 
                         "sp1", "r1", 
                         "sp2", "r2")
  # The response is fixed and equal to 1
  dat.gam[,1] <- rep(1, nrow(first_records))
  
  set.seed(seed)
  # For each FR:
  for (i in 1:nrow(first_records)){
    
    ### INFORMATION CONCERNING THE EVENT ####
    # year of the invasion event
    dat.gam[i,2] <- year <- first_records[i,"year"]
    # invading species
    dat.gam[i,3] <- s.ev <- first_records[i,"species"]
    # invaded country
    dat.gam[i,4] <- r.ev <- first_records[i,"region"]
    
    ### POSSIBLE EVENTS ####
    # Events occurred at the same time of the considered event
    # are removed from the risk set
    sub_stp <- first_records[first_records$year==year,
                             c("species", "region")]
    ni <- nrow(sub_stp)
    for (j in 1:ni){
      # Mark these (species, region) pairs as not at risk
      alien.occ[sub_stp[j,1],sub_stp[j,2]] <- 0
    }
    at.risk <- c(at.risk, sum(alien.occ==1))
    
    ### SAMPLING THE NON-EVENT ####
    sr.nv<-sample(which(alien.occ!=0),1)
    # species non-event
    dat.gam[i,5] <- s.nv <- spec[(sr.nv-1)%%s+1]
    # region non-event
    dat.gam[i,6] <- r.nv <- reg.lf[(sr.nv-1)%/%s+1]
  }
  
  return(dat.gam)
}

## --------------------------------------------------------------------------------------------------------------------------
dat.gam <- creating_case_control_dataset(first_records, 
                                         spec, reg.lf, reg)

## --------------------------------------------------------------------------------------------------------------------------
dat.gam$sp1.num <- match(dat.gam$sp1, spec)
dat.gam$r1.num <- match(dat.gam$r1, reg)
dat.gam$sp2.num <- match(dat.gam$sp2, spec)
dat.gam$r2.num <- match(dat.gam$r2, reg)

## --------------------------------------------------------------------------------------------------------------------------
head(dat.gam[c("year","sp1", "r1","sp2", "r2")])
head(dat.gam[c("year","sp1.num", "r1.num","sp2.num", "r2.num")])

## --------------------------------------------------------------------------------------------------------------------------
#save.image("01-Data/01-Inputs/input01.RData")

## --------------------------------------------------------------------------------------------------------------------------
load(file="01-Data/01-Inputs/input01.RData")

## --------------------------------------------------------------------------------------------------------------------------
climatic_dissimilarity <- function(sp.n, r.n, y, native, first_records, 
                                   reg, data_temperature){
  
  # Convert input arguments to numeric type if not already
  sp.n <- as.numeric(sp.n)
  r.n <- as.numeric(r.n)
  y <- as.numeric(y)
  
  # Find regions invaded by the species before the current time
  inv <- invaded.regions(sp.n = sp.n, 
                         r.n = r.n, 
                         y = y, 
                         native = native, 
                         first_records = first_records)
  
  # Consider the minimum absolute difference in temperature
  # between the region of interest and those already invaded
  avg_temp_invaded <- data_temperature[data_temperature[,1] %in% reg[inv], 2]
  # Get the average temperature value for the region of interest
  avg_temp_interest <- data_temperature[r.n, 2]
  # Find the minimum absolute difference in temperature
  dt.value <- min(abs(avg_temp_invaded - avg_temp_interest))
  
  # Return the calculated climatic dissimilarity value
  return(dt.value)
}

## --------------------------------------------------------------------------------------------------------------------------
dat.gam$dt1 <- apply(dat.gam[,c("sp1.num", "r1.num", "year")], 1, 
                            function(x) climatic_dissimilarity(sp.n = x[1], 
                                        r.n = x[2], 
                                        y = x[3], 
                                        native = native, 
                                        first_records = first_records, 
                                        reg = reg, 
                                        data_temperature = data_temperature))
dat.gam$dt2 <- apply(dat.gam[,c("sp2.num", "r2.num", "year")], 1, 
                            function(x) climatic_dissimilarity(sp.n = x[1], 
                                        r.n = x[2], 
                                        y = x[3], 
                                        native = native, 
                                        first_records = first_records,
                                        reg = reg, 
                                        data_temperature = data_temperature))
dat.gam$dt = dat.gam$dt1 - dat.gam$dt2

## --------------------------------------------------------------------------------------------------------------------------
gam_dt.only <- gam(y ~ dt - 1,
    family="binomial"(link = 'logit'), data=dat.gam)

## --------------------------------------------------------------------------------------------------------------------------
summary(gam_dt.only)

## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
#save(gam_dt.only, file="01-Data/02-Gam-Fits/gam_dt.only.RData")


## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
knitr::include_graphics("03-Images/trade-covariate-computation.pdf")

## --------------------------------------------------------------------------------------------------------------------------
t <- which(data_trade$transfer < 0)
data_trade$transfer[t] <- 0

## --------------------------------------------------------------------------------------------------------------------------
trade.funct <- function(inv, r.n, y, reg, data_trade){
  
  # Convert input arguments to numeric type if not already
  r.n <- as.numeric(r.n)
  y <- as.numeric(y)
  
  # Check if there are already invaded countries
  if(length(inv)!=0){
    
    # Find rows that involve invaded countries as sending trade
    u <- which(data_trade$FromRegion %in% reg[inv])
    # Find rows that involve region of interest as receiving trade
    v <- which(data_trade$ToRegion == reg[r.n])
    # Find the intersection of the two sets
    w <- intersect(u,v)
    x <- data_trade[w,]
    # Consider the trade instances occurred before or at the time of interest
    x <- x[x$year<=y,]
    trade_value <- NULL
    # If there are rows in the filtered dataset
    if(nrow(x)>0){
      # For each invaded country, the maximum year is recorded
      o <- aggregate(x$year, list(x$FromRegion), FUN=max)
      # For each of them, the corresponding transfer is stored
      for (o.i in 1:nrow(o)){
        trade_value <- c(trade_value, 
                         x$transfer[x$FromRegion==o[o.i,1] & 
                                      x$year==o[o.i,2]])}
    }
  } else { 
  # If there are not already invaded countries, trade is set equal to 0
    trade_value <- 0
  }
  # Compute the log-transformed sum of trade values (with an added constant 1)
  log_trade.value <- ifelse(length(trade_value)>0, 
                            log(sum(trade_value, na.rm =T)+1),0)
  
  # Return the computed log-transformed trade value
  return(log_trade.value)
}

## --------------------------------------------------------------------------------------------------------------------------
log_trade <- function(sp.n, r.n, y, native, first_records, reg, data_trade){
  
  inv <- invaded.regions(sp.n = sp.n, 
                         r.n = r.n, 
                         y = y, 
                         native = native, 
                         first_records = first_records)
  
  log_trade.value <- ifelse(r.n==match("USACanada", reg), 
                            mean(trade.funct(inv = inv,
                                             r.n = match("United States",reg),
                                             y = y, 
                                             reg = reg,
                                             data_trade = data_trade),
                                 trade.funct(inv = inv,
                                             r.n = match("Canada",reg),
                                             y = y, 
                                             reg = reg,
                                             data_trade = data_trade)),
                            trade.funct(inv = inv,
                                        r.n = r.n,
                                        y = y,
                                        reg = reg, 
                                        data_trade = data_trade))
  
  return(log_trade.value)
}

## --------------------------------------------------------------------------------------------------------------------------
dat.gam$tr1 <- apply(dat.gam[,c("sp1.num", "r1.num", "year")], 1, 
                            function(x) log_trade(sp.n = x[1], 
                                                  r.n = x[2], 
                                                  y = x[3], 
                                                  native = native, 
                                                  first_records = 
                                                    first_records,
                                                  reg = reg,
                                                  data_trade = data_trade))
dat.gam$tr2 <- apply(dat.gam[,c("sp2.num", "r2.num", "year")], 1, 
                            function(x) log_trade(x[1], x[2], x[3], 
                                                  native = native, 
                                                     first_records = 
                                                       first_records,
                                                  reg = reg,
                                                  data_trade = data_trade))
dat.gam$tr = dat.gam$tr1 - dat.gam$tr2

## --------------------------------------------------------------------------------------------------------------------------
# nx1 vector of trade covariate for events
x.ev <- dat.gam$tr1
# nx1 vector of trade covariate for non-events
x.nv <- dat.gam$tr2
x <- x.ev - x.nv # difference between ev and nv

stp <- dat.gam$year ## time 

gam_tr.only <- gam(y ~ s(stp, by=x) - 1,
    family="binomial"(link = 'logit'), data=dat.gam)

plot(gam_tr.only)

## --------------------------------------------------------------------------------------------------------------------------
filter_data <- which(dat.gam$tr != 0)
lp_matrix <- predict.gam(gam_tr.only, type="lpmatrix")[filter_data,]
basis_fcn_tr <- lp_matrix/dat.gam$tr[filter_data]


plot(dat.gam$year[filter_data],basis_fcn_tr[,1], type="l", lwd=1.5, ylim=c(-3.5, 4.5),
     xlab="Time",
     ylab="b(Time)", col = colors[1])
for (l in 2:10) {
  lines(y = basis_fcn_tr[,l],
         x = dat.gam$year[filter_data],
         lwd = 0.8,
         col = colors[l])
}
legend("topright",
       legend=c(sapply(1:10, function(x) paste("Basis",
                                               as.character(x)))),
       col=colors,
       lwd=rep(0.8, 10),
       cex=0.45)

## --------------------------------------------------------------------------------------------------------------------------
coefficients(gam_tr.only)

## --------------------------------------------------------------------------------------------------------------------------
predicted_effect_tr <- as.vector(coefficients(gam_tr.only) %*%t(lp_matrix/dat.gam$tr[filter_data]))
data_effect_tr <- data.frame(x = stp[filter_data],
                             y = predicted_effect_tr)

plot(data_effect_tr, type="l", lwd=1.5, ylim=c(-0.5, 1.2),
     xlab="Time",
     ylab="Contribution to the log-hazard")
for (l in 1:10) {
  lines(y = coefficients(gam_tr.only)[l] *
          lp_matrix[,l]/dat.gam$tr[filter_data],
         x = dat.gam$year[filter_data],
         lwd = 0.8,
         col = colors[l])
}
legend("topright",
       legend=c("Time-varying effect",
                sapply(1:10, function(x) paste("Contr. basis",
                                               as.character(x)))),
       col=c(1, colors),
       lwd=c(1.5,rep(0.8, 10)),
       cex=0.45)

## --------------------------------------------------------------------------------------------------------------------------
#save(gam_tr.only, file="01-Data/02-Gam-fits/gam_tr.only.RData")

## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
knitr::include_graphics("03-Images/distance-covariate-computation.pdf")

## --------------------------------------------------------------------------------------------------------------------------
log_distance <- function(sp.n, r.n, y, native, first_records, data_distance){
  
  # Convert input arguments to numeric type if not already
  sp.n <- as.numeric(sp.n)
  r.n <- as.numeric(r.n)
  y <- as.numeric(y)
  
  # Find regions invaded by the species before the current time
  inv <- invaded.regions(sp.n = sp.n, 
                         r.n = r.n, 
                         y = y, 
                         native = native, 
                         first_records = first_records)
  
  # Consider the logarithm of the minimum distance
  # between the region of interest and those already invaded
  log_dist.value <- log(min(data_distance[r.n, inv])+1)
  
  # Return the calculated distance
  return(log_dist.value)
}

## --------------------------------------------------------------------------------------------------------------------------
dat.gam$d1 <- apply(dat.gam[,c("sp1.num", "r1.num", "year")], 1, 
                            function(x) log_distance(x[1], x[2], x[3], 
                                                     native = native, 
                                                     first_records = 
                                                       first_records,
                                                     data_distance = data_distance))
dat.gam$d2 <- apply(dat.gam[,c("sp2.num", "r2.num", "year")], 1, 
                            function(x) log_distance(x[1], x[2], x[3], 
                                                     native = native, 
                                                     first_records = 
                                                       first_records,
                                                     data_distance = data_distance))
dat.gam$d = dat.gam$d1 - dat.gam$d2

## ---- echo=FALSE-----------------------------------------------------------------------------------------------------------
knitr::include_graphics("03-Images/distance-covariate-computation.pdf")

## --------------------------------------------------------------------------------------------------------------------------
x.ev <- dat.gam$d1
x.nv <- dat.gam$d2
unit <- rep(1, nrow(dat.gam))

X = cbind(x.ev,x.nv)
I = cbind(unit,-unit)	

gam_d.only <- gam(y ~ s(X, by=I) - 1,
    family="binomial"(link = 'logit'), data=dat.gam)

plot(gam_d.only)


## --------------------------------------------------------------------------------------------------------------------------
lp_matrix <- predict.gam(gam_d.only, type="lpmatrix",
                         newdata = data.frame(X=dat.gam$d1,
                                              I=1))
predicted_effect_d <- as.vector(coefficients(gam_d.only) %*% t(lp_matrix))
data_effect_d <- data.frame(x = sort(exp(dat.gam$d1)-1),
                            y = predicted_effect_d[order(exp(dat.gam$d1)-1)])

plot(data_effect_d, lwd=1.5,
     xlab="Distance",
     ylab="Contribution to the log-hazard",
     ylim=c(-2, 3), type = "l",
     col=0)
for (l in 1:9) {
  lines(y = coefficients(gam_d.only)[l]*lp_matrix[,l][order(exp(dat.gam$d1)-1)],
         x = sort(exp(dat.gam$d1)-1),
         lwd = 0.8,
         cex=0.2,
         col = colors[l])
}
lines(data_effect_d, cex=0.4,
       col=1)
legend("topright",
       legend=c("Non-linear effect",
                sapply(1:9, function(x) paste("Contr. basis",
                                               as.character(x)))),
       col=c(1, colors),
       lwd=c(1.5,rep(0.8, 9)),
       cex=0.45)

## --------------------------------------------------------------------------------------------------------------------------
#save(gam_d.only, file="01-Data/02-Gam-Fits/gam_d.only.RData")

## --------------------------------------------------------------------------------------------------------------------------
#save.image("01-Data/01-Inputs/input02.RData")

## --------------------------------------------------------------------------------------------------------------------------
load("01-Data/01-Inputs/input02.RData")

## --------------------------------------------------------------------------------------------------------------------------
sp1 <- dat.gam$sp1
sp2 <- dat.gam$sp2
sp <- factor(c(sp1,sp2))
dim(sp) <- c(length(sp1),2)

## --------------------------------------------------------------------------------------------------------------------------
unit <- rep(1, nrow(dat.gam))
I = cbind(unit,-unit)	
gam_sp.only <- gam(y ~ s(sp, by=I, bs="re") - 1,
    family="binomial"(link = 'logit'), data=dat.gam)

re.species <- coefficients(gam_sp.only)
names(re.species) <- levels(sp)

# 5 most invasive species
sort(re.species, decreasing = TRUE)[1:5]
# 5 least invasive species
sort(re.species)[1:5]

## --------------------------------------------------------------------------------------------------------------------------
#save(gam_sp.only, file="01-Data/02-Gam-Fits/gam_sp.only.RData")

## --------------------------------------------------------------------------------------------------------------------------
load("01-Data/01-Inputs/input02.RData")

## --------------------------------------------------------------------------------------------------------------------------
unit <- rep(1, nrow(dat.gam))
stp = dat.gam$year
X = cbind(dat.gam$d1,dat.gam$d2)
I = cbind(unit,-unit)		

## --------------------------------------------------------------------------------------------------------------------------
gam_complete <- gam(y ~ dt + 
                      s(stp, by=tr) +
                      s(X, by=I) +
                      s(sp, by=I, bs="re") - 1,
    family="binomial"(link = 'logit'), data=dat.gam)

## --------------------------------------------------------------------------------------------------------------------------
load(file="01-Data/02-Gam-Fits/gam_dt.only.RData")
load(file="01-Data/02-Gam-Fits/gam_tr.only.RData")
load(file="01-Data/02-Gam-Fits/gam_d.only.RData")
load(file="01-Data/02-Gam-Fits/gam_sp.only.RData")

## --------------------------------------------------------------------------------------------------------------------------
AIC(gam_dt.only)
AIC(gam_tr.only)
AIC(gam_d.only)
AIC(gam_sp.only)
AIC(gam_complete)

