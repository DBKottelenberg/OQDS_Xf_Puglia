###########
##       ##
## Setup ##
##       ##
###########

rm(list = ls())
setwd("C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R")

library(readxl)
library(reshape2)
library(geosphere)
library(ggplot2)
library(emdbook)
library(tidyverse)
library(dplyr)
library(xlsx)

# Set ggplot graph theme
theme_set(theme_bw())

dat_base <- read_xlsx("xylelladata.xlsx")

###################################################
###################################################
######                                       ######
######           LOGISTIC FUNCTION           ######
######                                       ######
###################################################
###################################################

#############
### YEARS ###
#############

###############################
##                           ##
## Data reading and morphing ##
##                           ##
###############################

# Read in data
dat <- dat_base

dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter((lon != 0) | (lat != 0))

# Set number of years, used for length of loops and such
years <- unique(dat$year)

gallipoli <- c(lon=17.992615,lat=40.055851)


################
##            ##
## Simulation ##
##            ##
################

cvector <- seq(from = 5, to = 16, by = 1)
x50vector <- seq(from = -40, to = -5, by = 5)

# Set the number of times every parameter set is simulated
n_sim <- 10

# Prepare data frame to save estimated means and sd's 
dat_accuracy <- tibble("c" = sort(rep(cvector, times = length(x50vector))), 
                       "x50" = rep(x50vector, times = length(cvector)), 
                       "accuracy_log_od_mean" = NA,
                       "accuracy_log_od_sd" = NA,
                       "accuracy_log_pc_mean" = NA,
                       "accuracy_log_pc_sd" = NA,
                       "accuracy_log_nc_mean" = NA,
                       "accuracy_log_nc_sd" = NA,
                       "accuracy_log_ncpc_mean" = NA,
                       "accuracy_log_ncpc_sd" = NA,
                       "accuracy_conexp_od_mean" = NA,
                       "accuracy_conexp_od_sd" = NA,
                       "accuracy_conexp_pc_mean" = NA,
                       "accuracy_conexp_pc_sd" = NA,
                       "accuracy_conexp_nc_mean" = NA,
                       "accuracy_conexp_nc_sd" = NA,
                       "accuracy_conexp_ncpc_mean" = NA,
                       "accuracy_conexp_ncpc_sd" = NA,
                       "accuracy_hill_od_mean" = NA,
                       "accuracy_hill_od_sd" = NA,
                       "accuracy_hill_pc_mean" = NA,
                       "accuracy_hill_pc_sd" = NA,
                       "accuracy_hill_nc_mean" = NA,
                       "accuracy_hill_nc_sd" = NA,
                       "accuracy_hill_ncpc_mean" = NA,
                       "accuracy_hill_ncpc_sd" = NA)

for(k in 1:nrow(dat_accuracy)){ # For every combination of parameters
  sim_par <- c("a" = 0.08, # Set the parameter set
               "x50" = dat_accuracy$x50[k], 
               "c" = dat_accuracy$c[k], 
               "theta" = 1) 
  
  # Prepare a list with the estimates for every simulation of this parameter set
  accuracy <- list("log_od" = numeric(n_sim),
                   "log_pc" = numeric(n_sim),
                   "log_nc" = numeric(n_sim),
                   "log_ncpc" = numeric(n_sim),
                   "conexp_od" = numeric(n_sim),
                   "conexp_pc" = numeric(n_sim),
                   "conexp_nc" = numeric(n_sim),
                   "conexp_ncpc" = numeric(n_sim),
                   "hill_od" = numeric(n_sim),
                   "hill_pc" = numeric(n_sim),
                   "hill_nc" = numeric(n_sim),
                   "hill_ncpc" = numeric(n_sim))
  
  for(l in 1:n_sim){ # For every iteration of this parameter set
    
    #####  Generating Sampling Data #####
    
    # Generate sampling dataset
    ## Setup new sampled dataset
    dat_sam <- dat %>% 
      dplyr::select(lon, lat, year) %>% 
      mutate("dist" = NA)
    
    # Calculate distance of sampled points to Gallipoli
    dat_sam <- dat_sam %>% 
      mutate(dist = round(distHaversine(tibble(lon = lon, lat = lat), 
                                        c(gallipoli[1], gallipoli[2]))/1000))
    
    dat_sam2 <- dat_sam %>% 
      dplyr::select(year, dist) %>% 
      mutate(n = 1)
    
    dat_newsam <- as_tibble(melt(tapply(dat_sam2$n, list("year" = dat_sam2$year, # Find the number of samples in each distance circle, in each year
                                                         "dist" = dat_sam2$dist), sum))) %>% 
      mutate("pos" = NA)
    
    names(dat_newsam)[3] <- "n" # Name the number of samples column
    
    # Get sample results based on probability for this location and a binomial distribution
    for(i in 1:nrow(dat_newsam)){
      db1 = dat_newsam$year[i] >= 2014 # db1 is TRUE for 2014 and above, otherwise FALSE
      db2 = dat_newsam$year[i] >= 2015 # db2 is TRUE for 2015 and above, otherwise FALSE
      db3 = dat_newsam$year[i] >= 2016 
      db4 = dat_newsam$year[i] >= 2017
      db5 = dat_newsam$year[i] >= 2018
      
      # Deterministic model
      prob <- (1/(1 + exp(sim_par["a"] * (dat_newsam$dist[i] - (sim_par["x50"] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * sim_par["c"])))))
      
      # Get result from beta-binomial distribution
      dat_newsam$pos[i] <- rbetabinom(n = 1, prob = prob, size = dat_newsam$n[i], theta = sim_par["theta"])
    }
    
    dat_newsam <- dat_newsam %>% 
      mutate(prop = pos/n)
    
    
    ######  Analysis of sampling ######
    
    ###                           
    ### Rate of spread estimation 
    ###       
    
    #
    # Logistic function
    #
    
    # Setup logistic function: y = 1 / (1 + exp(a * (x - t*c))) where c is dependent on the year
    mlsp_fun_log <- function(par, x, z, n, db1, db2, db3, db4, db5){
      mu <- (1/(1 + exp(par[1] *                                                                 # par[1] is parameter 'a' (declining slope), par[2] is parameter 'x50', par[3] is parameter 'c',
                          (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))))   # 'dbx' can be TRUE or FALSE, the number of1's dependent on the 'dbx''s is the value of 't'.
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
      return(nll)
    }

    ##
    ## Original data rate of spread analysis
    ##
    
    dat_rs_an_od <- dat_newsam %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_od <- c(a = 0.05, x0 = -10, c = 10, theta = 1) # Starting parameters
    
    # Model fit
    mlsp_opt_log_od <- optim(par = mlsp_par_log_od, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_od$dist, 
                             z = dat_rs_an_od$pos, 
                             n = dat_rs_an_od$n, 
                             db1 = dat_rs_an_od$year >= 2014, # db1 is TRUE for 2014 and above
                             db2 = dat_rs_an_od$year >= 2015,  
                             db3 = dat_rs_an_od$year >= 2016, 
                             db4 = dat_rs_an_od$year >= 2017, 
                             db5 = dat_rs_an_od$year >= 2018, 
                             control = list(maxit = 10000), method = "Nelder-Mead") 
    mlsp_coefs_log_od <- mlsp_opt_log_od$par # Set parameters
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_pc <- dat_newsam %>% # Set new data frame to create cumulative data
      replace_na(list(n = 0, pos = 0, prop = 0)) %>%   # Set NA values to 0, so positives can later be added from earlier years
      arrange(year, dist)
      
    for(i in 1:nrow(dat_rs_an_pc)){
      yr = dat_rs_an_pc$year[i] # The year of this datapoint
      dst = dat_rs_an_pc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr + 1)] &
                                 dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr + 1)] &
                                 dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_pc$n[i] > 0, ndup/dat_rs_an_pc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous year
      pos_prev <- dat_rs_an_pc$pos[dat_rs_an_pc$year == (dat_rs_an_pc$year[i] - 1) &
                                        dat_rs_an_pc$dist == dat_rs_an_pc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous year (this value can be 0)
        dat_rs_an_pc$pos[i] = dat_rs_an_pc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous year in this distance circle
        dat_rs_an_pc$n[i] = dat_rs_an_pc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_pc <- dat_rs_an_pc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_pc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_pc <- optim(par = mlsp_par_log_pc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_pc$dist, 
                             z = dat_rs_an_pc$pos, 
                             n = dat_rs_an_pc$n, 
                             db1 = dat_rs_an_pc$year >= 2014, 
                             db2 = dat_rs_an_pc$year >= 2015,
                             db3 = dat_rs_an_pc$year >= 2016, 
                             db4 = dat_rs_an_pc$year >= 2017,
                             db5 = dat_rs_an_pc$year >= 2018,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_pc <- mlsp_opt_log_pc$par
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    dat_rs_an_nc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(year), dist)

    for(i in 1:nrow(dat_rs_an_nc)){
      yr = dat_rs_an_nc$year[i] # The year of this datapoint
      dst = dat_rs_an_nc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr - 1)] &
                             dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_nc$n[i] > 0, ndup/dat_rs_an_nc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_nc$n[i] - dat_rs_an_nc$pos[i]
      
      # Add the number of negatives of this year to the previous year, in the same distance circle, corrected for duplicate samples
      dat_rs_an_nc$n[dat_rs_an_nc$year == (yr - 1) & dat_rs_an_nc$dist == dst] <-  
        dat_rs_an_nc$n[dat_rs_an_nc$year == (yr - 1) & dat_rs_an_nc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_nc <- dat_rs_an_nc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_nc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_nc <- optim(par = mlsp_par_log_nc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_nc$dist, 
                             z = dat_rs_an_nc$pos, 
                             n = dat_rs_an_nc$n, 
                             db1 = dat_rs_an_nc$year >= 2014, 
                             db2 = dat_rs_an_nc$year >= 2015, 
                             db3 = dat_rs_an_nc$year >= 2016, 
                             db4 = dat_rs_an_nc$year >= 2017,
                             db5 = dat_rs_an_nc$year >= 2018,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_nc <- mlsp_opt_log_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_ncpc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(year), dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      yr = dat_rs_an_ncpc$year[i] # The year of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr - 1)] &
                             dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_ncpc$n[i] - dat_rs_an_ncpc$pos[i]
      
      # Add the number of negatives of this year to the previous year, in the same distance circle, corrected for duplicate samples
      dat_rs_an_ncpc$n[dat_rs_an_ncpc$year == (yr - 1) & dat_rs_an_ncpc$dist == dst] <-  
        dat_rs_an_ncpc$n[dat_rs_an_ncpc$year == (yr - 1) & dat_rs_an_ncpc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      arrange(year, dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      yr = dat_rs_an_ncpc$year[i] # The year of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr + 1)] &
                             dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr + 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous year
      pos_prev <- dat_rs_an_ncpc$pos[dat_rs_an_ncpc$year == (dat_rs_an_ncpc$year[i] - 1) &
                                        dat_rs_an_ncpc$dist == dat_rs_an_ncpc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous year (this value can be 0)
        dat_rs_an_ncpc$pos[i] = dat_rs_an_ncpc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous year in this distance circle
        dat_rs_an_ncpc$n[i] = dat_rs_an_ncpc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_ncpc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_ncpc <- optim(par = mlsp_par_log_ncpc, 
                               fn = mlsp_fun_log, 
                               x = dat_rs_an_ncpc$dist, 
                               z = dat_rs_an_ncpc$pos, 
                               n = dat_rs_an_ncpc$n, 
                               db1 = dat_rs_an_ncpc$year >= 2014, 
                               db2 = dat_rs_an_ncpc$year >= 2015, 
                               db3 = dat_rs_an_ncpc$year >= 2016, 
                               db4 = dat_rs_an_ncpc$year >= 2017,
                               db5 = dat_rs_an_ncpc$year >= 2018,
                               control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_ncpc <- mlsp_opt_log_ncpc$par 

    
    #                                           
    # Constrained Negative Exponential function 
    #     
    
    # Setup CNE function
    mlsp_fun_conexp <- function(par, x, z, n, db1, db2, db3, db4, db5){
      mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x) < 0.999,
                   (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x), 0.999) 
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_od <- c(b = 0.1, x100 = -5, c = 12, theta = 1) 
    
    # Model fit
    mlsp_opt_conexp_od <- optim(par = mlsp_par_conexp_od, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_od$dist, 
                                z = dat_rs_an_od$pos, 
                                n = dat_rs_an_od$n, 
                                db1 = dat_rs_an_od$year >= 2014, 
                                db2 = dat_rs_an_od$year >= 2015, 
                                db3 = dat_rs_an_od$year >= 2016, 
                                db4 = dat_rs_an_od$year >= 2017,
                                db5 = dat_rs_an_od$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_od <- mlsp_opt_conexp_od$par 
    
    
    ##
    ## Positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_pc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_pc <- optim(par = mlsp_par_conexp_pc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_pc$dist, 
                                z = dat_rs_an_pc$pos, 
                                n = dat_rs_an_pc$n, 
                                db1 = dat_rs_an_pc$year >= 2014, 
                                db2 = dat_rs_an_pc$year >= 2015, 
                                db3 = dat_rs_an_pc$year >= 2016, 
                                db4 = dat_rs_an_pc$year >= 2017,
                                db5 = dat_rs_an_pc$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_pc <- mlsp_opt_conexp_pc$par 
    
    
    ##
    ## Negative clairvoyance data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_nc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_nc <- optim(par = mlsp_par_conexp_nc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_nc$dist, 
                                z = dat_rs_an_nc$pos, 
                                n = dat_rs_an_nc$n, 
                                db1 = dat_rs_an_nc$year >= 2014,
                                db2 = dat_rs_an_nc$year >= 2015,
                                db3 = dat_rs_an_nc$year >= 2016, 
                                db4 = dat_rs_an_nc$year >= 2017,
                                db5 = dat_rs_an_nc$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_nc <- mlsp_opt_conexp_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_ncpc <- c(b = 0.1, x100 = -5, c = 14, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_ncpc <- optim(par = mlsp_par_conexp_ncpc, 
                                  fn = mlsp_fun_conexp, 
                                  x = dat_rs_an_ncpc$dist, 
                                  z = dat_rs_an_ncpc$pos, 
                                  n = dat_rs_an_ncpc$n, 
                                  db1 = dat_rs_an_ncpc$year >= 2014,
                                  db2 = dat_rs_an_ncpc$year >= 2015,
                                  db3 = dat_rs_an_ncpc$year >= 2016, 
                                  db4 = dat_rs_an_ncpc$year >= 2017,
                                  db5 = dat_rs_an_ncpc$year >= 2018,
                                  control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_ncpc <- mlsp_opt_conexp_ncpc$par 
    
    
    #               
    # Hill function 
    #               
    
    # Setup hill function
    mlsp_fun_hill <- function(par, x, z, n, db1, db2, db3, db4, db5){
      mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), par[1], 
                   (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2])/
                              (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2]))))
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_od <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1) 
    
    # Model fit
    mlsp_opt_hill_od <- optim(par = mlsp_par_hill_od, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_od$dist, 
                              z = dat_rs_an_od$pos, 
                              n = dat_rs_an_od$n, 
                              db1 = dat_rs_an_od$year >= 2014, 
                              db2 = dat_rs_an_od$year >= 2015,
                              db3 = dat_rs_an_od$year >= 2016, 
                              db4 = dat_rs_an_od$year >= 2017, 
                              db5 = dat_rs_an_od$year >= 2018, 
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_od <- mlsp_opt_hill_od$par
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_pc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_pc <- optim(par = mlsp_par_hill_pc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_pc$dist, 
                              z = dat_rs_an_pc$pos, 
                              n = dat_rs_an_pc$n, 
                              db1 = dat_rs_an_pc$year >= 2014, 
                              db2 = dat_rs_an_pc$year >= 2015, 
                              db3 = dat_rs_an_pc$year >= 2016, 
                              db4 = dat_rs_an_pc$year >= 2017,
                              db5 = dat_rs_an_pc$year >= 2018,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_pc <- mlsp_opt_hill_pc$par 
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_nc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_nc <- optim(par = mlsp_par_hill_nc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_nc$dist, 
                              z = dat_rs_an_nc$pos, 
                              n = dat_rs_an_nc$n, 
                              db1 = dat_rs_an_nc$year >= 2014,
                              db2 = dat_rs_an_nc$year >= 2015, 
                              db3 = dat_rs_an_nc$year >= 2016, 
                              db4 = dat_rs_an_nc$year >= 2017,
                              db5 = dat_rs_an_nc$year >= 2018,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_nc <- mlsp_opt_hill_nc$par 
    
    
    ##
    ## Positive carry-over & negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_ncpc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_ncpc <- optim(par = mlsp_par_hill_ncpc, 
                                fn = mlsp_fun_hill, 
                                x = dat_rs_an_ncpc$dist, 
                                z = dat_rs_an_ncpc$pos, 
                                n = dat_rs_an_ncpc$n, 
                                db1 = dat_rs_an_ncpc$year >= 2014, 
                                db2 = dat_rs_an_ncpc$year >= 2015, 
                                db3 = dat_rs_an_ncpc$year >= 2016, 
                                db4 = dat_rs_an_ncpc$year >= 2017,
                                db5 = dat_rs_an_ncpc$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_ncpc <- mlsp_opt_hill_ncpc$par 
    
    accuracy$log_od[l] = mlsp_coefs_log_od["c"] - sim_par["c"]  
    accuracy$log_pc[l] = mlsp_coefs_log_pc["c"] - sim_par["c"]  
    accuracy$log_nc[l] = mlsp_coefs_log_nc["c"] - sim_par["c"]  
    accuracy$log_ncpc[l] = mlsp_coefs_log_ncpc["c"] - sim_par["c"]  
    
    accuracy$conexp_od[l] = mlsp_coefs_conexp_od["c"] - sim_par["c"]  
    accuracy$conexp_pc[l] = mlsp_coefs_conexp_pc["c"] - sim_par["c"]  
    accuracy$conexp_nc[l] = mlsp_coefs_conexp_nc["c"] - sim_par["c"]  
    accuracy$conexp_ncpc[l] = mlsp_coefs_conexp_ncpc["c"] - sim_par["c"]  
    
    accuracy$hill_od[l] = mlsp_coefs_hill_od["c"] - sim_par["c"]  
    accuracy$hill_pc[l] = mlsp_coefs_hill_pc["c"] - sim_par["c"]  
    accuracy$hill_nc[l] = mlsp_coefs_hill_nc["c"] - sim_par["c"]  
    accuracy$hill_ncpc[l] = mlsp_coefs_hill_ncpc["c"] - sim_par["c"] 
    
    # Calculate progress percentage, to keep track of how far along the simulation is
    pro_perc <- ((k-1)*n_sim + l)*((100/(nrow(dat_accuracy)*n_sim)))
    print(pro_perc)
    
  }
  
  dat_accuracy$accuracy_log_od_mean[k] = mean(accuracy$log_od)
  dat_accuracy$accuracy_log_od_sd[k] = sd(accuracy$log_od)
  dat_accuracy$accuracy_log_pc_mean[k] = mean(accuracy$log_pc)
  dat_accuracy$accuracy_log_pc_sd[k] = sd(accuracy$log_pc)
  dat_accuracy$accuracy_log_nc_mean[k] = mean(accuracy$log_nc)
  dat_accuracy$accuracy_log_nc_sd[k] = sd(accuracy$log_nc)
  dat_accuracy$accuracy_log_ncpc_mean[k] = mean(accuracy$log_ncpc)
  dat_accuracy$accuracy_log_ncpc_sd[k] = sd(accuracy$log_ncpc)
  dat_accuracy$accuracy_conexp_od_mean[k] = mean(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_od_sd[k] = sd(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_pc_mean[k] = mean(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_pc_sd[k] = sd(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_nc_mean[k] = mean(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_nc_sd[k] = sd(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_ncpc_mean[k] = mean(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_conexp_ncpc_sd[k] = sd(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_hill_od_mean[k] = mean(accuracy$hill_od)
  dat_accuracy$accuracy_hill_od_sd[k] = sd(accuracy$hill_od)
  dat_accuracy$accuracy_hill_pc_mean[k] = mean(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_pc_sd[k] = sd(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_nc_mean[k] = mean(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_nc_sd[k] = sd(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_ncpc_mean[k] = mean(accuracy$hill_ncpc)
  dat_accuracy$accuracy_hill_ncpc_sd[k] = sd(accuracy$hill_ncpc)
}

write.xlsx(dat_accuracy, "C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy.xlsx", sheetName = "Years")
dat_accuracy_years <- read.xlsx("C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy.xlsx", sheetName = "Years")


###############
### SEASONS ###
###############

###############################
##                           ##
## Data reading and morphing ##
##                           ##
###############################

dat <- dat_base

dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter((lon != 0) | (lat != 0))

dat <- dat %>% 
  mutate(season = ifelse((dat$year == 2013) | (dat$year == 2014 & dat$month %in% 1:3), 1, 
                         ifelse((dat$year == 2014 & dat$month %in% 4:12) | (dat$year == 2015 & dat$month %in% 1:3), 2, 
                                ifelse((dat$year == 2015 & dat$month %in% 4:12) | (dat$year == 2016 & dat$month %in% 1:3), 3, 
                                       ifelse((dat$year == 2016 & dat$month %in% 4:12) | (dat$year == 2017 & dat$month %in% 1:3), 4, 5)))))


# Set number of seasons, used for length of loops and such
seasons <- unique(dat$season)

gallipoli <- c(lon=17.992615,lat=40.055851)


################
##            ##
## Simulation ##
##            ##
################

cvector <- seq(from = 5, to = 16, by = 1)
x50vector <- seq(from = -40, to = -5, by = 5)

n_sim <- 10

dat_accuracy <- tibble("c" = sort(rep(cvector, times = length(x50vector))), 
                       "x50" = rep(x50vector, times = length(cvector)), 
                       "accuracy_log_od_mean" = NA,
                       "accuracy_log_od_sd" = NA,
                       "accuracy_log_pc_mean" = NA,
                       "accuracy_log_pc_sd" = NA,
                       "accuracy_log_nc_mean" = NA,
                       "accuracy_log_nc_sd" = NA,
                       "accuracy_log_ncpc_mean" = NA,
                       "accuracy_log_ncpc_sd" = NA,
                       "accuracy_conexp_od_mean" = NA,
                       "accuracy_conexp_od_sd" = NA,
                       "accuracy_conexp_pc_mean" = NA,
                       "accuracy_conexp_pc_sd" = NA,
                       "accuracy_conexp_nc_mean" = NA,
                       "accuracy_conexp_nc_sd" = NA,
                       "accuracy_conexp_ncpc_mean" = NA,
                       "accuracy_conexp_ncpc_sd" = NA,
                       "accuracy_hill_od_mean" = NA,
                       "accuracy_hill_od_sd" = NA,
                       "accuracy_hill_pc_mean" = NA,
                       "accuracy_hill_pc_sd" = NA,
                       "accuracy_hill_nc_mean" = NA,
                       "accuracy_hill_nc_sd" = NA,
                       "accuracy_hill_ncpc_mean" = NA,
                       "accuracy_hill_ncpc_sd" = NA)

for(k in 1:nrow(dat_accuracy)){
  sim_par <- c("a" = 0.08, "x50" = dat_accuracy$x50[k], "c" = dat_accuracy$c[k], "theta" = 1)
  
  accuracy <- list("log_od" = numeric(n_sim),
                   "log_pc" = numeric(n_sim),
                   "log_nc" = numeric(n_sim),
                   "log_ncpc" = numeric(n_sim),
                   "conexp_od" = numeric(n_sim),
                   "conexp_pc" = numeric(n_sim),
                   "conexp_nc" = numeric(n_sim),
                   "conexp_ncpc" = numeric(n_sim),
                   "hill_od" = numeric(n_sim),
                   "hill_pc" = numeric(n_sim),
                   "hill_nc" = numeric(n_sim),
                   "hill_ncpc" = numeric(n_sim))
  
  for(l in 1:n_sim){
    
    #####  Generating Sampling Data #####
    
    # Generate sampling dataset
    ## Setup new sampled dataset
    dat_sam <- dat %>% 
      dplyr::select(lon, lat, season) %>% 
      mutate("dist" = NA)
    
    # Calculate distance of sampled points to Gallipoli
    dat_sam <- dat_sam %>% 
      mutate(dist = round(distHaversine(tibble(lon = lon, lat = lat), 
                                        c(gallipoli[1], gallipoli[2]))/1000))
    
    dat_sam2 <- dat_sam %>% 
      dplyr::select(season, dist) %>% 
      mutate(n = 1)
    
    dat_newsam <- as_tibble(melt(tapply(dat_sam2$n, list("season" = dat_sam2$season, # Find the number of samples in each distance circle, in each year
                                                         "dist" = dat_sam2$dist), sum))) %>% 
      mutate("pos" = NA)
    
    names(dat_newsam)[3] <- "n" # Name the number of samples column
    
    # Get sample results based on probability for this location and a binomial distribution
    for(i in 1:nrow(dat_newsam)){
      db1 = dat_newsam$season[i] >= 2 # db1 is TRUE for 2 and above, otherwise FALSE
      db2 = dat_newsam$season[i] >= 3 # db2 is TRUE for 3 and above, otherwise FALSE
      db3 = dat_newsam$season[i] >= 4 
      db4 = dat_newsam$season[i] >= 5
      
      # Deterministic model
      prob <- (1/(1 + exp(sim_par["a"] * (dat_newsam$dist[i] - (sim_par["x50"] + (1*db1 + 1*db2 + 1*db3 + 1*db4) * sim_par["c"])))))
      
      # Get result from beta-binomial distribution
      dat_newsam$pos[i] <- rbetabinom(n = 1, prob = prob, size = dat_newsam$n[i], theta = sim_par["theta"])
    }
    
    dat_newsam <- dat_newsam %>% 
      mutate(prop = pos/n)

    ######  Analysis of sampling ######
    
    ###                           
    ### Rate of spread estimation 
    ###       
    
    #
    # Logistic function
    #
    
    # Setup logistic function: y = 1 / (1 + exp(a * (x - t*c))) where c is dependent on the year
    mlsp_fun_log <- function(par, x, z, n, db1, db2, db3, db4){
      mu <- (1/(1 + exp(par[1] *                                                                 # par[1] is parameter 'a' (declining slope), par[2] is parameter 'x50', par[3] is parameter 'c',
                          (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))))   # 'dbx' can be TRUE or FALSE, the number of1's dependent on the 'dbx''s is the value of 't'.
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
      return(nll)
    }
    
    ##
    ## Original data rate of spread analysis
    ##
    
    dat_rs_an_od <- dat_newsam %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_od <- c(a = 0.05, x0 = -10, c = 10, theta = 1) # Starting parameters
    
    # Model fit
    mlsp_opt_log_od <- optim(par = mlsp_par_log_od, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_od$dist, 
                             z = dat_rs_an_od$pos, 
                             n = dat_rs_an_od$n, 
                             db1 = dat_rs_an_od$season >= 2, # db1 is TRUE for 2014 and above
                             db2 = dat_rs_an_od$season >= 3,  
                             db3 = dat_rs_an_od$season >= 4, 
                             db4 = dat_rs_an_od$season >= 5, 
                             control = list(maxit = 10000), method = "Nelder-Mead") 
    mlsp_coefs_log_od <- mlsp_opt_log_od$par # Set parameters
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_pc <- dat_newsam %>% # Set new data frame to create cumulative data
      replace_na(list(n = 0, pos = 0, prop = 0)) %>%   # Set NA values to 0, so positives can later be added from earlier seasons
      arrange(season, dist)
    
    for(i in 1:nrow(dat_rs_an_pc)){
      sn = dat_rs_an_pc$season[i] # The season of this datapoint
      dst = dat_rs_an_pc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn + 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn + 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_pc$n[i] > 0, ndup/dat_rs_an_pc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous season
      pos_prev <- dat_rs_an_pc$pos[dat_rs_an_pc$season == (dat_rs_an_pc$season[i] - 1) &
                                     dat_rs_an_pc$dist == dat_rs_an_pc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous season (this value can be 0)
        dat_rs_an_pc$pos[i] = dat_rs_an_pc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous season in this distance circle
        dat_rs_an_pc$n[i] = dat_rs_an_pc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_pc <- dat_rs_an_pc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_pc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_pc <- optim(par = mlsp_par_log_pc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_pc$dist, 
                             z = dat_rs_an_pc$pos, 
                             n = dat_rs_an_pc$n, 
                             db1 = dat_rs_an_pc$season >= 2, 
                             db2 = dat_rs_an_pc$season >= 3,
                             db3 = dat_rs_an_pc$season >= 4, 
                             db4 = dat_rs_an_pc$season >= 5,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_pc <- mlsp_opt_log_pc$par
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    dat_rs_an_nc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(season), dist)
    
    for(i in 1:nrow(dat_rs_an_nc)){
      sn = dat_rs_an_nc$season[i] # The season of this datapoint
      dst = dat_rs_an_nc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn - 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_nc$n[i] > 0, ndup/dat_rs_an_nc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_nc$n[i] - dat_rs_an_nc$pos[i]
      
      # Add the number of negatives of this season to the previous season, in the same distance circle, corrected for duplicate samples
      dat_rs_an_nc$n[dat_rs_an_nc$season == (sn - 1) & dat_rs_an_nc$dist == dst] <-  
        dat_rs_an_nc$n[dat_rs_an_nc$season == (sn - 1) & dat_rs_an_nc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_nc <- dat_rs_an_nc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_nc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_nc <- optim(par = mlsp_par_log_nc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_nc$dist, 
                             z = dat_rs_an_nc$pos, 
                             n = dat_rs_an_nc$n, 
                             db1 = dat_rs_an_nc$season >= 2, 
                             db2 = dat_rs_an_nc$season >= 3, 
                             db3 = dat_rs_an_nc$season >= 4, 
                             db4 = dat_rs_an_nc$season >= 5,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_nc <- mlsp_opt_log_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_ncpc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(season), dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      sn = dat_rs_an_ncpc$season[i] # The season of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn - 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_ncpc$n[i] - dat_rs_an_ncpc$pos[i]
      
      # Add the number of negatives of this season to the previous season, in the same distance circle, corrected for duplicate samples
      dat_rs_an_ncpc$n[dat_rs_an_ncpc$season == (sn - 1) & dat_rs_an_ncpc$dist == dst] <-  
        dat_rs_an_ncpc$n[dat_rs_an_ncpc$season == (sn - 1) & dat_rs_an_ncpc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      arrange(season, dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      sn = dat_rs_an_ncpc$season[i] # The season of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn + 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn + 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous season
      pos_prev <- dat_rs_an_ncpc$pos[dat_rs_an_ncpc$season == (dat_rs_an_ncpc$season[i] - 1) &
                                       dat_rs_an_ncpc$dist == dat_rs_an_ncpc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous season (this value can be 0)
        dat_rs_an_ncpc$pos[i] = dat_rs_an_ncpc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous season in this distance circle
        dat_rs_an_ncpc$n[i] = dat_rs_an_ncpc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_ncpc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_ncpc <- optim(par = mlsp_par_log_ncpc, 
                               fn = mlsp_fun_log, 
                               x = dat_rs_an_ncpc$dist, 
                               z = dat_rs_an_ncpc$pos, 
                               n = dat_rs_an_ncpc$n, 
                               db1 = dat_rs_an_ncpc$season >= 2, 
                               db2 = dat_rs_an_ncpc$season >= 3, 
                               db3 = dat_rs_an_ncpc$season >= 4, 
                               db4 = dat_rs_an_ncpc$season >= 5,
                               control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_ncpc <- mlsp_opt_log_ncpc$par 
    
    
    #                                           
    # Constrained Negative Exponential function 
    #     
    
    # Setup CNE function
    mlsp_fun_conexp <- function(par, x, z, n, db1, db2, db3, db4){
      mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))) * exp(-par[1] * x) < 0.999,
                   (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))) * exp(-par[1] * x), 0.999) 
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_od <- c(b = 0.1, x100 = -5, c = 12, theta = 1) 
    
    # Model fit
    mlsp_opt_conexp_od <- optim(par = mlsp_par_conexp_od, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_od$dist, 
                                z = dat_rs_an_od$pos, 
                                n = dat_rs_an_od$n, 
                                db1 = dat_rs_an_od$season >= 2, 
                                db2 = dat_rs_an_od$season >= 3, 
                                db3 = dat_rs_an_od$season >= 4, 
                                db4 = dat_rs_an_od$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_od <- mlsp_opt_conexp_od$par 
    
    
    ##
    ## Positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_pc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_pc <- optim(par = mlsp_par_conexp_pc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_pc$dist, 
                                z = dat_rs_an_pc$pos, 
                                n = dat_rs_an_pc$n, 
                                db1 = dat_rs_an_pc$season >= 2, 
                                db2 = dat_rs_an_pc$season >= 3, 
                                db3 = dat_rs_an_pc$season >= 4, 
                                db4 = dat_rs_an_pc$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_pc <- mlsp_opt_conexp_pc$par 
    
    
    ##
    ## Negative clairvoyance data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_nc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_nc <- optim(par = mlsp_par_conexp_nc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_nc$dist, 
                                z = dat_rs_an_nc$pos, 
                                n = dat_rs_an_nc$n, 
                                db1 = dat_rs_an_nc$season >= 2,
                                db2 = dat_rs_an_nc$season >= 3,
                                db3 = dat_rs_an_nc$season >= 4, 
                                db4 = dat_rs_an_nc$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_nc <- mlsp_opt_conexp_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_ncpc <- c(b = 0.1, x100 = -5, c = 14, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_ncpc <- optim(par = mlsp_par_conexp_ncpc, 
                                  fn = mlsp_fun_conexp, 
                                  x = dat_rs_an_ncpc$dist, 
                                  z = dat_rs_an_ncpc$pos, 
                                  n = dat_rs_an_ncpc$n, 
                                  db1 = dat_rs_an_ncpc$season >= 2,
                                  db2 = dat_rs_an_ncpc$season >= 3,
                                  db3 = dat_rs_an_ncpc$season >= 4, 
                                  db4 = dat_rs_an_ncpc$season >= 5,
                                  control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_ncpc <- mlsp_opt_conexp_ncpc$par 
    
    
    #               
    # Hill function 
    #               
    
    # Setup hill function
    mlsp_fun_hill <- function(par, x, z, n, db1, db2, db3, db4){
      mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), par[1], 
                   (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2])/
                              (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2]))))
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_od <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1) 
    
    # Model fit
    mlsp_opt_hill_od <- optim(par = mlsp_par_hill_od, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_od$dist, 
                              z = dat_rs_an_od$pos, 
                              n = dat_rs_an_od$n, 
                              db1 = dat_rs_an_od$season >= 2, 
                              db2 = dat_rs_an_od$season >= 3,
                              db3 = dat_rs_an_od$season >= 4, 
                              db4 = dat_rs_an_od$season >= 5, 
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_od <- mlsp_opt_hill_od$par
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_pc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_pc <- optim(par = mlsp_par_hill_pc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_pc$dist, 
                              z = dat_rs_an_pc$pos, 
                              n = dat_rs_an_pc$n, 
                              db1 = dat_rs_an_pc$season >= 2, 
                              db2 = dat_rs_an_pc$season >= 3, 
                              db3 = dat_rs_an_pc$season >= 4, 
                              db4 = dat_rs_an_pc$season >= 5,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_pc <- mlsp_opt_hill_pc$par 
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_nc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_nc <- optim(par = mlsp_par_hill_nc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_nc$dist, 
                              z = dat_rs_an_nc$pos, 
                              n = dat_rs_an_nc$n, 
                              db1 = dat_rs_an_nc$season >= 2,
                              db2 = dat_rs_an_nc$season >= 3, 
                              db3 = dat_rs_an_nc$season >= 4, 
                              db4 = dat_rs_an_nc$season >= 5,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_nc <- mlsp_opt_hill_nc$par 
    
    
    ##
    ## Positive carry-over & negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_ncpc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_ncpc <- optim(par = mlsp_par_hill_ncpc, 
                                fn = mlsp_fun_hill, 
                                x = dat_rs_an_ncpc$dist, 
                                z = dat_rs_an_ncpc$pos, 
                                n = dat_rs_an_ncpc$n, 
                                db1 = dat_rs_an_ncpc$season >= 2, 
                                db2 = dat_rs_an_ncpc$season >= 3, 
                                db3 = dat_rs_an_ncpc$season >= 4, 
                                db4 = dat_rs_an_ncpc$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_ncpc <- mlsp_opt_hill_ncpc$par 
    
    accuracy$log_od[l] = mlsp_coefs_log_od["c"] - sim_par["c"]  
    accuracy$log_pc[l] = mlsp_coefs_log_pc["c"] - sim_par["c"]  
    accuracy$log_nc[l] = mlsp_coefs_log_nc["c"] - sim_par["c"]  
    accuracy$log_ncpc[l] = mlsp_coefs_log_ncpc["c"] - sim_par["c"]  
    
    accuracy$conexp_od[l] = mlsp_coefs_conexp_od["c"] - sim_par["c"]  
    accuracy$conexp_pc[l] = mlsp_coefs_conexp_pc["c"] - sim_par["c"]  
    accuracy$conexp_nc[l] = mlsp_coefs_conexp_nc["c"] - sim_par["c"]  
    accuracy$conexp_ncpc[l] = mlsp_coefs_conexp_ncpc["c"] - sim_par["c"]  
    
    accuracy$hill_od[l] = mlsp_coefs_hill_od["c"] - sim_par["c"]  
    accuracy$hill_pc[l] = mlsp_coefs_hill_pc["c"] - sim_par["c"]  
    accuracy$hill_nc[l] = mlsp_coefs_hill_nc["c"] - sim_par["c"]  
    accuracy$hill_ncpc[l] = mlsp_coefs_hill_ncpc["c"] - sim_par["c"] 
    
    pro_perc <- ((k-1)*n_sim + l)*((100/(nrow(dat_accuracy)*n_sim)))
    print(pro_perc)
    
  }
  
  dat_accuracy$accuracy_log_od_mean[k] = mean(accuracy$log_od)
  dat_accuracy$accuracy_log_od_sd[k] = sd(accuracy$log_od)
  dat_accuracy$accuracy_log_pc_mean[k] = mean(accuracy$log_pc)
  dat_accuracy$accuracy_log_pc_sd[k] = sd(accuracy$log_pc)
  dat_accuracy$accuracy_log_nc_mean[k] = mean(accuracy$log_nc)
  dat_accuracy$accuracy_log_nc_sd[k] = sd(accuracy$log_nc)
  dat_accuracy$accuracy_log_ncpc_mean[k] = mean(accuracy$log_ncpc)
  dat_accuracy$accuracy_log_ncpc_sd[k] = sd(accuracy$log_ncpc)
  dat_accuracy$accuracy_conexp_od_mean[k] = mean(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_od_sd[k] = sd(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_pc_mean[k] = mean(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_pc_sd[k] = sd(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_nc_mean[k] = mean(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_nc_sd[k] = sd(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_ncpc_mean[k] = mean(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_conexp_ncpc_sd[k] = sd(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_hill_od_mean[k] = mean(accuracy$hill_od)
  dat_accuracy$accuracy_hill_od_sd[k] = sd(accuracy$hill_od)
  dat_accuracy$accuracy_hill_pc_mean[k] = mean(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_pc_sd[k] = sd(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_nc_mean[k] = mean(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_nc_sd[k] = sd(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_ncpc_mean[k] = mean(accuracy$hill_ncpc)
  dat_accuracy$accuracy_hill_ncpc_sd[k] = sd(accuracy$hill_ncpc)
}

write.xlsx(dat_accuracy, "C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy.xlsx", sheetName = "Seasons", append = TRUE)
dat_accuracy_seasons <- read.xlsx("C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy.xlsx", sheetName = "Seasons")


##############################################
##############################################
######                                  ######
######           CNE FUNCTION           ######
######                                  ######
##############################################
##############################################

#############
### YEARS ###
#############

###############################
##                           ##
## Data reading and morphing ##
##                           ##
###############################

# Read in data
dat <- dat_base

dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter((lon != 0) | (lat != 0))

# Set number of years, used for length of loops and such
years <- unique(dat$year)

gallipoli <- c(lon=17.992615,lat=40.055851)


################
##            ##
## Simulation ##
##            ##
################

# The length of the vectors and the value of n_sim determine how many times the simulation is run. 15 simulations take about 1 hour!
cvector <- seq(from = 10, to = 10, by = 1)
x100vector <- seq(from = -20, to = -20, by = 5)

n_sim <- 2

dat_accuracy <- tibble("c" = sort(rep(cvector, times = length(x100vector))), 
                       "x100" = rep(x100vector, times = length(cvector)), 
                       "accuracy_log_od_mean" = NA,
                       "accuracy_log_od_sd" = NA,
                       "accuracy_log_pc_mean" = NA,
                       "accuracy_log_pc_sd" = NA,
                       "accuracy_log_nc_mean" = NA,
                       "accuracy_log_nc_sd" = NA,
                       "accuracy_log_ncpc_mean" = NA,
                       "accuracy_log_ncpc_sd" = NA,
                       "accuracy_conexp_od_mean" = NA,
                       "accuracy_conexp_od_sd" = NA,
                       "accuracy_conexp_pc_mean" = NA,
                       "accuracy_conexp_pc_sd" = NA,
                       "accuracy_conexp_nc_mean" = NA,
                       "accuracy_conexp_nc_sd" = NA,
                       "accuracy_conexp_ncpc_mean" = NA,
                       "accuracy_conexp_ncpc_sd" = NA,
                       "accuracy_hill_od_mean" = NA,
                       "accuracy_hill_od_sd" = NA,
                       "accuracy_hill_pc_mean" = NA,
                       "accuracy_hill_pc_sd" = NA,
                       "accuracy_hill_nc_mean" = NA,
                       "accuracy_hill_nc_sd" = NA,
                       "accuracy_hill_ncpc_mean" = NA,
                       "accuracy_hill_ncpc_sd" = NA)

for(k in 1:nrow(dat_accuracy)){
  sim_par <- c("a" = 0.08, "x100" = dat_accuracy$x100[k], "c" = dat_accuracy$c[k], "theta" = 1)
  
  accuracy <- list("log_od" = numeric(n_sim),
                   "log_pc" = numeric(n_sim),
                   "log_nc" = numeric(n_sim),
                   "log_ncpc" = numeric(n_sim),
                   "conexp_od" = numeric(n_sim),
                   "conexp_pc" = numeric(n_sim),
                   "conexp_nc" = numeric(n_sim),
                   "conexp_ncpc" = numeric(n_sim),
                   "hill_od" = numeric(n_sim),
                   "hill_pc" = numeric(n_sim),
                   "hill_nc" = numeric(n_sim),
                   "hill_ncpc" = numeric(n_sim))
  
  for(l in 1:n_sim){
    ###################################################################################
    ############################                           ############################
    ############################  Generating Sampling Data ############################
    ############################                           ############################
    ###################################################################################
    
    #####  Generating Sampling Data #####
    
    # Generate sampling dataset
    ## Setup new sampled dataset
    dat_sam <- dat %>% 
      dplyr::select(lon, lat, year) %>% 
      mutate("dist" = NA)
    
    # Calculate distance of sampled points to Gallipoli
    dat_sam <- dat_sam %>% 
      mutate(dist = round(distHaversine(tibble(lon = lon, lat = lat), 
                                        c(gallipoli[1], gallipoli[2]))/1000))
    
    dat_sam2 <- dat_sam %>% 
      dplyr::select(year, dist) %>% 
      mutate(n = 1)
    
    dat_newsam <- as_tibble(melt(tapply(dat_sam2$n, list("year" = dat_sam2$year, # Find the number of samples in each distance circle, in each year
                                                         "dist" = dat_sam2$dist), sum))) %>% 
      mutate("pos" = NA)
    
    names(dat_newsam)[3] <- "n" # Name the number of samples column
    
    # Get sample results based on probability for this location and a binomial distribution
    for(i in 1:nrow(dat_newsam)){
      db1 = dat_newsam$year[i] >= 2014 # db1 is TRUE for 2014 and above, otherwise FALSE
      db2 = dat_newsam$year[i] >= 2015 # db2 is TRUE for 2015 and above, otherwise FALSE
      db3 = dat_newsam$year[i] >= 2016 
      db4 = dat_newsam$year[i] >= 2017
      db5 = dat_newsam$year[i] >= 2018
      
      # Deterministic model
      prob <- ifelse((1/(exp(-sim_par["a"]*(sim_par["x100"] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*sim_par["c"]))))*exp(-sim_par["a"]*dat_newsam$dist[i]) < 0.999,
                   (1/(exp(-sim_par["a"]*(sim_par["x100"] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*sim_par["c"]))))*exp(-sim_par["a"]*dat_newsam$dist[i]), 0.999)
      
      # Get result from beta-binomial distribution
      dat_newsam$pos[i] <- rbetabinom(n = 1, prob = prob, size = dat_newsam$n[i], theta = sim_par["theta"])
    }
    
    dat_newsam <- dat_newsam %>% 
      mutate(prop = pos/n)
      
      
    ######  Analysis of sampling ######
    
    ###                           
    ### Rate of spread estimation 
    ###       
    
    #
    # Logistic function
    #
    
    # Setup logistic function: y = 1 / (1 + exp(a * (x - t*c))) where c is dependent on the year
    mlsp_fun_log <- function(par, x, z, n, db1, db2, db3, db4, db5){
      mu <- (1/(1 + exp(par[1] *                                                                 # par[1] is parameter 'a' (declining slope), par[2] is parameter 'x50', par[3] is parameter 'c',
                          (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))))   # 'dbx' can be TRUE or FALSE, the number of1's dependent on the 'dbx''s is the value of 't'.
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
      return(nll)
    }
    
    ##
    ## Original data rate of spread analysis
    ##
    
    dat_rs_an_od <- dat_newsam %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_od <- c(a = 0.05, x0 = -10, c = 10, theta = 1) # Starting parameters
    
    # Model fit
    mlsp_opt_log_od <- optim(par = mlsp_par_log_od, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_od$dist, 
                             z = dat_rs_an_od$pos, 
                             n = dat_rs_an_od$n, 
                             db1 = dat_rs_an_od$year >= 2014, # db1 is TRUE for 2014 and above
                             db2 = dat_rs_an_od$year >= 2015,  
                             db3 = dat_rs_an_od$year >= 2016, 
                             db4 = dat_rs_an_od$year >= 2017, 
                             db5 = dat_rs_an_od$year >= 2018, 
                             control = list(maxit = 10000), method = "Nelder-Mead") 
    mlsp_coefs_log_od <- mlsp_opt_log_od$par # Set parameters
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_pc <- dat_newsam %>% # Set new data frame to create cumulative data
      replace_na(list(n = 0, pos = 0, prop = 0)) %>%   # Set NA values to 0, so positives can later be added from earlier years
      arrange(year, dist)
    
    for(i in 1:nrow(dat_rs_an_pc)){
      yr = dat_rs_an_pc$year[i] # The year of this datapoint
      dst = dat_rs_an_pc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr + 1)] &
                             dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr + 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_pc$n[i] > 0, ndup/dat_rs_an_pc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous year
      pos_prev <- dat_rs_an_pc$pos[dat_rs_an_pc$year == (dat_rs_an_pc$year[i] - 1) &
                                     dat_rs_an_pc$dist == dat_rs_an_pc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous year (this value can be 0)
        dat_rs_an_pc$pos[i] = dat_rs_an_pc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous year in this distance circle
        dat_rs_an_pc$n[i] = dat_rs_an_pc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_pc <- dat_rs_an_pc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_pc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_pc <- optim(par = mlsp_par_log_pc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_pc$dist, 
                             z = dat_rs_an_pc$pos, 
                             n = dat_rs_an_pc$n, 
                             db1 = dat_rs_an_pc$year >= 2014, 
                             db2 = dat_rs_an_pc$year >= 2015,
                             db3 = dat_rs_an_pc$year >= 2016, 
                             db4 = dat_rs_an_pc$year >= 2017,
                             db5 = dat_rs_an_pc$year >= 2018,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_pc <- mlsp_opt_log_pc$par
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    dat_rs_an_nc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(year), dist)
    
    for(i in 1:nrow(dat_rs_an_nc)){
      yr = dat_rs_an_nc$year[i] # The year of this datapoint
      dst = dat_rs_an_nc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr - 1)] &
                             dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_nc$n[i] > 0, ndup/dat_rs_an_nc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_nc$n[i] - dat_rs_an_nc$pos[i]
      
      # Add the number of negatives of this year to the previous year, in the same distance circle, corrected for duplicate samples
      dat_rs_an_nc$n[dat_rs_an_nc$year == (yr - 1) & dat_rs_an_nc$dist == dst] <-  
        dat_rs_an_nc$n[dat_rs_an_nc$year == (yr - 1) & dat_rs_an_nc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_nc <- dat_rs_an_nc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_nc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_nc <- optim(par = mlsp_par_log_nc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_nc$dist, 
                             z = dat_rs_an_nc$pos, 
                             n = dat_rs_an_nc$n, 
                             db1 = dat_rs_an_nc$year >= 2014, 
                             db2 = dat_rs_an_nc$year >= 2015, 
                             db3 = dat_rs_an_nc$year >= 2016, 
                             db4 = dat_rs_an_nc$year >= 2017,
                             db5 = dat_rs_an_nc$year >= 2018,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_nc <- mlsp_opt_log_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_ncpc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(year), dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      yr = dat_rs_an_ncpc$year[i] # The year of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr - 1)] &
                             dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_ncpc$n[i] - dat_rs_an_ncpc$pos[i]
      
      # Add the number of negatives of this year to the previous year, in the same distance circle, corrected for duplicate samples
      dat_rs_an_ncpc$n[dat_rs_an_ncpc$year == (yr - 1) & dat_rs_an_ncpc$dist == dst] <-  
        dat_rs_an_ncpc$n[dat_rs_an_ncpc$year == (yr - 1) & dat_rs_an_ncpc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      arrange(year, dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      yr = dat_rs_an_ncpc$year[i] # The year of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the year before
      ndup <- length(which(dat_sam$lon[dat_sam$year == yr] == dat_sam$lon[dat_sam$year == (yr + 1)] &
                             dat_sam$lat[dat_sam$year == yr] == dat_sam$lat[dat_sam$year == (yr + 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous year
      pos_prev <- dat_rs_an_ncpc$pos[dat_rs_an_ncpc$year == (dat_rs_an_ncpc$year[i] - 1) &
                                       dat_rs_an_ncpc$dist == dat_rs_an_ncpc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous year (this value can be 0)
        dat_rs_an_ncpc$pos[i] = dat_rs_an_ncpc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous year in this distance circle
        dat_rs_an_ncpc$n[i] = dat_rs_an_ncpc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_ncpc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_ncpc <- optim(par = mlsp_par_log_ncpc, 
                               fn = mlsp_fun_log, 
                               x = dat_rs_an_ncpc$dist, 
                               z = dat_rs_an_ncpc$pos, 
                               n = dat_rs_an_ncpc$n, 
                               db1 = dat_rs_an_ncpc$year >= 2014, 
                               db2 = dat_rs_an_ncpc$year >= 2015, 
                               db3 = dat_rs_an_ncpc$year >= 2016, 
                               db4 = dat_rs_an_ncpc$year >= 2017,
                               db5 = dat_rs_an_ncpc$year >= 2018,
                               control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_ncpc <- mlsp_opt_log_ncpc$par 
    
    
    #                                           
    # Constrained Negative Exponential function 
    #     
    
    # Setup CNE function
    mlsp_fun_conexp <- function(par, x, z, n, db1, db2, db3, db4, db5){
      mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x) < 0.999,
                   (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x), 0.999) 
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_od <- c(b = 0.1, x100 = -5, c = 12, theta = 1) 
    
    # Model fit
    mlsp_opt_conexp_od <- optim(par = mlsp_par_conexp_od, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_od$dist, 
                                z = dat_rs_an_od$pos, 
                                n = dat_rs_an_od$n, 
                                db1 = dat_rs_an_od$year >= 2014, 
                                db2 = dat_rs_an_od$year >= 2015, 
                                db3 = dat_rs_an_od$year >= 2016, 
                                db4 = dat_rs_an_od$year >= 2017,
                                db5 = dat_rs_an_od$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_od <- mlsp_opt_conexp_od$par 
    
    
    ##
    ## Positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_pc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_pc <- optim(par = mlsp_par_conexp_pc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_pc$dist, 
                                z = dat_rs_an_pc$pos, 
                                n = dat_rs_an_pc$n, 
                                db1 = dat_rs_an_pc$year >= 2014, 
                                db2 = dat_rs_an_pc$year >= 2015, 
                                db3 = dat_rs_an_pc$year >= 2016, 
                                db4 = dat_rs_an_pc$year >= 2017,
                                db5 = dat_rs_an_pc$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_pc <- mlsp_opt_conexp_pc$par 
    
    
    ##
    ## Negative clairvoyance data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_nc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_nc <- optim(par = mlsp_par_conexp_nc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_nc$dist, 
                                z = dat_rs_an_nc$pos, 
                                n = dat_rs_an_nc$n, 
                                db1 = dat_rs_an_nc$year >= 2014,
                                db2 = dat_rs_an_nc$year >= 2015,
                                db3 = dat_rs_an_nc$year >= 2016, 
                                db4 = dat_rs_an_nc$year >= 2017,
                                db5 = dat_rs_an_nc$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_nc <- mlsp_opt_conexp_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_ncpc <- c(b = 0.1, x100 = -5, c = 14, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_ncpc <- optim(par = mlsp_par_conexp_ncpc, 
                                  fn = mlsp_fun_conexp, 
                                  x = dat_rs_an_ncpc$dist, 
                                  z = dat_rs_an_ncpc$pos, 
                                  n = dat_rs_an_ncpc$n, 
                                  db1 = dat_rs_an_ncpc$year >= 2014,
                                  db2 = dat_rs_an_ncpc$year >= 2015,
                                  db3 = dat_rs_an_ncpc$year >= 2016, 
                                  db4 = dat_rs_an_ncpc$year >= 2017,
                                  db5 = dat_rs_an_ncpc$year >= 2018,
                                  control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_ncpc <- mlsp_opt_conexp_ncpc$par 
    
    
    #               
    # Hill function 
    #               
    
    # Setup hill function
    mlsp_fun_hill <- function(par, x, z, n, db1, db2, db3, db4, db5){
      mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), par[1], 
                   (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2])/
                              (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2]))))
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_od <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1) 
    
    # Model fit
    mlsp_opt_hill_od <- optim(par = mlsp_par_hill_od, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_od$dist, 
                              z = dat_rs_an_od$pos, 
                              n = dat_rs_an_od$n, 
                              db1 = dat_rs_an_od$year >= 2014, 
                              db2 = dat_rs_an_od$year >= 2015,
                              db3 = dat_rs_an_od$year >= 2016, 
                              db4 = dat_rs_an_od$year >= 2017, 
                              db5 = dat_rs_an_od$year >= 2018, 
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_od <- mlsp_opt_hill_od$par
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_pc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_pc <- optim(par = mlsp_par_hill_pc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_pc$dist, 
                              z = dat_rs_an_pc$pos, 
                              n = dat_rs_an_pc$n, 
                              db1 = dat_rs_an_pc$year >= 2014, 
                              db2 = dat_rs_an_pc$year >= 2015, 
                              db3 = dat_rs_an_pc$year >= 2016, 
                              db4 = dat_rs_an_pc$year >= 2017,
                              db5 = dat_rs_an_pc$year >= 2018,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_pc <- mlsp_opt_hill_pc$par 
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_nc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_nc <- optim(par = mlsp_par_hill_nc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_nc$dist, 
                              z = dat_rs_an_nc$pos, 
                              n = dat_rs_an_nc$n, 
                              db1 = dat_rs_an_nc$year >= 2014,
                              db2 = dat_rs_an_nc$year >= 2015, 
                              db3 = dat_rs_an_nc$year >= 2016, 
                              db4 = dat_rs_an_nc$year >= 2017,
                              db5 = dat_rs_an_nc$year >= 2018,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_nc <- mlsp_opt_hill_nc$par 
    
    
    ##
    ## Positive carry-over & negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_ncpc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_ncpc <- optim(par = mlsp_par_hill_ncpc, 
                                fn = mlsp_fun_hill, 
                                x = dat_rs_an_ncpc$dist, 
                                z = dat_rs_an_ncpc$pos, 
                                n = dat_rs_an_ncpc$n, 
                                db1 = dat_rs_an_ncpc$year >= 2014, 
                                db2 = dat_rs_an_ncpc$year >= 2015, 
                                db3 = dat_rs_an_ncpc$year >= 2016, 
                                db4 = dat_rs_an_ncpc$year >= 2017,
                                db5 = dat_rs_an_ncpc$year >= 2018,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_ncpc <- mlsp_opt_hill_ncpc$par 
    
    accuracy$log_od[l] = mlsp_coefs_log_od["c"] - sim_par["c"]  
    accuracy$log_pc[l] = mlsp_coefs_log_pc["c"] - sim_par["c"]  
    accuracy$log_nc[l] = mlsp_coefs_log_nc["c"] - sim_par["c"]  
    accuracy$log_ncpc[l] = mlsp_coefs_log_ncpc["c"] - sim_par["c"]  
    
    accuracy$conexp_od[l] = mlsp_coefs_conexp_od["c"] - sim_par["c"]  
    accuracy$conexp_pc[l] = mlsp_coefs_conexp_pc["c"] - sim_par["c"]  
    accuracy$conexp_nc[l] = mlsp_coefs_conexp_nc["c"] - sim_par["c"]  
    accuracy$conexp_ncpc[l] = mlsp_coefs_conexp_ncpc["c"] - sim_par["c"]  
    
    accuracy$hill_od[l] = mlsp_coefs_hill_od["c"] - sim_par["c"]  
    accuracy$hill_pc[l] = mlsp_coefs_hill_pc["c"] - sim_par["c"]  
    accuracy$hill_nc[l] = mlsp_coefs_hill_nc["c"] - sim_par["c"]  
    accuracy$hill_ncpc[l] = mlsp_coefs_hill_ncpc["c"] - sim_par["c"] 
    
    # Calculate progress percentage, to keep track of how far along the simulation is
    pro_perc <- ((k-1)*n_sim + l)*((100/(nrow(dat_accuracy)*n_sim)))
    print(pro_perc)
    
  }
  
  dat_accuracy$accuracy_log_od_mean[k] = mean(accuracy$log_od)
  dat_accuracy$accuracy_log_od_sd[k] = sd(accuracy$log_od)
  dat_accuracy$accuracy_log_pc_mean[k] = mean(accuracy$log_pc)
  dat_accuracy$accuracy_log_pc_sd[k] = sd(accuracy$log_pc)
  dat_accuracy$accuracy_log_nc_mean[k] = mean(accuracy$log_nc)
  dat_accuracy$accuracy_log_nc_sd[k] = sd(accuracy$log_nc)
  dat_accuracy$accuracy_log_ncpc_mean[k] = mean(accuracy$log_ncpc)
  dat_accuracy$accuracy_log_ncpc_sd[k] = sd(accuracy$log_ncpc)
  dat_accuracy$accuracy_conexp_od_mean[k] = mean(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_od_sd[k] = sd(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_pc_mean[k] = mean(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_pc_sd[k] = sd(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_nc_mean[k] = mean(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_nc_sd[k] = sd(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_ncpc_mean[k] = mean(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_conexp_ncpc_sd[k] = sd(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_hill_od_mean[k] = mean(accuracy$hill_od)
  dat_accuracy$accuracy_hill_od_sd[k] = sd(accuracy$hill_od)
  dat_accuracy$accuracy_hill_pc_mean[k] = mean(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_pc_sd[k] = sd(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_nc_mean[k] = mean(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_nc_sd[k] = sd(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_ncpc_mean[k] = mean(accuracy$hill_ncpc)
  dat_accuracy$accuracy_hill_ncpc_sd[k] = sd(accuracy$hill_ncpc)
}

write.xlsx(dat_accuracy, "C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy_conexp.xlsx", sheetName = "Years")
dat_accuracy_seasons <- read.xlsx("C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy_conexp.xlsx", sheetName = "Years")


###############
### SEASONS ###
###############

###############################
##                           ##
## Data reading and morphing ##
##                           ##
###############################

dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter((lon != 0) | (lat != 0))

dat <- dat %>% 
  mutate(season = ifelse((dat$year == 2013) | (dat$year == 2014 & dat$month %in% 1:3), 1, 
                         ifelse((dat$year == 2014 & dat$month %in% 4:12) | (dat$year == 2015 & dat$month %in% 1:3), 2, 
                                ifelse((dat$year == 2015 & dat$month %in% 4:12) | (dat$year == 2016 & dat$month %in% 1:3), 3, 
                                       ifelse((dat$year == 2016 & dat$month %in% 4:12) | (dat$year == 2017 & dat$month %in% 1:3), 4, 5)))))


# Set number of seasons, used for length of loops and such
seasons <- unique(dat$season)

gallipoli <- c(lon=17.992615,lat=40.055851)


################
##            ##
## Simulation ##
##            ##
################

cvector <- seq(from = 10, to = 10, by = 1)
x100vector <- seq(from = -20, to = -20, by = 5)

n_sim <- 2

dat_accuracy <- tibble("c" = sort(rep(cvector, times = length(x100vector))), 
                       "x100" = rep(x100vector, times = length(cvector)), 
                       "accuracy_log_od_mean" = NA,
                       "accuracy_log_od_sd" = NA,
                       "accuracy_log_pc_mean" = NA,
                       "accuracy_log_pc_sd" = NA,
                       "accuracy_log_nc_mean" = NA,
                       "accuracy_log_nc_sd" = NA,
                       "accuracy_log_ncpc_mean" = NA,
                       "accuracy_log_ncpc_sd" = NA,
                       "accuracy_conexp_od_mean" = NA,
                       "accuracy_conexp_od_sd" = NA,
                       "accuracy_conexp_pc_mean" = NA,
                       "accuracy_conexp_pc_sd" = NA,
                       "accuracy_conexp_nc_mean" = NA,
                       "accuracy_conexp_nc_sd" = NA,
                       "accuracy_conexp_ncpc_mean" = NA,
                       "accuracy_conexp_ncpc_sd" = NA,
                       "accuracy_hill_od_mean" = NA,
                       "accuracy_hill_od_sd" = NA,
                       "accuracy_hill_pc_mean" = NA,
                       "accuracy_hill_pc_sd" = NA,
                       "accuracy_hill_nc_mean" = NA,
                       "accuracy_hill_nc_sd" = NA,
                       "accuracy_hill_ncpc_mean" = NA,
                       "accuracy_hill_ncpc_sd" = NA)

for(k in 1:nrow(dat_accuracy)){
  sim_par <- c("a" = 0.08, "x100" = dat_accuracy$x100[k], "c" = dat_accuracy$c[k], "theta" = 1)
  
  accuracy <- list("log_od" = numeric(n_sim),
                   "log_pc" = numeric(n_sim),
                   "log_nc" = numeric(n_sim),
                   "log_ncpc" = numeric(n_sim),
                   "conexp_od" = numeric(n_sim),
                   "conexp_pc" = numeric(n_sim),
                   "conexp_nc" = numeric(n_sim),
                   "conexp_ncpc" = numeric(n_sim),
                   "hill_od" = numeric(n_sim),
                   "hill_pc" = numeric(n_sim),
                   "hill_nc" = numeric(n_sim),
                   "hill_ncpc" = numeric(n_sim))
  
  for(l in 1:n_sim){

    #####  Generating Sampling Data #####
    
    # Generate sampling dataset
    ## Setup new sampled dataset
    dat_sam <- dat %>% 
      dplyr::select(lon, lat, season) %>% 
      mutate("dist" = NA)
    
    # Calculate distance of sampled points to Gallipoli
    dat_sam <- dat_sam %>% 
      mutate(dist = round(distHaversine(tibble(lon = lon, lat = lat), 
                                        c(gallipoli[1], gallipoli[2]))/1000))
    
    dat_sam2 <- dat_sam %>% 
      dplyr::select(season, dist) %>% 
      mutate(n = 1)
    
    dat_newsam <- as_tibble(melt(tapply(dat_sam2$n, list("season" = dat_sam2$season, # Find the number of samples in each distance circle, in each year
                                                         "dist" = dat_sam2$dist), sum))) %>% 
      mutate("pos" = NA)
    
    names(dat_newsam)[3] <- "n" # Name the number of samples column
    
    # Get sample results based on probability for this location and a binomial distribution
    for(i in 1:nrow(dat_newsam)){
      db1 = dat_newsam$season[i] >= 2 # db1 is TRUE for 2 and above, otherwise FALSE
      db2 = dat_newsam$season[i] >= 3 # db2 is TRUE for 3 and above, otherwise FALSE
      db3 = dat_newsam$season[i] >= 4 
      db4 = dat_newsam$season[i] >= 5
      
      # Deterministic model
      prob <- ifelse((1/(exp(-sim_par["a"]*(sim_par["x100"] + (1*db1 + 1*db2 + 1*db3 + 1*db4)*sim_par["c"]))))*exp(-sim_par["a"]*dat_newsam$dist[i]) < 0.999,
                     (1/(exp(-sim_par["a"]*(sim_par["x100"] + (1*db1 + 1*db2 + 1*db3 + 1*db4)*sim_par["c"]))))*exp(-sim_par["a"]*dat_newsam$dist[i]), 0.999)
      
      # Get result from beta-binomial distribution
      dat_newsam$pos[i] <- rbetabinom(n = 1, prob = prob, size = dat_newsam$n[i], theta = sim_par["theta"])
    }
    
    dat_newsam <- dat_newsam %>% 
      mutate(prop = pos/n)
      
    
    ######  Analysis of sampling ######
    
    ###                           
    ### Rate of spread estimation 
    ###       
    
    #
    # Logistic function
    #
    
    # Setup logistic function: y = 1 / (1 + exp(a * (x - t*c))) where c is dependent on the year
    mlsp_fun_log <- function(par, x, z, n, db1, db2, db3, db4){
      mu <- (1/(1 + exp(par[1] *                                                                 # par[1] is parameter 'a' (declining slope), par[2] is parameter 'x50', par[3] is parameter 'c',
                          (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))))   # 'dbx' can be TRUE or FALSE, the number of1's dependent on the 'dbx''s is the value of 't'.
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
      return(nll)
    }
    
    ##
    ## Original data rate of spread analysis
    ##
    
    dat_rs_an_od <- dat_newsam %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_od <- c(a = 0.05, x0 = -10, c = 10, theta = 1) # Starting parameters
    
    # Model fit
    mlsp_opt_log_od <- optim(par = mlsp_par_log_od, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_od$dist, 
                             z = dat_rs_an_od$pos, 
                             n = dat_rs_an_od$n, 
                             db1 = dat_rs_an_od$season >= 2, # db1 is TRUE for 2014 and above
                             db2 = dat_rs_an_od$season >= 3,  
                             db3 = dat_rs_an_od$season >= 4, 
                             db4 = dat_rs_an_od$season >= 5, 
                             control = list(maxit = 10000), method = "Nelder-Mead") 
    mlsp_coefs_log_od <- mlsp_opt_log_od$par # Set parameters
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_pc <- dat_newsam %>% # Set new data frame to create cumulative data
      replace_na(list(n = 0, pos = 0, prop = 0)) %>%   # Set NA values to 0, so positives can later be added from earlier seasons
      arrange(season, dist)
    
    for(i in 1:nrow(dat_rs_an_pc)){
      sn = dat_rs_an_pc$season[i] # The season of this datapoint
      dst = dat_rs_an_pc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn + 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn + 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_pc$n[i] > 0, ndup/dat_rs_an_pc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous season
      pos_prev <- dat_rs_an_pc$pos[dat_rs_an_pc$season == (dat_rs_an_pc$season[i] - 1) &
                                     dat_rs_an_pc$dist == dat_rs_an_pc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous season (this value can be 0)
        dat_rs_an_pc$pos[i] = dat_rs_an_pc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous season in this distance circle
        dat_rs_an_pc$n[i] = dat_rs_an_pc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_pc <- dat_rs_an_pc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_pc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_pc <- optim(par = mlsp_par_log_pc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_pc$dist, 
                             z = dat_rs_an_pc$pos, 
                             n = dat_rs_an_pc$n, 
                             db1 = dat_rs_an_pc$season >= 2, 
                             db2 = dat_rs_an_pc$season >= 3,
                             db3 = dat_rs_an_pc$season >= 4, 
                             db4 = dat_rs_an_pc$season >= 5,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_pc <- mlsp_opt_log_pc$par
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    dat_rs_an_nc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(season), dist)
    
    for(i in 1:nrow(dat_rs_an_nc)){
      sn = dat_rs_an_nc$season[i] # The season of this datapoint
      dst = dat_rs_an_nc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn - 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_nc$n[i] > 0, ndup/dat_rs_an_nc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_nc$n[i] - dat_rs_an_nc$pos[i]
      
      # Add the number of negatives of this season to the previous season, in the same distance circle, corrected for duplicate samples
      dat_rs_an_nc$n[dat_rs_an_nc$season == (sn - 1) & dat_rs_an_nc$dist == dst] <-  
        dat_rs_an_nc$n[dat_rs_an_nc$season == (sn - 1) & dat_rs_an_nc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_nc <- dat_rs_an_nc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_nc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_nc <- optim(par = mlsp_par_log_nc, 
                             fn = mlsp_fun_log, 
                             x = dat_rs_an_nc$dist, 
                             z = dat_rs_an_nc$pos, 
                             n = dat_rs_an_nc$n, 
                             db1 = dat_rs_an_nc$season >= 2, 
                             db2 = dat_rs_an_nc$season >= 3, 
                             db3 = dat_rs_an_nc$season >= 4, 
                             db4 = dat_rs_an_nc$season >= 5,
                             control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_nc <- mlsp_opt_log_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data rate of spread analysis
    ##
    
    dat_rs_an_ncpc <- dat_newsam %>% 
      replace_na(list(n = 0, pos = 0, prop = 0)) %>% 
      arrange(desc(season), dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      sn = dat_rs_an_ncpc$season[i] # The season of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn - 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn - 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of negatives in the distance circle of the current datapoint
      neg_curr <- dat_rs_an_ncpc$n[i] - dat_rs_an_ncpc$pos[i]
      
      # Add the number of negatives of this season to the previous season, in the same distance circle, corrected for duplicate samples
      dat_rs_an_ncpc$n[dat_rs_an_ncpc$season == (sn - 1) & dat_rs_an_ncpc$dist == dst] <-  
        dat_rs_an_ncpc$n[dat_rs_an_ncpc$season == (sn - 1) & dat_rs_an_ncpc$dist == dst] + round((1 - p_ndup)*neg_curr)      
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      arrange(season, dist)
    
    for(i in 1:nrow(dat_rs_an_ncpc)){
      sn = dat_rs_an_ncpc$season[i] # The season of this datapoint
      dst = dat_rs_an_ncpc$dist[i] # The distance circle of this datapoint
      
      # Find the number of samples in this distance circle that were also sampled the season before
      ndup <- length(which(dat_sam$lon[dat_sam$season == sn] == dat_sam$lon[dat_sam$season == (sn + 1)] &
                             dat_sam$lat[dat_sam$season == sn] == dat_sam$lat[dat_sam$season == (sn + 1)] &
                             dat_sam$dist == dst))
      
      # Set the proportion of samples that are resampled
      p_ndup <- ifelse(dat_rs_an_ncpc$n[i] > 0, ndup/dat_rs_an_ncpc$n[i], 0)
      
      # Find the number of positives in the distance circle of the current datapoint, sampled in the previous season
      pos_prev <- dat_rs_an_ncpc$pos[dat_rs_an_ncpc$season == (dat_rs_an_ncpc$season[i] - 1) &
                                       dat_rs_an_ncpc$dist == dat_rs_an_ncpc$dist[i]]
      
      if(length(pos_prev) > 0){ # If pos_prev is not empty, so there is a value of positives of this distance circle in a previous season (this value can be 0)
        dat_rs_an_ncpc$pos[i] = dat_rs_an_ncpc$pos[i] + round((1 - p_ndup)*pos_prev) # Add the positives to pos. The number of positives to carry over is the proportion of new (not resampled) samples times the number of positives from the previous season in this distance circle
        dat_rs_an_ncpc$n[i] = dat_rs_an_ncpc$n[i] + round((1 - p_ndup)*pos_prev) # Add the number of positives to the number of samples
      }
    }
    
    dat_rs_an_ncpc <- dat_rs_an_ncpc %>% 
      mutate(prop = pos/n) %>% 
      drop_na()
    
    # Setup parameter starting values
    mlsp_par_log_ncpc <- c(a = 0.05, x0 = -10, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_log_ncpc <- optim(par = mlsp_par_log_ncpc, 
                               fn = mlsp_fun_log, 
                               x = dat_rs_an_ncpc$dist, 
                               z = dat_rs_an_ncpc$pos, 
                               n = dat_rs_an_ncpc$n, 
                               db1 = dat_rs_an_ncpc$season >= 2, 
                               db2 = dat_rs_an_ncpc$season >= 3, 
                               db3 = dat_rs_an_ncpc$season >= 4, 
                               db4 = dat_rs_an_ncpc$season >= 5,
                               control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_log_ncpc <- mlsp_opt_log_ncpc$par 
    
    
    #                                           
    # Constrained Negative Exponential function 
    #     
    
    # Setup CNE function
    mlsp_fun_conexp <- function(par, x, z, n, db1, db2, db3, db4){
      mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))) * exp(-par[1] * x) < 0.999,
                   (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))) * exp(-par[1] * x), 0.999) 
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_od <- c(b = 0.1, x100 = -5, c = 12, theta = 1) 
    
    # Model fit
    mlsp_opt_conexp_od <- optim(par = mlsp_par_conexp_od, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_od$dist, 
                                z = dat_rs_an_od$pos, 
                                n = dat_rs_an_od$n, 
                                db1 = dat_rs_an_od$season >= 2, 
                                db2 = dat_rs_an_od$season >= 3, 
                                db3 = dat_rs_an_od$season >= 4, 
                                db4 = dat_rs_an_od$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_od <- mlsp_opt_conexp_od$par 
    
    
    ##
    ## Positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_pc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_pc <- optim(par = mlsp_par_conexp_pc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_pc$dist, 
                                z = dat_rs_an_pc$pos, 
                                n = dat_rs_an_pc$n, 
                                db1 = dat_rs_an_pc$season >= 2, 
                                db2 = dat_rs_an_pc$season >= 3, 
                                db3 = dat_rs_an_pc$season >= 4, 
                                db4 = dat_rs_an_pc$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_pc <- mlsp_opt_conexp_pc$par 
    
    
    ##
    ## Negative clairvoyance data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_nc <- c(b = 0.1, x100 = -5, c = 12, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_nc <- optim(par = mlsp_par_conexp_nc, 
                                fn = mlsp_fun_conexp, 
                                x = dat_rs_an_nc$dist, 
                                z = dat_rs_an_nc$pos, 
                                n = dat_rs_an_nc$n, 
                                db1 = dat_rs_an_nc$season >= 2,
                                db2 = dat_rs_an_nc$season >= 3,
                                db3 = dat_rs_an_nc$season >= 4, 
                                db4 = dat_rs_an_nc$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_conexp_nc <- mlsp_opt_conexp_nc$par 
    
    
    ##
    ## Negative clairvoyance & positive carry-over data speed of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_conexp_ncpc <- c(b = 0.1, x100 = -5, c = 14, theta = 1)
    
    # Model fit
    mlsp_opt_conexp_ncpc <- optim(par = mlsp_par_conexp_ncpc, 
                                  fn = mlsp_fun_conexp, 
                                  x = dat_rs_an_ncpc$dist, 
                                  z = dat_rs_an_ncpc$pos, 
                                  n = dat_rs_an_ncpc$n, 
                                  db1 = dat_rs_an_ncpc$season >= 2,
                                  db2 = dat_rs_an_ncpc$season >= 3,
                                  db3 = dat_rs_an_ncpc$season >= 4, 
                                  db4 = dat_rs_an_ncpc$season >= 5,
                                  control = list(maxit = 5000), method = "Nelder-Mead")
    mlsp_coefs_conexp_ncpc <- mlsp_opt_conexp_ncpc$par 
    
    
    #               
    # Hill function 
    #               
    
    # Setup hill function
    mlsp_fun_hill <- function(par, x, z, n, db1, db2, db3, db4){
      mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), par[1], 
                   (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2])/
                              (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2]))))
      nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) 
      return(nll)
    }
    
    ##
    ## Original data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_od <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1) 
    
    # Model fit
    mlsp_opt_hill_od <- optim(par = mlsp_par_hill_od, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_od$dist, 
                              z = dat_rs_an_od$pos, 
                              n = dat_rs_an_od$n, 
                              db1 = dat_rs_an_od$season >= 2, 
                              db2 = dat_rs_an_od$season >= 3,
                              db3 = dat_rs_an_od$season >= 4, 
                              db4 = dat_rs_an_od$season >= 5, 
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_od <- mlsp_opt_hill_od$par
    
    
    ##
    ## Positive carry-over data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_pc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_pc <- optim(par = mlsp_par_hill_pc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_pc$dist, 
                              z = dat_rs_an_pc$pos, 
                              n = dat_rs_an_pc$n, 
                              db1 = dat_rs_an_pc$season >= 2, 
                              db2 = dat_rs_an_pc$season >= 3, 
                              db3 = dat_rs_an_pc$season >= 4, 
                              db4 = dat_rs_an_pc$season >= 5,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_pc <- mlsp_opt_hill_pc$par 
    
    
    ##
    ## Negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_nc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_nc <- optim(par = mlsp_par_hill_nc, 
                              fn = mlsp_fun_hill, 
                              x = dat_rs_an_nc$dist, 
                              z = dat_rs_an_nc$pos, 
                              n = dat_rs_an_nc$n, 
                              db1 = dat_rs_an_nc$season >= 2,
                              db2 = dat_rs_an_nc$season >= 3, 
                              db3 = dat_rs_an_nc$season >= 4, 
                              db4 = dat_rs_an_nc$season >= 5,
                              control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_nc <- mlsp_opt_hill_nc$par 
    
    
    ##
    ## Positive carry-over & negative clairvoyance data rate of spread analysis
    ##
    
    # Setup parameter starting values
    mlsp_par_hill_ncpc <- c(a = 0.99, n = 2, h = 30, c = 10, theta = 1)
    
    # Model fit
    mlsp_opt_hill_ncpc <- optim(par = mlsp_par_hill_ncpc, 
                                fn = mlsp_fun_hill, 
                                x = dat_rs_an_ncpc$dist, 
                                z = dat_rs_an_ncpc$pos, 
                                n = dat_rs_an_ncpc$n, 
                                db1 = dat_rs_an_ncpc$season >= 2, 
                                db2 = dat_rs_an_ncpc$season >= 3, 
                                db3 = dat_rs_an_ncpc$season >= 4, 
                                db4 = dat_rs_an_ncpc$season >= 5,
                                control = list(maxit = 5000), method = "Nelder-Mead") 
    mlsp_coefs_hill_ncpc <- mlsp_opt_hill_ncpc$par 
    
    accuracy$log_od[l] = mlsp_coefs_log_od["c"] - sim_par["c"]  
    accuracy$log_pc[l] = mlsp_coefs_log_pc["c"] - sim_par["c"]  
    accuracy$log_nc[l] = mlsp_coefs_log_nc["c"] - sim_par["c"]  
    accuracy$log_ncpc[l] = mlsp_coefs_log_ncpc["c"] - sim_par["c"]  
    
    accuracy$conexp_od[l] = mlsp_coefs_conexp_od["c"] - sim_par["c"]  
    accuracy$conexp_pc[l] = mlsp_coefs_conexp_pc["c"] - sim_par["c"]  
    accuracy$conexp_nc[l] = mlsp_coefs_conexp_nc["c"] - sim_par["c"]  
    accuracy$conexp_ncpc[l] = mlsp_coefs_conexp_ncpc["c"] - sim_par["c"]  
    
    accuracy$hill_od[l] = mlsp_coefs_hill_od["c"] - sim_par["c"]  
    accuracy$hill_pc[l] = mlsp_coefs_hill_pc["c"] - sim_par["c"]  
    accuracy$hill_nc[l] = mlsp_coefs_hill_nc["c"] - sim_par["c"]  
    accuracy$hill_ncpc[l] = mlsp_coefs_hill_ncpc["c"] - sim_par["c"] 
    
    pro_perc <- ((k-1)*n_sim + l)*((100/(nrow(dat_accuracy)*n_sim)))
    print(pro_perc)
    
    
  }
  
  
  dat_accuracy$accuracy_log_od_mean[k] = mean(accuracy$log_od)
  dat_accuracy$accuracy_log_od_sd[k] = sd(accuracy$log_od)
  dat_accuracy$accuracy_log_pc_mean[k] = mean(accuracy$log_pc)
  dat_accuracy$accuracy_log_pc_sd[k] = sd(accuracy$log_pc)
  dat_accuracy$accuracy_log_nc_mean[k] = mean(accuracy$log_nc)
  dat_accuracy$accuracy_log_nc_sd[k] = sd(accuracy$log_nc)
  dat_accuracy$accuracy_log_ncpc_mean[k] = mean(accuracy$log_ncpc)
  dat_accuracy$accuracy_log_ncpc_sd[k] = sd(accuracy$log_ncpc)
  dat_accuracy$accuracy_conexp_od_mean[k] = mean(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_od_sd[k] = sd(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_pc_mean[k] = mean(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_pc_sd[k] = sd(accuracy$conexp_pc)
  dat_accuracy$accuracy_conexp_nc_mean[k] = mean(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_nc_sd[k] = sd(accuracy$conexp_nc)
  dat_accuracy$accuracy_conexp_ncpc_mean[k] = mean(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_conexp_ncpc_sd[k] = sd(accuracy$conexp_ncpc)
  dat_accuracy$accuracy_hill_od_mean[k] = mean(accuracy$hill_od)
  dat_accuracy$accuracy_hill_od_sd[k] = sd(accuracy$hill_od)
  dat_accuracy$accuracy_hill_pc_mean[k] = mean(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_pc_sd[k] = sd(accuracy$hill_pc)
  dat_accuracy$accuracy_hill_nc_mean[k] = mean(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_nc_sd[k] = sd(accuracy$hill_nc)
  dat_accuracy$accuracy_hill_ncpc_mean[k] = mean(accuracy$hill_ncpc)
  dat_accuracy$accuracy_hill_ncpc_sd[k] = sd(accuracy$hill_ncpc)
  
}

write.xlsx(dat_accuracy, "C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy_conexp.xlsx", sheetName = "Seasons", append = TRUE)
dat_accuracy_seasons <- read.xlsx("C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy_conexp.xlsx", sheetName = "Seasons")
