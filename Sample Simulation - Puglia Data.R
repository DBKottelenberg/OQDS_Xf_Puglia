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
                       "accuracy_conexp_od_mean" = NA,
                       "accuracy_conexp_od_sd" = NA)

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
    ## rate of spread analysis
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
    ## speed of spread analysis
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
    
    
    accuracy$log_od[l] = mlsp_coefs_log_od["c"] - sim_par["c"]  
    accuracy$conexp_od[l] = mlsp_coefs_conexp_od["c"] - sim_par["c"]  
    
    # Calculate progress percentage, to keep track of how far along the simulation is
    pro_perc <- ((k-1)*n_sim + l)*((100/(nrow(dat_accuracy)*n_sim)))
    print(pro_perc)
    
  }
  
  dat_accuracy$accuracy_log_od_mean[k] = mean(accuracy$log_od)
  dat_accuracy$accuracy_log_od_sd[k] = sd(accuracy$log_od)*0.9 # Account for one more degree of freedom than the `sd()` function accounts for

  dat_accuracy$accuracy_conexp_od_mean[k] = mean(accuracy$conexp_od)
  dat_accuracy$accuracy_conexp_od_sd[k] = sd(accuracy$conexp_od)*0.9
}

write.xlsx(dat_accuracy, "C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy.xlsx", sheetName = "Years")
dat_accuracy_years <- read.xlsx("C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R/Simulation/dat_accuracy.xlsx", sheetName = "Years")
