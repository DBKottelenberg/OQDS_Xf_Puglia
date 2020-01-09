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

# Set ggplot graph theme
theme_set(theme_bw())


#######################################
#######################################
######                           ######
######           YEARS           ######
######                           ######
#######################################
#######################################

###############################
##                           ##
## Data reading and morphing ##
##                           ##
###############################

# Read in data
dat_base <- read_xlsx("xylelladata.xlsx")
dat <- dat_base

# Set date factor correctly, drop na and incorrect lon/lat values
dat <- dat_base
dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter(lon != 0 & lat != 0)

# Set number of years, used for length of loops and such
years <- unique(dat$year)

# Assumed origin
gallipoli <- c(lon=17.992615,lat=40.055851)

# Calculate distance of data points to origin
dat <- dat %>%
  mutate(dist = round((distHaversine(tibble(lon = dat$lon, lat = dat$lat), 
                                     c(gallipoli[1], gallipoli[2])))/1000))

# Find the maximum distance measured for the distance circels and the maximum distance where a positive is found for visualization
max_dist <- max(dat$dist)
max_pos <- max(dat[dat$result == 1,]$dist)

# Split the data for each year
dat <- dat %>% 
  split(dat$year)


###################
##               ##
## Model fitting ##
##               ##
###################

#################
# Original data #
#################

#
# Data morphing and distance circle creation
#

# Set distance circles
dc <- 1:(max_dist) # Number of distance circles is the maximum distance a measurement is done.
pos <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for number of positives in each dc
n <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for total number of measurements in each dc

dat_2 <- list() # Prepare data list

for(i in 1:length(years)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j] <- sum((dat[[i]]$result[dat[[i]]$dist >= j & (dat[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j] <- length((dat[[i]]$result[dat[[i]]$dist >= j & (dat[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_2[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
}

# Combine years and create year column for year 
dat_3 <- list()
for(i in 1:length(years)){
  dat_2[[i]] <- dat_2[[i]] %>% 
    mutate(year = (2012 + i))
  dat_3 <- rbind(dat_3, dat_2[[i]])
}

# Plot dat_3
ggplot() +
  geom_point(data = dat_3[dat_3$year == 2013,], aes(x = dist, y = prop, color = "2013"), 
             shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2014,], aes(x = dist, y = prop, color = "2014"), 
             shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2015,], aes(x = dist, y = prop, color = "2015"), 
             shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2016,], aes(x = dist, y = prop, color = "2016"), 
             shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2017,], aes(x = dist, y = prop, color = "2017"), 
             shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2018,], aes(x = dist, y = prop, color = "2018"), 
             shape = 1) +
  labs(
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C1. Data points of the original data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", 
                                               "magenta3")) +
  scale_size_manual(values = 1.05) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function: y = 1 / (1 + exp(a * (x - c*t))) where c is dependent on the year
mlsp_fun_od <- function(par, x, z, n, db1, db2, db3, db4, db5){ # maximum likelihood speed of spread function for original data (mlsp_fun_od)
  mu <- (1/(1 + exp(par[1] *                                                               # par[1] is parameter 'a' (declining slope)
                      (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3]))))) # par[2] is parameter x50, par[3] is parameter c, dbx can be TRUE or FALSE, 
                                                                                           # the number of 1's dependent on the dbx's is the value of t.
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_od <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_od <- optim(par = mlsp_par_od, 
                     fn = mlsp_fun_od, 
                     x = dat_3$dist, 
                     z = dat_3$pos, 
                     n = dat_3$n, 
                     db1 = dat_3$year >= 2014, # db1 is TRUE for 2014 and above
                     db2 = dat_3$year >= 2015, # db2 is TRUE for 2015 and above
                     db3 = dat_3$year >= 2016, 
                     db4 = dat_3$year >= 2017, 
                     db5 = dat_3$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_od <- mlsp_opt_od$par # Set parameters

mlsp_coefs_od

# Plot the lines
ggplot() +
  geom_point(data = dat_3[dat_3$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 0*mlsp_coefs_od[3]))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 1*mlsp_coefs_od[3]))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 2*mlsp_coefs_od[3]))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 3*mlsp_coefs_od[3]))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 4*mlsp_coefs_od[3]))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 5*mlsp_coefs_od[3]))))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C2. Logistic function through the original data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_od <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_od <- function(par, x, z, n, x50, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_od <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_od <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_od = c(mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4]) # Starting parameters are the estimated parameters of the optimized function
avec_od = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_od = numeric(100) # Prepare profile of likelihood values

pars_x50_od = c(mlsp_coefs_od[1], mlsp_coefs_od[3], mlsp_coefs_od[4])
x50vec_od = seq(-120, 100, length = 100)
x50prof_od = numeric(100)

pars_c_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[4])
cvec_od = seq(0.1, 50, length = 100)
cprof_od = numeric(100)

pars_theta_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3])
thetavec_od = seq(0.01, 10, length = 100)
thetaprof_od = numeric(100)

# Optimize original function without fixed parameters
pars_2_od <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters

opt_log_bbinom_od <- optim(par = pars_2_od, 
                           fn = mlsp_fun_od, 
                           x = dat_3$dist, 
                           z = dat_3$pos, 
                           n = dat_3$n,
                           db1 = dat_3$year >= 2014, 
                           db2 = dat_3$year >= 2015, 
                           db3 = dat_3$year >= 2016, 
                           db4 = dat_3$year >= 2017, 
                           db5 = dat_3$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_od <- list(rep(list(), times = 100))
x50_coefs_vec_od <- list(rep(list(), times = 100))
c_coefs_vec_od <- list(rep(list(), times = 100))
theta_coefs_vec_od <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_od <- optim(par = pars_a_od, 
                    fn = ci_a_fun_od, 
                    x = dat_3$dist, 
                    z = dat_3$pos, 
                    n = dat_3$n, 
                    a = avec_od[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_3$year >= 2014, 
                    db2 = dat_3$year >= 2015, 
                    db3 = dat_3$year >= 2016, 
                    db4 = dat_3$year >= 2017, 
                    db5 = dat_3$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_od[[j]] = opt_a_od$par # Store parameters for control
  aprof_od[j] <- opt_a_od$value # Set likelihood value for profile
  
  opt_x50_od <- optim(par = pars_x50_od, 
                      fn = ci_x50_fun_od, 
                      x = dat_3$dist, 
                      z = dat_3$pos, 
                      n = dat_3$n, 
                      x50 = x50vec_od[j],
                      db1 = dat_3$year >= 2014, 
                      db2 = dat_3$year >= 2015, 
                      db3 = dat_3$year >= 2016, 
                      db4 = dat_3$year >= 2017, 
                      db5 = dat_3$year >= 2018, 
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_od[j] <- opt_x50_od$value
  x50_coefs_vec_od[[j]] = opt_x50_od$par
  
  opt_c_od <- optim(par = pars_c_od, 
                    fn = ci_c_fun_od, 
                    x = dat_3$dist, 
                    z = dat_3$pos, 
                    n = dat_3$n, 
                    c = cvec_od[j],
                    db1 = dat_3$year >= 2014, 
                    db2 = dat_3$year >= 2015, 
                    db3 = dat_3$year >= 2016, 
                    db4 = dat_3$year >= 2017, 
                    db5 = dat_3$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_od[j] <- opt_c_od$value
  c_coefs_vec_od[[j]] = opt_c_od$par
  
  opt_theta_od <- optim(par = pars_theta_od, 
                        fn = ci_theta_fun_od, 
                        x = dat_3$dist, 
                        z = dat_3$pos, 
                        n = dat_3$n, 
                        theta = thetavec_od[j], 
                        db1 = dat_3$year >= 2014, 
                        db2 = dat_3$year >= 2015, 
                        db3 = dat_3$year >= 2016, 
                        db4 = dat_3$year >= 2017, 
                        db5 = dat_3$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_od[j] <- opt_theta_od$value
  theta_coefs_vec_od[[j]] = opt_theta_od$par
}

# Set the 95% confidence limits
aprof_lower_od <- aprof_od[1:which.min(aprof_od)] # Likelihood values below the best estimate
avec_lower_od <- avec_od[1:which.min(aprof_od)] # Parameter values below the best estimate
aprof_higher_od <- aprof_od[which.min(aprof_od):length(aprof_od)] # Likelihood values above the best estimate
avec_higher_od <- avec_od[which.min(aprof_od):length(aprof_od)] # Parameter values above the best estimate

## Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
l_a_ci_od <- approx(aprof_lower_od, avec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) 
## Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df
r_a_ci_od <- approx(aprof_higher_od, avec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) 

x50prof_lower_od <- x50prof_od[1:which.min(x50prof_od)]
x50vec_lower_od <- x50vec_od[1:which.min(x50prof_od)]
x50prof_higher_od <- x50prof_od[which.min(x50prof_od):length(x50prof_od)]
x50vec_higher_od <- x50vec_od[which.min(x50prof_od):length(x50prof_od)]
l_x50_ci_od <- approx(x50prof_lower_od, x50vec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_x50_ci_od <- approx(x50prof_higher_od, x50vec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

cprof_lower_od <- cprof_od[1:which.min(cprof_od)]
cvec_lower_od <- cvec_od[1:which.min(cprof_od)]
cprof_higher_od <- cprof_od[which.min(cprof_od):length(cprof_od)]
cvec_higher_od <- cvec_od[which.min(cprof_od):length(cprof_od)]
l_c_ci_od <- approx(cprof_lower_od, cvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_c_ci_od <- approx(cprof_higher_od, cvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

thetaprof_lower_od <- thetaprof_od[1:which.min(thetaprof_od)]
thetavec_lower_od <- thetavec_od[1:which.min(thetaprof_od)]
thetaprof_higher_od <- thetaprof_od[which.min(thetaprof_od):length(thetaprof_od)]
thetavec_higher_od <- thetavec_od[which.min(thetaprof_od):length(thetaprof_od)]
l_theta_ci_od <- approx(thetaprof_lower_od, thetavec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_theta_ci_od <- approx(thetaprof_higher_od, thetavec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_od_log_yr <- c(opt_log_bbinom_od$par[1], l_a_ci_od$y, r_a_ci_od$y) 
x50_ci_od_log_yr <- c(opt_log_bbinom_od$par[2], l_x50_ci_od$y, r_x50_ci_od$y)
c_ci_od_log_yr <- c(opt_log_bbinom_od$par[3], l_c_ci_od$y, r_c_ci_od$y)
theta_ci_od_log_yr <- c(opt_log_bbinom_od$par[4], l_theta_ci_od$y, r_theta_ci_od$y)

# Display the parameters and their CI's
a_ci_od_log_yr
x50_ci_od_log_yr
c_ci_od_log_yr
theta_ci_od_log_yr

##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function: y = 1/(exp(-b*(x100 + c*t)))*exp(-b*x)
mlsp_fun_od <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_od <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_od <- optim(par = mlsp_par_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, 
                     db1 = dat_3$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_3$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_3$year >= 2016, 
                     db4 = dat_3$year >= 2017, 
                     db5 = dat_3$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_od <- mlsp_opt_od$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_3[dat_3$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((0) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((0) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((1) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((1) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((2) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((2) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((3) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((3) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((4) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((4) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((5) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((5) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C3. CNE function through the original data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_od <- function(par, x, z, n, b, db1, db2, db3, db4, db5){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_od <- function(par, x, z, n, x100, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_od <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_od <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_od = c(mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_od = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_od = numeric(100) # Prepare profile of likelihood values

pars_x100_od = c(mlsp_coefs_od[1], mlsp_coefs_od[3], mlsp_coefs_od[4])
x100vec_od = seq(-120, 100, length = 100)
x100prof_od = numeric(100)

pars_c_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[4])
cvec_od = seq(0.1, 50, length = 100)
cprof_od = numeric(100)

pars_theta_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3])
thetavec_od = seq(0.01, 10, length = 100)
thetaprof_od = numeric(100)

# Optimize original function without fixed parameters
pars_2_od <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_od <- optim(par = pars_2_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n,
                           db1 = dat_3$year >= 2014, 
                           db2 = dat_3$year >= 2015, 
                           db3 = dat_3$year >= 2016, 
                           db4 = dat_3$year >= 2017, 
                           db5 = dat_3$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_od <- list(rep(list(), times = 100))
x100_coefs_vec_od <- list(rep(list(), times = 100))
c_coefs_vec_od <- list(rep(list(), times = 100))
theta_coefs_vec_od <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_od <- optim(par = pars_b_od, fn = ci_b_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, b = bvec_od[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_3$year >= 2014, 
                    db2 = dat_3$year >= 2015, 
                    db3 = dat_3$year >= 2016, 
                    db4 = dat_3$year >= 2017, 
                    db5 = dat_3$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_od[[j]] = opt_b_od$par # Store parameters for control
  # print(opt_a$convergence)
  bprof_od[j] <- opt_b_od$value # Set likelihood value for profile
  
  opt_x100_od <- optim(par = pars_x100_od, fn = ci_x100_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, x100 = x100vec_od[j],
                       db1 = dat_3$year >= 2014, 
                       db2 = dat_3$year >= 2015, 
                       db3 = dat_3$year >= 2016, 
                       db4 = dat_3$year >= 2017, 
                       db5 = dat_3$year >= 2018, 
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_od[j] <- opt_x100_od$value
  x100_coefs_vec_od[[j]] = opt_x100_od$par
  
  opt_c_od <- optim(par = pars_c_od, fn = ci_c_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, c = cvec_od[j],
                    db1 = dat_3$year >= 2014, 
                    db2 = dat_3$year >= 2015, 
                    db3 = dat_3$year >= 2016, 
                    db4 = dat_3$year >= 2017, 
                    db5 = dat_3$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_od[j] <- opt_c_od$value
  c_coefs_vec_od[[j]] = opt_c_od$par
  
  opt_theta_od <- optim(par = pars_theta_od, fn = ci_theta_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, theta = thetavec_od[j], 
                        db1 = dat_3$year >= 2014, 
                        db2 = dat_3$year >= 2015, 
                        db3 = dat_3$year >= 2016, 
                        db4 = dat_3$year >= 2017, 
                        db5 = dat_3$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_od[j] <- opt_theta_od$value
  theta_coefs_vec_od[[j]] = opt_theta_od$par
}

# Set the 95% confidence limits
bprof_lower_od <- bprof_od[1:which.min(bprof_od)] # Likelihood values below the best estimate
bvec_lower_od <- bvec_od[1:which.min(bprof_od)] # Parameter values below the best estimate
bprof_higher_od <- bprof_od[which.min(bprof_od):length(bprof_od)] # Likelihood values above the best estimate
bvec_higher_od <- bvec_od[which.min(bprof_od):length(bprof_od)] # Parameter values above the best estimate
l_b_ci_od <- approx(bprof_lower_od, bvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_od <- approx(bprof_higher_od, bvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_od <- x100prof_od[1:which.min(x100prof_od)]
x100vec_lower_od <- x100vec_od[1:which.min(x100prof_od)]
x100prof_higher_od <- x100prof_od[which.min(x100prof_od):length(x100prof_od)]
x100vec_higher_od <- x100vec_od[which.min(x100prof_od):length(x100prof_od)]
l_x100_ci_od <- approx(x100prof_lower_od, x100vec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_x100_ci_od <- approx(x100prof_higher_od, x100vec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

cprof_lower_od <- cprof_od[1:which.min(cprof_od)]
cvec_lower_od <- cvec_od[1:which.min(cprof_od)]
cprof_higher_od <- cprof_od[which.min(cprof_od):length(cprof_od)]
cvec_higher_od <- cvec_od[which.min(cprof_od):length(cprof_od)]
l_c_ci_od <- approx(cprof_lower_od, cvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_c_ci_od <- approx(cprof_higher_od, cvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

thetaprof_lower_od <- thetaprof_od[1:which.min(thetaprof_od)]
thetavec_lower_od <- thetavec_od[1:which.min(thetaprof_od)]
thetaprof_higher_od <- thetaprof_od[which.min(thetaprof_od):length(thetaprof_od)]
thetavec_higher_od <- thetavec_od[which.min(thetaprof_od):length(thetaprof_od)]
l_theta_ci_od <- approx(thetaprof_lower_od, thetavec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_theta_ci_od <- approx(thetaprof_higher_od, thetavec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_od_cne_yr <- c(opt_log_bbinom_od$par[1], l_b_ci_od$y, r_b_ci_od$y) 
x100_ci_od_cne_yr <- c(opt_log_bbinom_od$par[2], l_x100_ci_od$y, r_x100_ci_od$y)
c_ci_od_cne_yr <- c(opt_log_bbinom_od$par[3], l_c_ci_od$y, r_c_ci_od$y)
theta_ci_od_cne_yr <- c(opt_log_bbinom_od$par[4], l_theta_ci_od$y, r_theta_ci_od$y)

# Display the parameters and their CI's
b_ci_od_cne_yr
x100_ci_od_cne_yr
c_ci_od_cne_yr
theta_ci_od_cne_yr


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function: -a*((x^n)/((h + c*t)^n + x^n)) + a
mlsp_fun_od <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_od <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_od <- optim(par = mlsp_par_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, 
                     db1 = dat_3$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_3$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_3$year >= 2016, 
                     db4 = dat_3$year >= 2017, 
                     db5 = dat_3$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_od <- mlsp_opt_od$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_3[dat_3$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_3[dat_3$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^mlsp_coefs_od[2])))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C5. Hill function through the original data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_od <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                          (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_od <- function(par, x, z, n, pn, db1, db2, db3, db4, db5){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_od <- function(par, x, z, n, h, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_od <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_od <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_od = c(mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4], mlsp_coefs_od[5]) # Starting parameters are the estimated parameters of the optimized function
avec_od = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_od = numeric(100) # Prepare profile of likelihood values

pars_pn_od = c(mlsp_coefs_od[1], mlsp_coefs_od[3], mlsp_coefs_od[4], mlsp_coefs_od[5])
pnvec_od = seq(0.1, 10, length = 100)
pnprof_od = numeric(100)

pars_h_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[4], mlsp_coefs_od[5])
hvec_od = seq(5, 100, length = 100)
hprof_od = numeric(100)

pars_c_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[5])
cvec_od = seq(0.1, 50, length = 100)
cprof_od = numeric(100)

pars_theta_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4])
thetavec_od = seq(0.01, 10, length = 100)
thetaprof_od = numeric(100)

# Optimize original function without fixed parameters
pars_2_od <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_od <- optim(par = pars_2_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n,
                           db1 = dat_3$year >= 2014, 
                           db2 = dat_3$year >= 2015, 
                           db3 = dat_3$year >= 2016, 
                           db4 = dat_3$year >= 2017, 
                           db5 = dat_3$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_od <- list(rep(list(), times = 100))
pn_coefs_vec_od <- list(rep(list(), times = 100))
h_coefs_vec_od <- list(rep(list(), times = 100))
c_coefs_vec_od <- list(rep(list(), times = 100))
theta_coefs_vec_od <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_od <- optim(par = pars_a_od, fn = ci_a_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, a = avec_od[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_3$year >= 2014, 
                    db2 = dat_3$year >= 2015, 
                    db3 = dat_3$year >= 2016, 
                    db4 = dat_3$year >= 2017, 
                    db5 = dat_3$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_od[[j]] = opt_a_od$par # Store parameters for control
  aprof_od[j] <- opt_a_od$value # Set likelihood value for profile
  
  opt_pn_od <- optim(par = pars_pn_od, fn = ci_pn_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, pn = pnvec_od[j],
                     db1 = dat_3$year >= 2014, 
                     db2 = dat_3$year >= 2015, 
                     db3 = dat_3$year >= 2016, 
                     db4 = dat_3$year >= 2017, 
                     db5 = dat_3$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_od[j] <- opt_pn_od$value
  pn_coefs_vec_od[[j]] = opt_pn_od$par
  
  opt_h_od <- optim(par = pars_h_od, fn = ci_h_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, h = hvec_od[j],
                    db1 = dat_3$year >= 2014, 
                    db2 = dat_3$year >= 2015, 
                    db3 = dat_3$year >= 2016, 
                    db4 = dat_3$year >= 2017, 
                    db5 = dat_3$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_od[j] <- opt_h_od$value
  h_coefs_vec_od[[j]] = opt_h_od$par
  
  opt_c_od <- optim(par = pars_c_od, fn = ci_c_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, c = cvec_od[j],
                    db1 = dat_3$year >= 2014, 
                    db2 = dat_3$year >= 2015, 
                    db3 = dat_3$year >= 2016, 
                    db4 = dat_3$year >= 2017, 
                    db5 = dat_3$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_od[j] <- opt_c_od$value
  c_coefs_vec_od[[j]] = opt_c_od$par
  
  opt_theta_od <- optim(par = pars_theta_od, fn = ci_theta_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, theta = thetavec_od[j], 
                        db1 = dat_3$year >= 2014, 
                        db2 = dat_3$year >= 2015, 
                        db3 = dat_3$year >= 2016, 
                        db4 = dat_3$year >= 2017, 
                        db5 = dat_3$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_od[j] <- opt_theta_od$value
  theta_coefs_vec_od[[j]] = opt_theta_od$par
}

# Set the 95% confidence limits
aprof_lower_od <- aprof_od[1:which.min(aprof_od)] # Likelihood values below the best estimate
avec_lower_od <- avec_od[1:which.min(aprof_od)] # Parameter values below the best estimate
aprof_higher_od <- aprof_od[which.min(aprof_od):length(aprof_od)] # Likelihood values above the best estimate
avec_higher_od <- avec_od[which.min(aprof_od):length(aprof_od)] # Parameter values above the best estimate
l_a_ci_od <- approx(aprof_lower_od, avec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_od <- approx(aprof_higher_od, avec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_od <- pnprof_od[1:which.min(pnprof_od)]
pnvec_lower_od <- pnvec_od[1:which.min(pnprof_od)]
pnprof_higher_od <- pnprof_od[which.min(pnprof_od):length(pnprof_od)]
pnvec_higher_od <- pnvec_od[which.min(pnprof_od):length(pnprof_od)]
l_pn_ci_od <- approx(pnprof_lower_od, pnvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_pn_ci_od <- approx(pnprof_higher_od, pnvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

hprof_lower_od <- hprof_od[1:which.min(hprof_od)]
hvec_lower_od <- hvec_od[1:which.min(hprof_od)]
hprof_higher_od <- hprof_od[which.min(hprof_od):length(hprof_od)]
hvec_higher_od <- hvec_od[which.min(hprof_od):length(hprof_od)]
l_h_ci_od <- approx(hprof_lower_od, hvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_h_ci_od <- approx(hprof_higher_od, hvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

cprof_lower_od <- cprof_od[1:which.min(cprof_od)]
cvec_lower_od <- cvec_od[1:which.min(cprof_od)]
cprof_higher_od <- cprof_od[which.min(cprof_od):length(cprof_od)]
cvec_higher_od <- cvec_od[which.min(cprof_od):length(cprof_od)]
l_c_ci_od <- approx(cprof_lower_od, cvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_c_ci_od <- approx(cprof_higher_od, cvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

thetaprof_lower_od <- thetaprof_od[1:which.min(thetaprof_od)]
thetavec_lower_od <- thetavec_od[1:which.min(thetaprof_od)]
thetaprof_higher_od <- thetaprof_od[which.min(thetaprof_od):length(thetaprof_od)]
thetavec_higher_od <- thetavec_od[which.min(thetaprof_od):length(thetaprof_od)]
l_theta_ci_od <- approx(thetaprof_lower_od, thetavec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_theta_ci_od <- approx(thetaprof_higher_od, thetavec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_od_hill_yr <- c(opt_log_bbinom_od$par[1], l_a_ci_od$y, r_a_ci_od$y) 
pn_ci_od_hill_yr <- c(opt_log_bbinom_od$par[2], l_pn_ci_od$y, r_pn_ci_od$y)
h_ci_od_hill_yr <- c(opt_log_bbinom_od$par[3], l_h_ci_od$y, r_h_ci_od$y)
c_ci_od_hill_yr <- c(opt_log_bbinom_od$par[4], l_c_ci_od$y, r_c_ci_od$y)
theta_ci_od_hill_yr <- c(opt_log_bbinom_od$par[5], l_theta_ci_od$y, r_theta_ci_od$y)

# Display the parameters and their CI's
a_ci_od_hill_yr
pn_ci_od_hill_yr
h_ci_od_hill_yr
c_ci_od_hill_yr
theta_ci_od_hill_yr

# Create graph of nprof
nprof <- data.frame(pnvec = pnvec_od, pnprof = pnprof_od)
ggplot() +
  geom_line(data = nprof, aes(x = pnvec, y = pnprof)) +
  labs(x = "n value",
       y = "Negative log-likelihood",
       caption = "Figure C6. The NLL profile of parameter 'n'.") +
  theme(plot.caption = element_text(hjust = 0))



############################
# Positive carry-over data #
############################

#
# Data morphing and distance circle creation
#

dat_4 <- dat

# make positives carry over
for(i in 2:length(years)){ # For-loop through every year from 2014
  co <- which(dat_4[[i-1]]$result == 1) # Which lines of the previous year have a result of 1
  dat_4[[i]] <- rbind(dat_4[[i]], dat_4[[i-1]][co,]) # Include the lines of the previous year with a result of 1 with the current year
  co_2 <- which(dat_4[[i]]$lon == dat_4[[i-1]]$lon[co] & dat_4[[i]]$lat == dat_4[[i-1]]$lat[co] & dat_4[[i]]$season == (2012 + i)) # Which lines of this year have the same coordinates of the positives last year
  # print(co_2)
  if(length(co_2) >= 1){ # If this is not 0, then there has been a measurement at these coordinates this year as well as last year, or multiple times in a single year (with these dat_4a, only multiple measurements in 2013)
    dat_4[[i]] <- dat_4[[i]][-co_2,] # Remove all the duplicates. If the duplicate measurement from this year is a 0, it will still be a 1, because I assume that once a disease is there, it will not leave. I assume that 0 to be a false negative.
  }
}

# Set distance circles
dc <- 1:(max_dist) # Number of distance circles is the maximum distance a measurement is done.
pos <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for number of positives in each dc
n <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for total number of measurements in each dc

dat_5 <- list() # Prepare data list

for(i in 1:length(years)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j] <- sum((dat_4[[i]]$result[dat_4[[i]]$dist >= j & (dat_4[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j] <- length((dat_4[[i]]$result[dat_4[[i]]$dist >= j & (dat_4[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_5[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
}

# Combine years and create year column for year 
dat_6 <- list()
for(i in 1:length(years)){
  dat_5[[i]] <- dat_5[[i]] %>% 
    mutate(year = (2012 + i))
  dat_6 <- rbind(dat_6, dat_5[[i]])
}

# Plot the data
ggplot() +
  geom_point(data = dat_6[dat_6$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  labs(
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C7. Data points of the positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))

##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function
mlsp_fun_pc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- (1/(1 + exp(par[1] * (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3]))))) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
  return(nll)
}

# Setup parameter starting values
mlsp_par_pc <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_pc <- optim(par = mlsp_par_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, 
                     db1 = dat_6$year >= 2014, 
                     db2 = dat_6$year >= 2015, 
                     db3 = dat_6$year >= 2016, 
                     db4 = dat_6$year >= 2017, 
                     db5 = dat_6$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_pc <- mlsp_opt_pc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_6[dat_6$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 0*mlsp_coefs_pc[3]))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 1*mlsp_coefs_pc[3]))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 2*mlsp_coefs_pc[3]))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 3*mlsp_coefs_pc[3]))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 4*mlsp_coefs_pc[3]))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 5*mlsp_coefs_pc[3]))))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C8. Logistic function through the positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_pc <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_pc <- function(par, x, z, n, x50, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_pc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_pc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_pc = c(mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4]) # Starting parameters are the estimated parameters of the optimized function
avec_pc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_pc = numeric(100) # Prepare profile of likelihood values

pars_x50_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[3], mlsp_coefs_pc[4])
x50vec_pc = seq(-120, 100, length = 100)
x50prof_pc = numeric(100)

pars_c_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[4])
cvec_pc = seq(0.1, 50, length = 100)
cprof_pc = numeric(100)

pars_theta_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3])
thetavec_pc = seq(0.01, 10, length = 100)
thetaprof_pc = numeric(100)

# Optimize original function without fixed parameters
pars_2_pc <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_pc <- optim(par = pars_2_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n,
                           db1 = dat_6$year >= 2014, 
                           db2 = dat_6$year >= 2015, 
                           db3 = dat_6$year >= 2016, 
                           db4 = dat_6$year >= 2017, 
                           db5 = dat_6$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_pc <- list(rep(list(), times = 100))
x50_coefs_vec_pc <- list(rep(list(), times = 100))
c_coefs_vec_pc <- list(rep(list(), times = 100))
theta_coefs_vec_pc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_pc <- optim(par = pars_a_pc, fn = ci_a_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, a = avec_pc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_6$year >= 2014, 
                    db2 = dat_6$year >= 2015, 
                    db3 = dat_6$year >= 2016, 
                    db4 = dat_6$year >= 2017, 
                    db5 = dat_6$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_pc[[j]] = opt_a_pc$par # Store parameters for control
  # print(opt_a$convergence)
  aprof_pc[j] <- opt_a_pc$value # Set likelihood value for profile
  
  opt_x50_pc <- optim(par = pars_x50_pc, fn = ci_x50_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, x50 = x50vec_pc[j],
                      db1 = dat_6$year >= 2014, 
                      db2 = dat_6$year >= 2015, 
                      db3 = dat_6$year >= 2016, 
                      db4 = dat_6$year >= 2017, 
                      db5 = dat_6$year >= 2018, 
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_pc[j] <- opt_x50_pc$value
  x50_coefs_vec_pc[[j]] = opt_x50_pc$par
  
  opt_c_pc <- optim(par = pars_c_pc, fn = ci_c_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, c = cvec_pc[j],
                    db1 = dat_6$year >= 2014, 
                    db2 = dat_6$year >= 2015, 
                    db3 = dat_6$year >= 2016, 
                    db4 = dat_6$year >= 2017, 
                    db5 = dat_6$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_pc[j] <- opt_c_pc$value
  c_coefs_vec_pc[[j]] = opt_c_pc$par
  
  opt_theta_pc <- optim(par = pars_theta_pc, fn = ci_theta_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, theta = thetavec_pc[j], 
                        db1 = dat_6$year >= 2014, 
                        db2 = dat_6$year >= 2015, 
                        db3 = dat_6$year >= 2016, 
                        db4 = dat_6$year >= 2017, 
                        db5 = dat_6$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_pc[j] <- opt_theta_pc$value
  theta_coefs_vec_pc[[j]] = opt_theta_pc$par
}

# Set the 95% confidence limits
aprof_lower_pc <- aprof_pc[1:which.min(aprof_pc)] # Likelihood values below the best estimate
avec_lower_pc <- avec_pc[1:which.min(aprof_pc)] # Parameter values below the best estimate
aprof_higher_pc <- aprof_pc[which.min(aprof_pc):length(aprof_pc)] # Likelihood values above the best estimate
avec_higher_pc <- avec_pc[which.min(aprof_pc):length(aprof_pc)] # Parameter values above the best estimate
l_a_ci_pc <- approx(aprof_lower_pc, avec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_pc <- approx(aprof_higher_pc, avec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x50prof_lower_pc <- x50prof_pc[1:which.min(x50prof_pc)]
x50vec_lower_pc <- x50vec_pc[1:which.min(x50prof_pc)]
x50prof_higher_pc <- x50prof_pc[which.min(x50prof_pc):length(x50prof_pc)]
x50vec_higher_pc <- x50vec_pc[which.min(x50prof_pc):length(x50prof_pc)]
l_x50_ci_pc <- approx(x50prof_lower_pc, x50vec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_x50_ci_pc <- approx(x50prof_higher_pc, x50vec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

cprof_lower_pc <- cprof_pc[1:which.min(cprof_pc)]
cvec_lower_pc <- cvec_pc[1:which.min(cprof_pc)]
cprof_higher_pc <- cprof_pc[which.min(cprof_pc):length(cprof_pc)]
cvec_higher_pc <- cvec_pc[which.min(cprof_pc):length(cprof_pc)]
l_c_ci_pc <- approx(cprof_lower_pc, cvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_c_ci_pc <- approx(cprof_higher_pc, cvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

thetaprof_lower_pc <- thetaprof_pc[1:which.min(thetaprof_pc)]
thetavec_lower_pc <- thetavec_pc[1:which.min(thetaprof_pc)]
thetaprof_higher_pc <- thetaprof_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
thetavec_higher_pc <- thetavec_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
l_theta_ci_pc <- approx(thetaprof_lower_pc, thetavec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_theta_ci_pc <- approx(thetaprof_higher_pc, thetavec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_pc_log_yr <- c(opt_log_bbinom_pc$par[1], l_a_ci_pc$y, r_a_ci_pc$y) 
x50_ci_pc_log_yr <- c(opt_log_bbinom_pc$par[2], l_x50_ci_pc$y, r_x50_ci_pc$y)
c_ci_pc_log_yr <- c(opt_log_bbinom_pc$par[3], l_c_ci_pc$y, r_c_ci_pc$y)
theta_ci_pc_log_yr <- c(opt_log_bbinom_pc$par[4], l_theta_ci_pc$y, r_theta_ci_pc$y)

a_ci_pc_log_yr
x50_ci_pc_log_yr
c_ci_pc_log_yr
theta_ci_pc_log_yr


##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function
mlsp_fun_pc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_pc <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_pc <- optim(par = mlsp_par_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, 
                     db1 = dat_6$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_6$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_6$year >= 2016, 
                     db4 = dat_6$year >= 2017, 
                     db5 = dat_6$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_pc <- mlsp_opt_pc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_6[dat_6$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((0) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((0) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((1) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((1) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((2) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((2) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((3) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((3) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((4) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((4) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((5) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((5) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C9. CNE function through the positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_pc <- function(par, x, z, n, b, db1, db2, db3, db4, db5){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_pc <- function(par, x, z, n, x100, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_pc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_pc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_pc = c(mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_pc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_pc = numeric(100) # Prepare profile of likelihood values

pars_x100_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[3], mlsp_coefs_pc[4])
x100vec_pc = seq(-120, 100, length = 100)
x100prof_pc = numeric(100)

pars_c_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[4])
cvec_pc = seq(0.1, 50, length = 100)
cprof_pc = numeric(100)

pars_theta_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3])
thetavec_pc = seq(0.01, 10, length = 100)
thetaprof_pc = numeric(100)

# Optimize original function without fixed parameters
pars_2_pc <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_pc <- optim(par = pars_2_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n,
                           db1 = dat_6$year >= 2014, 
                           db2 = dat_6$year >= 2015, 
                           db3 = dat_6$year >= 2016, 
                           db4 = dat_6$year >= 2017, 
                           db5 = dat_6$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_pc <- list(rep(list(), times = 100))
x100_coefs_vec_pc <- list(rep(list(), times = 100))
c_coefs_vec_pc <- list(rep(list(), times = 100))
theta_coefs_vec_pc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_pc <- optim(par = pars_b_pc, fn = ci_b_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, b = bvec_pc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_6$year >= 2014, 
                    db2 = dat_6$year >= 2015, 
                    db3 = dat_6$year >= 2016, 
                    db4 = dat_6$year >= 2017, 
                    db5 = dat_6$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_pc[[j]] = opt_b_pc$par # Store parameters for control
  # print(opt_a$convergence)
  bprof_pc[j] <- opt_b_pc$value # Set likelihood value for profile
  
  opt_x100_pc <- optim(par = pars_x100_pc, fn = ci_x100_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, x100 = x100vec_pc[j],
                       db1 = dat_6$year >= 2014, 
                       db2 = dat_6$year >= 2015, 
                       db3 = dat_6$year >= 2016, 
                       db4 = dat_6$year >= 2017, 
                       db5 = dat_6$year >= 2018, 
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_pc[j] <- opt_x100_pc$value
  x100_coefs_vec_pc[[j]] = opt_x100_pc$par
  
  opt_c_pc <- optim(par = pars_c_pc, fn = ci_c_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, c = cvec_pc[j],
                    db1 = dat_6$year >= 2014, 
                    db2 = dat_6$year >= 2015, 
                    db3 = dat_6$year >= 2016, 
                    db4 = dat_6$year >= 2017, 
                    db5 = dat_6$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_pc[j] <- opt_c_pc$value
  c_coefs_vec_pc[[j]] = opt_c_pc$par
  
  opt_theta_pc <- optim(par = pars_theta_pc, fn = ci_theta_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, theta = thetavec_pc[j], 
                        db1 = dat_6$year >= 2014, 
                        db2 = dat_6$year >= 2015, 
                        db3 = dat_6$year >= 2016, 
                        db4 = dat_6$year >= 2017, 
                        db5 = dat_6$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_pc[j] <- opt_theta_pc$value
  theta_coefs_vec_pc[[j]] = opt_theta_pc$par
}

# Set the 95% confidence limits
bprof_lower_pc <- bprof_pc[1:which.min(bprof_pc)] # Likelihood values below the best estimate
bvec_lower_pc <- bvec_pc[1:which.min(bprof_pc)] # Parameter values below the best estimate
bprof_higher_pc <- bprof_pc[which.min(bprof_pc):length(bprof_pc)] # Likelihood values above the best estimate
bvec_higher_pc <- bvec_pc[which.min(bprof_pc):length(bprof_pc)] # Parameter values above the best estimate
l_b_ci_pc <- approx(bprof_lower_pc, bvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_pc <- approx(bprof_higher_pc, bvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_pc <- x100prof_pc[1:which.min(x100prof_pc)]
x100vec_lower_pc <- x100vec_pc[1:which.min(x100prof_pc)]
x100prof_higher_pc <- x100prof_pc[which.min(x100prof_pc):length(x100prof_pc)]
x100vec_higher_pc <- x100vec_pc[which.min(x100prof_pc):length(x100prof_pc)]
l_x100_ci_pc <- approx(x100prof_lower_pc, x100vec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_x100_ci_pc <- approx(x100prof_higher_pc, x100vec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

cprof_lower_pc <- cprof_pc[1:which.min(cprof_pc)]
cvec_lower_pc <- cvec_pc[1:which.min(cprof_pc)]
cprof_higher_pc <- cprof_pc[which.min(cprof_pc):length(cprof_pc)]
cvec_higher_pc <- cvec_pc[which.min(cprof_pc):length(cprof_pc)]
l_c_ci_pc <- approx(cprof_lower_pc, cvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_c_ci_pc <- approx(cprof_higher_pc, cvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

thetaprof_lower_pc <- thetaprof_pc[1:which.min(thetaprof_pc)]
thetavec_lower_pc <- thetavec_pc[1:which.min(thetaprof_pc)]
thetaprof_higher_pc <- thetaprof_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
thetavec_higher_pc <- thetavec_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
l_theta_ci_pc <- approx(thetaprof_lower_pc, thetavec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_theta_ci_pc <- approx(thetaprof_higher_pc, thetavec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_pc_cne_yr <- c(opt_log_bbinom_pc$par[1], l_b_ci_pc$y, r_b_ci_pc$y) 
x100_ci_pc_cne_yr <- c(opt_log_bbinom_pc$par[2], l_x100_ci_pc$y, r_x100_ci_pc$y)
c_ci_pc_cne_yr <- c(opt_log_bbinom_pc$par[3], l_c_ci_pc$y, r_c_ci_pc$y)
theta_ci_pc_cne_yr <- c(opt_log_bbinom_pc$par[4], l_theta_ci_pc$y, r_theta_ci_pc$y)

b_ci_pc_cne_yr
x100_ci_pc_cne_yr
c_ci_pc_cne_yr
theta_ci_pc_cne_yr


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function
mlsp_fun_pc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_pc <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_pc <- optim(par = mlsp_par_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, 
                     db1 = dat_6$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_6$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_6$year >= 2016, 
                     db4 = dat_6$year >= 2017, 
                     db5 = dat_6$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_pc <- mlsp_opt_pc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_6[dat_6$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(0)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(0))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(0))^mlsp_coefs_pc[2])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(1)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(1))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(1))^mlsp_coefs_pc[2])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(2)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(2))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(2))^mlsp_coefs_pc[2])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(3)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(3))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(3))^mlsp_coefs_pc[2])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(4)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(4))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(4))^mlsp_coefs_pc[2])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(5)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(5))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(5))^mlsp_coefs_pc[2])))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C10. Hill function through the positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_pc <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_pc <- function(par, x, z, n, pn, db1, db2, db3, db4, db5){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_pc <- function(par, x, z, n, h, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_pc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_pc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_pc = c(mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4], mlsp_coefs_pc[5]) # Starting parameters are the estimated parameters of the optimized function
avec_pc = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_pc = numeric(100) # Prepare profile of likelihood values

pars_pn_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[3], mlsp_coefs_pc[4], mlsp_coefs_pc[5])
pnvec_pc = seq(0.1, 10, length = 100)
pnprof_pc = numeric(100)

pars_h_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[4], mlsp_coefs_pc[5])
hvec_pc = seq(2, 100, length = 100)
hprof_pc = numeric(100)

pars_c_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[5])
cvec_pc = seq(0.1, 50, length = 100)
cprof_pc = numeric(100)

pars_theta_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4])
thetavec_pc = seq(0.01, 10, length = 100)
thetaprof_pc = numeric(100)

# Optimize original function without fixed parameters
pars_2_pc <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_pc <- optim(par = pars_2_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n,
                           db1 = dat_6$year >= 2014, 
                           db2 = dat_6$year >= 2015, 
                           db3 = dat_6$year >= 2016, 
                           db4 = dat_6$year >= 2017, 
                           db5 = dat_6$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_pc <- list(rep(list(), times = 100))
pn_coefs_vec_pc <- list(rep(list(), times = 100))
h_coefs_vec_pc <- list(rep(list(), times = 100))
c_coefs_vec_pc <- list(rep(list(), times = 100))
theta_coefs_vec_pc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_pc <- optim(par = pars_a_pc, fn = ci_a_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, a = avec_pc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_6$year >= 2014, 
                    db2 = dat_6$year >= 2015, 
                    db3 = dat_6$year >= 2016, 
                    db4 = dat_6$year >= 2017, 
                    db5 = dat_6$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_pc[[j]] = opt_a_pc$par # Store parameters for control
  aprof_pc[j] <- opt_a_pc$value # Set likelihood value for profile
  
  opt_pn_pc <- optim(par = pars_pn_pc, fn = ci_pn_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, pn = pnvec_pc[j],
                     db1 = dat_6$year >= 2014, 
                     db2 = dat_6$year >= 2015, 
                     db3 = dat_6$year >= 2016, 
                     db4 = dat_6$year >= 2017, 
                     db5 = dat_6$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_pc[j] <- opt_pn_pc$value
  pn_coefs_vec_pc[[j]] = opt_pn_pc$par
  
  opt_h_pc <- optim(par = pars_h_pc, fn = ci_h_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, h = hvec_pc[j],
                    db1 = dat_6$year >= 2014, 
                    db2 = dat_6$year >= 2015, 
                    db3 = dat_6$year >= 2016, 
                    db4 = dat_6$year >= 2017, 
                    db5 = dat_6$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_pc[j] <- opt_h_pc$value
  h_coefs_vec_pc[[j]] = opt_h_pc$par
  
  opt_c_pc <- optim(par = pars_c_pc, fn = ci_c_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, c = cvec_pc[j],
                    db1 = dat_6$year >= 2014, 
                    db2 = dat_6$year >= 2015, 
                    db3 = dat_6$year >= 2016, 
                    db4 = dat_6$year >= 2017, 
                    db5 = dat_6$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_pc[j] <- opt_c_pc$value
  c_coefs_vec_pc[[j]] = opt_c_pc$par
  
  opt_theta_pc <- optim(par = pars_theta_pc, fn = ci_theta_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, theta = thetavec_pc[j], 
                        db1 = dat_6$year >= 2014, 
                        db2 = dat_6$year >= 2015, 
                        db3 = dat_6$year >= 2016, 
                        db4 = dat_6$year >= 2017, 
                        db5 = dat_6$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_pc[j] <- opt_theta_pc$value
  theta_coefs_vec_pc[[j]] = opt_theta_pc$par
}

# Set the 95% confidence limits
aprof_lower_pc <- aprof_pc[1:which.min(aprof_pc)] # Likelihood values below the best estimate
avec_lower_pc <- avec_pc[1:which.min(aprof_pc)] # Parameter values below the best estimate
aprof_higher_pc <- aprof_pc[which.min(aprof_pc):length(aprof_pc)] # Likelihood values above the best estimate
avec_higher_pc <- avec_pc[which.min(aprof_pc):length(aprof_pc)] # Parameter values above the best estimate
l_a_ci_pc <- approx(aprof_lower_pc, avec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_pc <- approx(aprof_higher_pc, avec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_pc <- pnprof_pc[1:which.min(pnprof_pc)]
pnvec_lower_pc <- pnvec_pc[1:which.min(pnprof_pc)]
pnprof_higher_pc <- pnprof_pc[which.min(pnprof_pc):length(pnprof_pc)]
pnvec_higher_pc <- pnvec_pc[which.min(pnprof_pc):length(pnprof_pc)]
l_pn_ci_pc <- approx(pnprof_lower_pc, pnvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_pn_ci_pc <- approx(pnprof_higher_pc, pnvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

hprof_lower_pc <- hprof_pc[1:which.min(hprof_pc)]
hvec_lower_pc <- hvec_pc[1:which.min(hprof_pc)]
hprof_higher_pc <- hprof_pc[which.min(hprof_pc):length(hprof_pc)]
hvec_higher_pc <- hvec_pc[which.min(hprof_pc):length(hprof_pc)]
l_h_ci_pc <- approx(hprof_lower_pc, hvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_h_ci_pc <- approx(hprof_higher_pc, hvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

cprof_lower_pc <- cprof_pc[1:which.min(cprof_pc)]
cvec_lower_pc <- cvec_pc[1:which.min(cprof_pc)]
cprof_higher_pc <- cprof_pc[which.min(cprof_pc):length(cprof_pc)]
cvec_higher_pc <- cvec_pc[which.min(cprof_pc):length(cprof_pc)]
l_c_ci_pc <- approx(cprof_lower_pc, cvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_c_ci_pc <- approx(cprof_higher_pc, cvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

thetaprof_lower_pc <- thetaprof_pc[1:which.min(thetaprof_pc)]
thetavec_lower_pc <- thetavec_pc[1:which.min(thetaprof_pc)]
thetaprof_higher_pc <- thetaprof_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
thetavec_higher_pc <- thetavec_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
l_theta_ci_pc <- approx(thetaprof_lower_pc, thetavec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_theta_ci_pc <- approx(thetaprof_higher_pc, thetavec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_pc_hill_yr <- c(opt_log_bbinom_pc$par[1], l_a_ci_pc$y, r_a_ci_pc$y) 
pn_ci_pc_hill_yr <- c(opt_log_bbinom_pc$par[2], l_pn_ci_pc$y, r_pn_ci_pc$y)
h_ci_pc_hill_yr <- c(opt_log_bbinom_pc$par[3], l_h_ci_pc$y, r_h_ci_pc$y)
c_ci_pc_hill_yr <- c(opt_log_bbinom_pc$par[4], l_c_ci_pc$y, r_c_ci_pc$y)
theta_ci_pc_hill_yr <- c(opt_log_bbinom_pc$par[5], l_theta_ci_pc$y, r_theta_ci_pc$y)

a_ci_pc_hill_yr
pn_ci_pc_hill_yr
h_ci_pc_hill_yr
c_ci_pc_hill_yr
theta_ci_pc_hill_yr


##############################
# Negative clairvoyance data #
##############################

#
# Data morphing and distance circle creation
#

dat_7 <- dat

# Make data clairvoyant
for(i in length(years):2){ # Work backwards, from 2018 to 2014
  co <- which(dat_7[[i]]$result == 0) # Which of the data have a result of 0
  dat_7[[i-1]] <- rbind(dat_7[[i-1]], dat_7[[i]][co,]) # Combine the previous year with the negative results of the later year
  co_2 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & dat_7[[i-1]]$season == (2012 + i)) # Check for multiple samples at the same location
  if(length(co_2) >= 1){ # If there are multiple samples at one location
    for(j in 1:length(co_2)){ # For every duplicate sample
      if(dat_7[[i-1]]$result[co_2[j]] == 1){ # If the duplicate from the year before has result 1, the result at this location should be a 1
        co_3 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & dat_7[[i-1]]$season == (2012 + i))
        dat_7[[i-1]] <- dat_7[[i-1]][-co_3[j],]
      }
      if(dat_7[[i-1]]$result[co_2[j]] == 0){ # If the result is a 0, it will stay a 0.
        dat_7[[i-1]] <- dat_7[[i-1]][-co_2[j],]
      }
    }
  } # Note that this should never be the case, since a 1 would have turned into a 0. This part of code is a safeguard to be used in future data or simulated data.
}

# Set distance circles
dc <- 1:(max_dist) # Number of distance circles is the maximum distance a measurement is done.
pos <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for number of positives in each dc
n <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for total number of measurements in each dc

dat_8 <- list() # Prepare data list

for(i in 1:length(years)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j] <- sum((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j] <- length((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_8[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
}

# Combine years and create year column for year 
dat_9 <- list()
for(i in 1:length(years)){
  dat_8[[i]] <- dat_8[[i]] %>% 
    mutate(year = (2012 + i))
  dat_9 <- rbind(dat_9, dat_8[[i]])
}

# Plot the data
ggplot() +
  geom_point(data = dat_9[dat_9$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  labs(x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C11. Data points of the negative clairvoyance data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function
mlsp_fun_nc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- (1/(1 + exp(par[1] * (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
  return(nll)
}

# Setup parameter starting values
mlsp_par_nc <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_nc <- optim(par = mlsp_par_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, 
                     db1 = dat_9$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_9$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_9$year >= 2016, 
                     db4 = dat_9$year >= 2017, 
                     db5 = dat_9$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_nc <- mlsp_opt_nc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_9[dat_9$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 0*mlsp_coefs_nc[3]))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 1*mlsp_coefs_nc[3]))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 2*mlsp_coefs_nc[3]))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 3*mlsp_coefs_nc[3]))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 4*mlsp_coefs_nc[3]))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 5*mlsp_coefs_nc[3]))))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C12. Logistic function through the negative clairvoyance data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_nc <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_nc <- function(par, x, z, n, x50, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_nc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_nc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_nc = c(mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4]) # Starting parameters are the estimated parameters of the optimized function
avec_nc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_nc = numeric(100) # Prepare profile of likelihood values

pars_x50_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[3], mlsp_coefs_nc[4])
x50vec_nc = seq(-120, 100, length = 100)
x50prof_nc = numeric(100)

pars_c_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[4])
cvec_nc = seq(0.1, 50, length = 100)
cprof_nc = numeric(100)

pars_theta_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3])
thetavec_nc = seq(0.01, 10, length = 100)
thetaprof_nc = numeric(100)

# Optimize original function without fixed parameters
pars_2_nc <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_nc <- optim(par = pars_2_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n,
                           db1 = dat_9$year >= 2014, 
                           db2 = dat_9$year >= 2015, 
                           db3 = dat_9$year >= 2016, 
                           db4 = dat_9$year >= 2017, 
                           db5 = dat_9$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_nc <- list(rep(list(), times = 100))
x50_coefs_vec_nc <- list(rep(list(), times = 100))
c_coefs_vec_nc <- list(rep(list(), times = 100))
theta_coefs_vec_nc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_nc <- optim(par = pars_a_nc, fn = ci_a_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, a = avec_nc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_9$year >= 2014, 
                    db2 = dat_9$year >= 2015, 
                    db3 = dat_9$year >= 2016, 
                    db4 = dat_9$year >= 2017, 
                    db5 = dat_9$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_nc[[j]] = opt_a_nc$par # Store parameters for control
  # print(opt_a$convergence)
  aprof_nc[j] <- opt_a_nc$value # Set likelihood value for profile
  
  opt_x50_nc <- optim(par = pars_x50_nc, fn = ci_x50_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, x50 = x50vec_nc[j],
                      db1 = dat_9$year >= 2014, 
                      db2 = dat_9$year >= 2015, 
                      db3 = dat_9$year >= 2016, 
                      db4 = dat_9$year >= 2017, 
                      db5 = dat_9$year >= 2018, 
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_nc[j] <- opt_x50_nc$value
  x50_coefs_vec_nc[[j]] = opt_x50_nc$par
  
  opt_c_nc <- optim(par = pars_c_nc, fn = ci_c_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, c = cvec_nc[j],
                    db1 = dat_9$year >= 2014, 
                    db2 = dat_9$year >= 2015, 
                    db3 = dat_9$year >= 2016, 
                    db4 = dat_9$year >= 2017, 
                    db5 = dat_9$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_nc[j] <- opt_c_nc$value
  c_coefs_vec_nc[[j]] = opt_c_nc$par
  
  opt_theta_nc <- optim(par = pars_theta_nc, fn = ci_theta_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, theta = thetavec_nc[j], 
                        db1 = dat_9$year >= 2014, 
                        db2 = dat_9$year >= 2015, 
                        db3 = dat_9$year >= 2016, 
                        db4 = dat_9$year >= 2017, 
                        db5 = dat_9$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_nc[j] <- opt_theta_nc$value
  theta_coefs_vec_nc[[j]] = opt_theta_nc$par
}

# Set the 95% confidence limits
aprof_lower_nc <- aprof_nc[1:which.min(aprof_nc)] # Likelihood values below the best estimate
avec_lower_nc <- avec_nc[1:which.min(aprof_nc)] # Parameter values below the best estimate
aprof_higher_nc <- aprof_nc[which.min(aprof_nc):length(aprof_nc)] # Likelihood values above the best estimate
avec_higher_nc <- avec_nc[which.min(aprof_nc):length(aprof_nc)] # Parameter values above the best estimate
l_a_ci_nc <- approx(aprof_lower_nc, avec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_nc <- approx(aprof_higher_nc, avec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x50prof_lower_nc <- x50prof_nc[1:which.min(x50prof_nc)]
x50vec_lower_nc <- x50vec_nc[1:which.min(x50prof_nc)]
x50prof_higher_nc <- x50prof_nc[which.min(x50prof_nc):length(x50prof_nc)]
x50vec_higher_nc <- x50vec_nc[which.min(x50prof_nc):length(x50prof_nc)]
l_x50_ci_nc <- approx(x50prof_lower_nc, x50vec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_x50_ci_nc <- approx(x50prof_higher_nc, x50vec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

cprof_lower_nc <- cprof_nc[1:which.min(cprof_nc)]
cvec_lower_nc <- cvec_nc[1:which.min(cprof_nc)]
cprof_higher_nc <- cprof_nc[which.min(cprof_nc):length(cprof_nc)]
cvec_higher_nc <- cvec_nc[which.min(cprof_nc):length(cprof_nc)]
l_c_ci_nc <- approx(cprof_lower_nc, cvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_c_ci_nc <- approx(cprof_higher_nc, cvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

thetaprof_lower_nc <- thetaprof_nc[1:which.min(thetaprof_nc)]
thetavec_lower_nc <- thetavec_nc[1:which.min(thetaprof_nc)]
thetaprof_higher_nc <- thetaprof_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
thetavec_higher_nc <- thetavec_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
l_theta_ci_nc <- approx(thetaprof_lower_nc, thetavec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_theta_ci_nc <- approx(thetaprof_higher_nc, thetavec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_nc_log_yr <- c(opt_log_bbinom_nc$par[1], l_a_ci_nc$y, r_a_ci_nc$y) 
x50_ci_nc_log_yr <- c(opt_log_bbinom_nc$par[2], l_x50_ci_nc$y, r_x50_ci_nc$y)
c_ci_nc_log_yr <- c(opt_log_bbinom_nc$par[3], l_c_ci_nc$y, r_c_ci_nc$y)
theta_ci_nc_log_yr <- c(opt_log_bbinom_nc$par[4], l_theta_ci_nc$y, r_theta_ci_nc$y)

a_ci_nc_log_yr
x50_ci_nc_log_yr
c_ci_nc_log_yr
theta_ci_nc_log_yr


##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function
mlsp_fun_nc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_nc <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_nc <- optim(par = mlsp_par_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, 
                     db1 = dat_9$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_9$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_9$year >= 2016, 
                     db4 = dat_9$year >= 2017, 
                     db5 = dat_9$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_nc <- mlsp_opt_nc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_9[dat_9$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((0) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((0) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((1) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((1) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((2) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((2) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((3) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((3) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((4) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((4) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((5) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((5) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C13. CNE function through the negative clairvoyance data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_nc <- function(par, x, z, n, b, db1, db2, db3, db4, db5){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_nc <- function(par, x, z, n, x100, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_nc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_nc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_nc = c(mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_nc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_nc = numeric(100) # Prepare profile of likelihood values

pars_x100_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[3], mlsp_coefs_nc[4])
x100vec_nc = seq(-120, 100, length = 100)
x100prof_nc = numeric(100)

pars_c_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[4])
cvec_nc = seq(0.1, 50, length = 100)
cprof_nc = numeric(100)

pars_theta_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3])
thetavec_nc = seq(0.01, 10, length = 100)
thetaprof_nc = numeric(100)

# Optimize original function without fixed parameters
pars_2_nc <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_nc <- optim(par = pars_2_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n,
                           db1 = dat_9$year >= 2014, 
                           db2 = dat_9$year >= 2015, 
                           db3 = dat_9$year >= 2016, 
                           db4 = dat_9$year >= 2017, 
                           db5 = dat_9$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_nc <- list(rep(list(), times = 100))
x100_coefs_vec_nc <- list(rep(list(), times = 100))
c_coefs_vec_nc <- list(rep(list(), times = 100))
theta_coefs_vec_nc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_nc <- optim(par = pars_b_nc, fn = ci_b_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, b = bvec_nc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_9$year >= 2014, 
                    db2 = dat_9$year >= 2015, 
                    db3 = dat_9$year >= 2016, 
                    db4 = dat_9$year >= 2017, 
                    db5 = dat_9$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_nc[[j]] = opt_b_nc$par # Store parameters for control
  bprof_nc[j] <- opt_b_nc$value # Set likelihood value for profile
  
  opt_x100_nc <- optim(par = pars_x100_nc, fn = ci_x100_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, x100 = x100vec_nc[j],
                       db1 = dat_9$year >= 2014, 
                       db2 = dat_9$year >= 2015, 
                       db3 = dat_9$year >= 2016, 
                       db4 = dat_9$year >= 2017, 
                       db5 = dat_9$year >= 2018, 
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_nc[j] <- opt_x100_nc$value
  x100_coefs_vec_nc[[j]] = opt_x100_nc$par
  
  opt_c_nc <- optim(par = pars_c_nc, fn = ci_c_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, c = cvec_nc[j],
                    db1 = dat_9$year >= 2014, 
                    db2 = dat_9$year >= 2015, 
                    db3 = dat_9$year >= 2016, 
                    db4 = dat_9$year >= 2017, 
                    db5 = dat_9$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_nc[j] <- opt_c_nc$value
  c_coefs_vec_nc[[j]] = opt_c_nc$par
  
  opt_theta_nc <- optim(par = pars_theta_nc, fn = ci_theta_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, theta = thetavec_nc[j], 
                        db1 = dat_9$year >= 2014, 
                        db2 = dat_9$year >= 2015, 
                        db3 = dat_9$year >= 2016, 
                        db4 = dat_9$year >= 2017, 
                        db5 = dat_9$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_nc[j] <- opt_theta_nc$value
  theta_coefs_vec_nc[[j]] = opt_theta_nc$par
}

# Set the 95% confidence limits
bprof_lower_nc <- bprof_nc[1:which.min(bprof_nc)] # Likelihood values below the best estimate
bvec_lower_nc <- bvec_nc[1:which.min(bprof_nc)] # Parameter values below the best estimate
bprof_higher_nc <- bprof_nc[which.min(bprof_nc):length(bprof_nc)] # Likelihood values above the best estimate
bvec_higher_nc <- bvec_nc[which.min(bprof_nc):length(bprof_nc)] # Parameter values above the best estimate
l_b_ci_nc <- approx(bprof_lower_nc, bvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_nc <- approx(bprof_higher_nc, bvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_nc <- x100prof_nc[1:which.min(x100prof_nc)]
x100vec_lower_nc <- x100vec_nc[1:which.min(x100prof_nc)]
x100prof_higher_nc <- x100prof_nc[which.min(x100prof_nc):length(x100prof_nc)]
x100vec_higher_nc <- x100vec_nc[which.min(x100prof_nc):length(x100prof_nc)]
l_x100_ci_nc <- approx(x100prof_lower_nc, x100vec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_x100_ci_nc <- approx(x100prof_higher_nc, x100vec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

cprof_lower_nc <- cprof_nc[1:which.min(cprof_nc)]
cvec_lower_nc <- cvec_nc[1:which.min(cprof_nc)]
cprof_higher_nc <- cprof_nc[which.min(cprof_nc):length(cprof_nc)]
cvec_higher_nc <- cvec_nc[which.min(cprof_nc):length(cprof_nc)]
l_c_ci_nc <- approx(cprof_lower_nc, cvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_c_ci_nc <- approx(cprof_higher_nc, cvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

thetaprof_lower_nc <- thetaprof_nc[1:which.min(thetaprof_nc)]
thetavec_lower_nc <- thetavec_nc[1:which.min(thetaprof_nc)]
thetaprof_higher_nc <- thetaprof_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
thetavec_higher_nc <- thetavec_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
l_theta_ci_nc <- approx(thetaprof_lower_nc, thetavec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_theta_ci_nc <- approx(thetaprof_higher_nc, thetavec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_nc_cne_yr <- c(opt_log_bbinom_nc$par[1], l_b_ci_nc$y, r_b_ci_nc$y) 
x100_ci_nc_cne_yr <- c(opt_log_bbinom_nc$par[2], l_x100_ci_nc$y, r_x100_ci_nc$y)
c_ci_nc_cne_yr <- c(opt_log_bbinom_nc$par[3], l_c_ci_nc$y, r_c_ci_nc$y)
theta_ci_nc_cne_yr <- c(opt_log_bbinom_nc$par[4], l_theta_ci_nc$y, r_theta_ci_nc$y)

b_ci_nc_cne_yr
x100_ci_nc_cne_yr
c_ci_nc_cne_yr
theta_ci_nc_cne_yr


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function
mlsp_fun_nc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_nc <- c(a = 0.1, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_nc <- optim(par = mlsp_par_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, 
                     db1 = dat_9$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_9$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_9$year >= 2016, 
                     db4 = dat_9$year >= 2017, 
                     db5 = dat_9$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_nc <- mlsp_opt_nc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_9[dat_9$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(0)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(0))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(0))^mlsp_coefs_nc[2])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(1)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(1))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(1))^mlsp_coefs_nc[2])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(2)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(2))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(2))^mlsp_coefs_nc[2])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(3)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(3))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(3))^mlsp_coefs_nc[2])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(4)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(4))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(4))^mlsp_coefs_nc[2])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(5)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(5))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(5))^mlsp_coefs_nc[2])))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C14. Hill function through the negative clairvoyance data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_nc <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_nc <- function(par, x, z, n, pn, db1, db2, db3, db4, db5){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_nc <- function(par, x, z, n, h, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_nc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_nc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_nc = c(mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4], mlsp_coefs_nc[5]) # Starting parameters are the estimated parameters of the optimized function
avec_nc = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_nc = numeric(100) # Prepare profile of likelihood values

pars_pn_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[3], mlsp_coefs_nc[4], mlsp_coefs_nc[5])
pnvec_nc = seq(0.1, 10, length = 100)
pnprof_nc = numeric(100)

pars_h_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[4], mlsp_coefs_nc[5])
hvec_nc = seq(2, 100, length = 100)
hprof_nc = numeric(100)

pars_c_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[5])
cvec_nc = seq(0.1, 50, length = 100)
cprof_nc = numeric(100)

pars_theta_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4])
thetavec_nc = seq(0.01, 10, length = 100)
thetaprof_nc = numeric(100)

# Optimize original function without fixed parameters
pars_2_nc <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_nc <- optim(par = pars_2_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n,
                           db1 = dat_9$year >= 2014, 
                           db2 = dat_9$year >= 2015, 
                           db3 = dat_9$year >= 2016, 
                           db4 = dat_9$year >= 2017, 
                           db5 = dat_9$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_nc <- list(rep(list(), times = 100))
pn_coefs_vec_nc <- list(rep(list(), times = 100))
h_coefs_vec_nc <- list(rep(list(), times = 100))
c_coefs_vec_nc <- list(rep(list(), times = 100))
theta_coefs_vec_nc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_nc <- optim(par = pars_a_nc, fn = ci_a_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, a = avec_nc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_9$year >= 2014, 
                    db2 = dat_9$year >= 2015, 
                    db3 = dat_9$year >= 2016, 
                    db4 = dat_9$year >= 2017, 
                    db5 = dat_9$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_nc[[j]] = opt_a_nc$par # Store parameters for control
  aprof_nc[j] <- opt_a_nc$value # Set likelihood value for profile
  
  opt_pn_nc <- optim(par = pars_pn_nc, fn = ci_pn_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, pn = pnvec_nc[j],
                     db1 = dat_9$year >= 2014, 
                     db2 = dat_9$year >= 2015, 
                     db3 = dat_9$year >= 2016, 
                     db4 = dat_9$year >= 2017, 
                     db5 = dat_9$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_nc[j] <- opt_pn_nc$value
  pn_coefs_vec_nc[[j]] = opt_pn_nc$par
  
  opt_h_nc <- optim(par = pars_h_nc, fn = ci_h_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, h = hvec_nc[j],
                    db1 = dat_9$year >= 2014, 
                    db2 = dat_9$year >= 2015, 
                    db3 = dat_9$year >= 2016, 
                    db4 = dat_9$year >= 2017, 
                    db5 = dat_9$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_nc[j] <- opt_h_nc$value
  h_coefs_vec_nc[[j]] = opt_h_nc$par
  
  opt_c_nc <- optim(par = pars_c_nc, fn = ci_c_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, c = cvec_nc[j],
                    db1 = dat_9$year >= 2014, 
                    db2 = dat_9$year >= 2015, 
                    db3 = dat_9$year >= 2016, 
                    db4 = dat_9$year >= 2017, 
                    db5 = dat_9$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_nc[j] <- opt_c_nc$value
  c_coefs_vec_nc[[j]] = opt_c_nc$par
  
  opt_theta_nc <- optim(par = pars_theta_nc, fn = ci_theta_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, theta = thetavec_nc[j], 
                        db1 = dat_9$year >= 2014, 
                        db2 = dat_9$year >= 2015, 
                        db3 = dat_9$year >= 2016, 
                        db4 = dat_9$year >= 2017, 
                        db5 = dat_9$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_nc[j] <- opt_theta_nc$value
  theta_coefs_vec_nc[[j]] = opt_theta_nc$par
}

# Set the 95% confidence limits
aprof_lower_nc <- aprof_nc[1:which.min(aprof_nc)] # Likelihood values below the best estimate
avec_lower_nc <- avec_nc[1:which.min(aprof_nc)] # Parameter values below the best estimate
aprof_higher_nc <- aprof_nc[which.min(aprof_nc):length(aprof_nc)] # Likelihood values above the best estimate
avec_higher_nc <- avec_nc[which.min(aprof_nc):length(aprof_nc)] # Parameter values above the best estimate
l_a_ci_nc <- approx(aprof_lower_nc, avec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_nc <- approx(aprof_higher_nc, avec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_nc <- pnprof_nc[1:which.min(pnprof_nc)]
pnvec_lower_nc <- pnvec_nc[1:which.min(pnprof_nc)]
pnprof_higher_nc <- pnprof_nc[which.min(pnprof_nc):length(pnprof_nc)]
pnvec_higher_nc <- pnvec_nc[which.min(pnprof_nc):length(pnprof_nc)]
l_pn_ci_nc <- approx(pnprof_lower_nc, pnvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_pn_ci_nc <- approx(pnprof_higher_nc, pnvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

hprof_lower_nc <- hprof_nc[1:which.min(hprof_nc)]
hvec_lower_nc <- hvec_nc[1:which.min(hprof_nc)]
hprof_higher_nc <- hprof_nc[which.min(hprof_nc):length(hprof_nc)]
hvec_higher_nc <- hvec_nc[which.min(hprof_nc):length(hprof_nc)]
l_h_ci_nc <- approx(hprof_lower_nc, hvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_h_ci_nc <- approx(hprof_higher_nc, hvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

cprof_lower_nc <- cprof_nc[1:which.min(cprof_nc)]
cvec_lower_nc <- cvec_nc[1:which.min(cprof_nc)]
cprof_higher_nc <- cprof_nc[which.min(cprof_nc):length(cprof_nc)]
cvec_higher_nc <- cvec_nc[which.min(cprof_nc):length(cprof_nc)]
l_c_ci_nc <- approx(cprof_lower_nc, cvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_c_ci_nc <- approx(cprof_higher_nc, cvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

thetaprof_lower_nc <- thetaprof_nc[1:which.min(thetaprof_nc)]
thetavec_lower_nc <- thetavec_nc[1:which.min(thetaprof_nc)]
thetaprof_higher_nc <- thetaprof_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
thetavec_higher_nc <- thetavec_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
l_theta_ci_nc <- approx(thetaprof_lower_nc, thetavec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_theta_ci_nc <- approx(thetaprof_higher_nc, thetavec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_nc_hill_yr <- c(opt_log_bbinom_nc$par[1], l_a_ci_nc$y, r_a_ci_nc$y) 
pn_ci_nc_hill_yr <- c(opt_log_bbinom_nc$par[2], l_pn_ci_nc$y, r_pn_ci_nc$y)
h_ci_nc_hill_yr <- c(opt_log_bbinom_nc$par[3], l_h_ci_nc$y, r_h_ci_nc$y)
c_ci_nc_hill_yr <- c(opt_log_bbinom_nc$par[4], l_c_ci_nc$y, r_c_ci_nc$y)
theta_ci_nc_hill_yr <- c(opt_log_bbinom_nc$par[5], l_theta_ci_nc$y, r_theta_ci_nc$y)

a_ci_nc_hill_yr
pn_ci_nc_hill_yr
h_ci_nc_hill_yr
c_ci_nc_hill_yr
theta_ci_nc_hill_yr


####################################################
# Negative clairvoyance & positive carry-over data #
####################################################

#
# Data morphing and distance circle creation
#

dat_10 <- dat

# Make data clairvoyant
for(i in length(years):2){ # Work backwards, from 2018 to 2014
  co <- which(dat_10[[i]]$result == 0) # Which of the data have a result of 0
  dat_10[[i-1]] <- rbind(dat_10[[i-1]], dat_10[[i]][co,]) # Combine the previous year with the negative results of the later year
  co_2 <- which(dat_10[[i-1]]$lon == dat_10[[i]]$lon[co] & dat_10[[i-1]]$lat == dat_10[[i]]$lat[co] & dat_10[[i-1]]$season == (2012 + i)) # Check for multiple samples at the same location
  if(length(co_2) >= 1){ # If there are multiple samples at one location
    for(j in 1:length(co_2)){ # For every duplicate sample
      if(dat_10[[i-1]]$result[co_2[j]] == 1){ # If the duplicate from the year before has result 1, the result at this location should be a 1
        co_3 <- which(dat_10[[i-1]]$lon == dat_10[[i]]$lon[co] & dat_10[[i-1]]$lat == dat_10[[i]]$lat[co] & dat_10[[i-1]]$season == (2012 + i))
        dat_10[[i-1]] <- dat_10[[i-1]][-co_3[j],]
      }
      if(dat_10[[i-1]]$result[co_2[j]] == 0){ # If the result is a 0, it will stay a 0.
        dat_10[[i-1]] <- dat_10[[i-1]][-co_2[j],]
      }
    }
  } # Note that this should never be the case, since a 1 would have turned into a 0. This part of code is a safeguard to be used in future data or simulated data.
}

# make positives carry over
for(i in 2:length(years)){ # For-loop through every year from 2014
  co <- which(dat_10[[i-1]]$result == 1) # Which lines of the previous year have a result of 1
  dat_10[[i]] <- rbind(dat_10[[i]], dat_10[[i-1]][co,]) # Include the lines of the previous year with a result of 1 with the current year
  co_2 <- which(dat_10[[i]]$lon == dat_10[[i-1]]$lon[co] & dat_10[[i]]$lat == dat_10[[i-1]]$lat[co] & dat_10[[i]]$season == (2012 + i)) # Which lines of this year have the same coordinates of the positives last year
  # print(co_2)
  if(length(co_2) >= 1){ # If this is not 0, then there has been a measurement at these coordinates this year as well as last year, or multiple times in a single year (with these data, only multiple measurements in 2013)
    dat_10[[i]] <- dat_10[[i]][-co_2,] # Remove all the duplicates. If the duplicate measurement from this year is a 0, it will still be a 1, because I assume that once a disease is there, it will not leave. I assume that 0 to be a false negative.
  }
}

# Set distance circles
dc <- 1:(max_dist) # Number of distance circles is the maximum distance a measurement is done.
pos <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for number of positives in each dc
n <- rep(list(rep(NA, times = length(dc))), times = length(years)) # Prepare list for total number of measurements in each dc

dat_11 <- list() # Prepare data list

for(i in 1:length(years)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j] <- sum((dat_10[[i]]$result[dat_10[[i]]$dist >= j & (dat_10[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j] <- length((dat_10[[i]]$result[dat_10[[i]]$dist >= j & (dat_10[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_11[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
}

# Combine years and create year column for year 
dat_12 <- list()
for(i in 1:length(years)){
  dat_11[[i]] <- dat_11[[i]] %>% 
    mutate(year = (2012 + i))
  dat_12 <- rbind(dat_12, dat_11[[i]])
}

# Plot dat_12
ggplot() +
  geom_point(data = dat_12[dat_12$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  labs(
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C15. Data points of the negative clairvoyance & positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))

##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function
mlsp_fun_ncpc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- (1/(1 + exp(par[1] * (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3]))))) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
  return(nll)
}

# Setup parameter starting values
mlsp_par_ncpc <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_ncpc <- optim(par = mlsp_par_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, 
                     db1 = dat_12$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_12$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_12$year >= 2016, 
                     db4 = dat_12$year >= 2017, 
                     db5 = dat_12$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_ncpc <- mlsp_opt_ncpc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_12[dat_12$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 0*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 1*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 2*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 3*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 4*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 5*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C16. Logistic function through the negative clairvoyance & positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_ncpc <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_ncpc <- function(par, x, z, n, x50, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_ncpc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_ncpc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_ncpc = c(mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4]) # Starting parameters are the estimated parameters of the optimized function
avec_ncpc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_ncpc = numeric(100) # Prepare profile of likelihood values

pars_x50_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4])
x50vec_ncpc = seq(-120, 100, length = 100)
x50prof_ncpc = numeric(100)

pars_c_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[4])
cvec_ncpc = seq(0.1, 50, length = 100)
cprof_ncpc = numeric(100)

pars_theta_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3])
thetavec_ncpc = seq(0.01, 10, length = 100)
thetaprof_ncpc = numeric(100)

# Optimize original function without fixed parameters
pars_2_ncpc <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_ncpc <- optim(par = pars_2_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n,
                           db1 = dat_12$year >= 2014, 
                           db2 = dat_12$year >= 2015, 
                           db3 = dat_12$year >= 2016, 
                           db4 = dat_12$year >= 2017, 
                           db5 = dat_12$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_ncpc <- list(rep(list(), times = 100))
x50_coefs_vec_ncpc <- list(rep(list(), times = 100))
c_coefs_vec_ncpc <- list(rep(list(), times = 100))
theta_coefs_vec_ncpc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_ncpc <- optim(par = pars_a_ncpc, fn = ci_a_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, a = avec_ncpc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_12$year >= 2014, 
                    db2 = dat_12$year >= 2015, 
                    db3 = dat_12$year >= 2016, 
                    db4 = dat_12$year >= 2017, 
                    db5 = dat_12$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_ncpc[[j]] = opt_a_ncpc$par # Store parameters for control
  # print(opt_a$convergence)
  aprof_ncpc[j] <- opt_a_ncpc$value # Set likelihood value for profile
  
  opt_x50_ncpc <- optim(par = pars_x50_ncpc, fn = ci_x50_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, x50 = x50vec_ncpc[j],
                      db1 = dat_12$year >= 2014, 
                      db2 = dat_12$year >= 2015, 
                      db3 = dat_12$year >= 2016, 
                      db4 = dat_12$year >= 2017, 
                      db5 = dat_12$year >= 2018, 
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_ncpc[j] <- opt_x50_ncpc$value
  x50_coefs_vec_ncpc[[j]] = opt_x50_ncpc$par
  
  opt_c_ncpc <- optim(par = pars_c_ncpc, fn = ci_c_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, c = cvec_ncpc[j],
                    db1 = dat_12$year >= 2014, 
                    db2 = dat_12$year >= 2015, 
                    db3 = dat_12$year >= 2016, 
                    db4 = dat_12$year >= 2017, 
                    db5 = dat_12$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_ncpc[j] <- opt_c_ncpc$value
  c_coefs_vec_ncpc[[j]] = opt_c_ncpc$par
  
  opt_theta_ncpc <- optim(par = pars_theta_ncpc, fn = ci_theta_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, theta = thetavec_ncpc[j], 
                        db1 = dat_12$year >= 2014, 
                        db2 = dat_12$year >= 2015, 
                        db3 = dat_12$year >= 2016, 
                        db4 = dat_12$year >= 2017, 
                        db5 = dat_12$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_ncpc[j] <- opt_theta_ncpc$value
  theta_coefs_vec_ncpc[[j]] = opt_theta_ncpc$par
}

# Set the 95% confidence limits
aprof_lower_ncpc <- aprof_ncpc[1:which.min(aprof_ncpc)] # Likelihood values below the best estimate
avec_lower_ncpc <- avec_ncpc[1:which.min(aprof_ncpc)] # Parameter values below the best estimate
aprof_higher_ncpc <- aprof_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Likelihood values above the best estimate
avec_higher_ncpc <- avec_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Parameter values above the best estimate
l_a_ci_ncpc <- approx(aprof_lower_ncpc, avec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_ncpc <- approx(aprof_higher_ncpc, avec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x50prof_lower_ncpc <- x50prof_ncpc[1:which.min(x50prof_ncpc)]
x50vec_lower_ncpc <- x50vec_ncpc[1:which.min(x50prof_ncpc)]
x50prof_higher_ncpc <- x50prof_ncpc[which.min(x50prof_ncpc):length(x50prof_ncpc)]
x50vec_higher_ncpc <- x50vec_ncpc[which.min(x50prof_ncpc):length(x50prof_ncpc)]
l_x50_ci_ncpc <- approx(x50prof_lower_ncpc, x50vec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_x50_ci_ncpc <- approx(x50prof_higher_ncpc, x50vec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

cprof_lower_ncpc <- cprof_ncpc[1:which.min(cprof_ncpc)]
cvec_lower_ncpc <- cvec_ncpc[1:which.min(cprof_ncpc)]
cprof_higher_ncpc <- cprof_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
cvec_higher_ncpc <- cvec_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
l_c_ci_ncpc <- approx(cprof_lower_ncpc, cvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_c_ci_ncpc <- approx(cprof_higher_ncpc, cvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

thetaprof_lower_ncpc <- thetaprof_ncpc[1:which.min(thetaprof_ncpc)]
thetavec_lower_ncpc <- thetavec_ncpc[1:which.min(thetaprof_ncpc)]
thetaprof_higher_ncpc <- thetaprof_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
thetavec_higher_ncpc <- thetavec_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
l_theta_ci_ncpc <- approx(thetaprof_lower_ncpc, thetavec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_theta_ci_ncpc <- approx(thetaprof_higher_ncpc, thetavec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_ncpc_log_yr <- c(opt_log_bbinom_ncpc$par[1], l_a_ci_ncpc$y, r_a_ci_ncpc$y) 
x50_ci_ncpc_log_yr <- c(opt_log_bbinom_ncpc$par[2], l_x50_ci_ncpc$y, r_x50_ci_ncpc$y)
c_ci_ncpc_log_yr <- c(opt_log_bbinom_ncpc$par[3], l_c_ci_ncpc$y, r_c_ci_ncpc$y)
theta_ci_ncpc_log_yr <- c(opt_log_bbinom_ncpc$par[4], l_theta_ci_ncpc$y, r_theta_ci_ncpc$y)

a_ci_ncpc_log_yr
x50_ci_ncpc_log_yr
c_ci_ncpc_log_yr
theta_ci_ncpc_log_yr


##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function
mlsp_fun_ncpc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_ncpc <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_ncpc <- optim(par = mlsp_par_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, 
                     db1 = dat_12$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_12$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_12$year >= 2016, 
                     db4 = dat_12$year >= 2017, 
                     db5 = dat_12$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_ncpc <- mlsp_opt_ncpc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_12[dat_12$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((0) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((0) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((1) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((1) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((2) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((2) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((3) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((3) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((4) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((4) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((5) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((5) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C17. CNE function through the negative clairvoyance & positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_ncpc <- function(par, x, z, n, b, db1, db2, db3, db4, db5){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_ncpc <- function(par, x, z, n, x100, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_ncpc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_ncpc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_ncpc = c(mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_ncpc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_ncpc = numeric(100) # Prepare profile of likelihood values

pars_x100_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4])
x100vec_ncpc = seq(-120, 100, length = 100)
x100prof_ncpc = numeric(100)

pars_c_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[4])
cvec_ncpc = seq(0.1, 50, length = 100)
cprof_ncpc = numeric(100)

pars_theta_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3])
thetavec_ncpc = seq(0.01, 10, length = 100)
thetaprof_ncpc = numeric(100)

# Optimize original function without fixed parameters
pars_2_ncpc <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_ncpc <- optim(par = pars_2_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n,
                           db1 = dat_12$year >= 2014, 
                           db2 = dat_12$year >= 2015, 
                           db3 = dat_12$year >= 2016, 
                           db4 = dat_12$year >= 2017, 
                           db5 = dat_12$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_ncpc <- list(rep(list(), times = 100))
x100_coefs_vec_ncpc <- list(rep(list(), times = 100))
c_coefs_vec_ncpc <- list(rep(list(), times = 100))
theta_coefs_vec_ncpc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_ncpc <- optim(par = pars_b_ncpc, fn = ci_b_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, b = bvec_ncpc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_12$year >= 2014, 
                    db2 = dat_12$year >= 2015, 
                    db3 = dat_12$year >= 2016, 
                    db4 = dat_12$year >= 2017, 
                    db5 = dat_12$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_ncpc[[j]] = opt_b_ncpc$par # Store parameters for control
  # print(opt_a$convergence)
  bprof_ncpc[j] <- opt_b_ncpc$value # Set likelihood value for profile
  
  opt_x100_ncpc <- optim(par = pars_x100_ncpc, fn = ci_x100_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, x100 = x100vec_ncpc[j],
                       db1 = dat_12$year >= 2014, 
                       db2 = dat_12$year >= 2015, 
                       db3 = dat_12$year >= 2016, 
                       db4 = dat_12$year >= 2017, 
                       db5 = dat_12$year >= 2018, 
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_ncpc[j] <- opt_x100_ncpc$value
  x100_coefs_vec_ncpc[[j]] = opt_x100_ncpc$par
  
  opt_c_ncpc <- optim(par = pars_c_ncpc, fn = ci_c_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, c = cvec_ncpc[j],
                    db1 = dat_12$year >= 2014, 
                    db2 = dat_12$year >= 2015, 
                    db3 = dat_12$year >= 2016, 
                    db4 = dat_12$year >= 2017, 
                    db5 = dat_12$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_ncpc[j] <- opt_c_ncpc$value
  c_coefs_vec_ncpc[[j]] = opt_c_ncpc$par
  
  opt_theta_ncpc <- optim(par = pars_theta_ncpc, fn = ci_theta_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, theta = thetavec_ncpc[j], 
                        db1 = dat_12$year >= 2014, 
                        db2 = dat_12$year >= 2015, 
                        db3 = dat_12$year >= 2016, 
                        db4 = dat_12$year >= 2017, 
                        db5 = dat_12$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_ncpc[j] <- opt_theta_ncpc$value
  theta_coefs_vec_ncpc[[j]] = opt_theta_ncpc$par
}

# Set the 95% confidence limits
bprof_lower_ncpc <- bprof_ncpc[1:which.min(bprof_ncpc)] # Likelihood values below the best estimate
bvec_lower_ncpc <- bvec_ncpc[1:which.min(bprof_ncpc)] # Parameter values below the best estimate
bprof_higher_ncpc <- bprof_ncpc[which.min(bprof_ncpc):length(bprof_ncpc)] # Likelihood values above the best estimate
bvec_higher_ncpc <- bvec_ncpc[which.min(bprof_ncpc):length(bprof_ncpc)] # Parameter values above the best estimate
l_b_ci_ncpc <- approx(bprof_lower_ncpc, bvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_ncpc <- approx(bprof_higher_ncpc, bvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_ncpc <- x100prof_ncpc[1:which.min(x100prof_ncpc)]
x100vec_lower_ncpc <- x100vec_ncpc[1:which.min(x100prof_ncpc)]
x100prof_higher_ncpc <- x100prof_ncpc[which.min(x100prof_ncpc):length(x100prof_ncpc)]
x100vec_higher_ncpc <- x100vec_ncpc[which.min(x100prof_ncpc):length(x100prof_ncpc)]
l_x100_ci_ncpc <- approx(x100prof_lower_ncpc, x100vec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_x100_ci_ncpc <- approx(x100prof_higher_ncpc, x100vec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

cprof_lower_ncpc <- cprof_ncpc[1:which.min(cprof_ncpc)]
cvec_lower_ncpc <- cvec_ncpc[1:which.min(cprof_ncpc)]
cprof_higher_ncpc <- cprof_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
cvec_higher_ncpc <- cvec_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
l_c_ci_ncpc <- approx(cprof_lower_ncpc, cvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_c_ci_ncpc <- approx(cprof_higher_ncpc, cvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

thetaprof_lower_ncpc <- thetaprof_ncpc[1:which.min(thetaprof_ncpc)]
thetavec_lower_ncpc <- thetavec_ncpc[1:which.min(thetaprof_ncpc)]
thetaprof_higher_ncpc <- thetaprof_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
thetavec_higher_ncpc <- thetavec_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
l_theta_ci_ncpc <- approx(thetaprof_lower_ncpc, thetavec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_theta_ci_ncpc <- approx(thetaprof_higher_ncpc, thetavec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_ncpc_cne_yr <- c(opt_log_bbinom_ncpc$par[1], l_b_ci_ncpc$y, r_b_ci_ncpc$y) 
x100_ci_ncpc_cne_yr <- c(opt_log_bbinom_ncpc$par[2], l_x100_ci_ncpc$y, r_x100_ci_ncpc$y)
c_ci_ncpc_cne_yr <- c(opt_log_bbinom_ncpc$par[3], l_c_ci_ncpc$y, r_c_ci_ncpc$y)
theta_ci_ncpc_cne_yr <- c(opt_log_bbinom_ncpc$par[4], l_theta_ci_ncpc$y, r_theta_ci_ncpc$y)

b_ci_ncpc_cne_yr
x100_ci_ncpc_cne_yr
c_ci_ncpc_cne_yr
theta_ci_ncpc_cne_yr


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function
mlsp_fun_ncpc <- function(par, x, z, n, db1, db2, db3, db4, db5){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_ncpc <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_ncpc <- optim(par = mlsp_par_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, 
                     db1 = dat_12$year >= 2014, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_12$year >= 2015, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_12$year >= 2016, 
                     db4 = dat_12$year >= 2017, 
                     db5 = dat_12$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_ncpc <- mlsp_opt_ncpc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_12[dat_12$year == 2013,], aes(x = dist, y = prop, color = "2013"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2014,], aes(x = dist, y = prop, color = "2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2015,], aes(x = dist, y = prop, color = "2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2016,], aes(x = dist, y = prop, color = "2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2017,], aes(x = dist, y = prop, color = "2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$year == 2018,], aes(x = dist, y = prop, color = "2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(0)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(0))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(0))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(1)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(1))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(1))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(2)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(2))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(2))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(3)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(3))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(3))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(4)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(4))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(4))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(5)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(5))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(5))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C18. Hill function through the negative clairvoyance & positive carry-over data with years.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_ncpc <- function(par, x, z, n, a, db1, db2, db3, db4, db5){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_ncpc <- function(par, x, z, n, pn, db1, db2, db3, db4, db5){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_ncpc <- function(par, x, z, n, h, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_ncpc <- function(par, x, z, n, c, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_ncpc <- function(par, x, z, n, theta, db1, db2, db3, db4, db5){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4 + 1*db5))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_ncpc = c(mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4], mlsp_coefs_ncpc[5]) # Starting parameters are the estimated parameters of the optimized function
avec_ncpc = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_ncpc = numeric(100) # Prepare profile of likelihood values

pars_pn_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4], mlsp_coefs_ncpc[5])
pnvec_ncpc = seq(1.5, 8, length = 100)
pnprof_ncpc = numeric(100)

pars_h_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[4], mlsp_coefs_ncpc[5])
hvec_ncpc = seq(2, 100, length = 100)
hprof_ncpc = numeric(100)

pars_c_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[5])
cvec_ncpc = seq(0.1, 50, length = 100)
cprof_ncpc = numeric(100)

pars_theta_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4])
thetavec_ncpc = seq(0.01, 10, length = 100)
thetaprof_ncpc = numeric(100)

# Optimize original function without fixed parameters
pars_2_ncpc <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_ncpc <- optim(par = pars_2_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n,
                           db1 = dat_12$year >= 2014, 
                           db2 = dat_12$year >= 2015, 
                           db3 = dat_12$year >= 2016, 
                           db4 = dat_12$year >= 2017, 
                           db5 = dat_12$year >= 2018, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_ncpc <- list(rep(list(), times = 100))
pn_coefs_vec_ncpc <- list(rep(list(), times = 100))
h_coefs_vec_ncpc <- list(rep(list(), times = 100))
c_coefs_vec_ncpc <- list(rep(list(), times = 100))
theta_coefs_vec_ncpc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_ncpc <- optim(par = pars_a_ncpc, fn = ci_a_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, a = avec_ncpc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_12$year >= 2014, 
                    db2 = dat_12$year >= 2015, 
                    db3 = dat_12$year >= 2016, 
                    db4 = dat_12$year >= 2017, 
                    db5 = dat_12$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_ncpc[[j]] = opt_a_ncpc$par # Store parameters for control
  aprof_ncpc[j] <- opt_a_ncpc$value # Set likelihood value for profile
  
  opt_pn_ncpc <- optim(par = pars_pn_ncpc, fn = ci_pn_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, pn = pnvec_ncpc[j],
                     db1 = dat_12$year >= 2014, 
                     db2 = dat_12$year >= 2015, 
                     db3 = dat_12$year >= 2016, 
                     db4 = dat_12$year >= 2017, 
                     db5 = dat_12$year >= 2018, 
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_ncpc[j] <- opt_pn_ncpc$value
  pn_coefs_vec_ncpc[[j]] = opt_pn_ncpc$par
  
  opt_h_ncpc <- optim(par = pars_h_ncpc, fn = ci_h_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, h = hvec_ncpc[j],
                    db1 = dat_12$year >= 2014, 
                    db2 = dat_12$year >= 2015, 
                    db3 = dat_12$year >= 2016, 
                    db4 = dat_12$year >= 2017, 
                    db5 = dat_12$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_ncpc[j] <- opt_h_ncpc$value
  h_coefs_vec_ncpc[[j]] = opt_h_ncpc$par
  
  opt_c_ncpc <- optim(par = pars_c_ncpc, fn = ci_c_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, c = cvec_ncpc[j],
                    db1 = dat_12$year >= 2014, 
                    db2 = dat_12$year >= 2015, 
                    db3 = dat_12$year >= 2016, 
                    db4 = dat_12$year >= 2017, 
                    db5 = dat_12$year >= 2018, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_ncpc[j] <- opt_c_ncpc$value
  c_coefs_vec_ncpc[[j]] = opt_c_ncpc$par
  
  opt_theta_ncpc <- optim(par = pars_theta_ncpc, fn = ci_theta_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, theta = thetavec_ncpc[j], 
                        db1 = dat_12$year >= 2014, 
                        db2 = dat_12$year >= 2015, 
                        db3 = dat_12$year >= 2016, 
                        db4 = dat_12$year >= 2017, 
                        db5 = dat_12$year >= 2018, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_ncpc[j] <- opt_theta_ncpc$value
  theta_coefs_vec_ncpc[[j]] = opt_theta_ncpc$par
}

# Set the 95% confidence limits
aprof_lower_ncpc <- aprof_ncpc[1:which.min(aprof_ncpc)] # Likelihood values below the best estimate
avec_lower_ncpc <- avec_ncpc[1:which.min(aprof_ncpc)] # Parameter values below the best estimate
aprof_higher_ncpc <- aprof_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Likelihood values above the best estimate
avec_higher_ncpc <- avec_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Parameter values above the best estimate
l_a_ci_ncpc <- approx(aprof_lower_ncpc, avec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_ncpc <- approx(aprof_higher_ncpc, avec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_ncpc <- pnprof_ncpc[1:which.min(pnprof_ncpc)]
pnvec_lower_ncpc <- pnvec_ncpc[1:which.min(pnprof_ncpc)]
pnprof_higher_ncpc <- pnprof_ncpc[which.min(pnprof_ncpc):length(pnprof_ncpc)]
pnvec_higher_ncpc <- pnvec_ncpc[which.min(pnprof_ncpc):length(pnprof_ncpc)]
l_pn_ci_ncpc <- approx(pnprof_lower_ncpc, pnvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_pn_ci_ncpc <- approx(pnprof_higher_ncpc, pnvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

hprof_lower_ncpc <- hprof_ncpc[1:which.min(hprof_ncpc)]
hvec_lower_ncpc <- hvec_ncpc[1:which.min(hprof_ncpc)]
hprof_higher_ncpc <- hprof_ncpc[which.min(hprof_ncpc):length(hprof_ncpc)]
hvec_higher_ncpc <- hvec_ncpc[which.min(hprof_ncpc):length(hprof_ncpc)]
l_h_ci_ncpc <- approx(hprof_lower_ncpc, hvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_h_ci_ncpc <- approx(hprof_higher_ncpc, hvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

cprof_lower_ncpc <- cprof_ncpc[1:which.min(cprof_ncpc)]
cvec_lower_ncpc <- cvec_ncpc[1:which.min(cprof_ncpc)]
cprof_higher_ncpc <- cprof_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
cvec_higher_ncpc <- cvec_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
l_c_ci_ncpc <- approx(cprof_lower_ncpc, cvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_c_ci_ncpc <- approx(cprof_higher_ncpc, cvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

thetaprof_lower_ncpc <- thetaprof_ncpc[1:which.min(thetaprof_ncpc)]
thetavec_lower_ncpc <- thetavec_ncpc[1:which.min(thetaprof_ncpc)]
thetaprof_higher_ncpc <- thetaprof_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
thetavec_higher_ncpc <- thetavec_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
l_theta_ci_ncpc <- approx(thetaprof_lower_ncpc, thetavec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_theta_ci_ncpc <- approx(thetaprof_higher_ncpc, thetavec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_ncpc_hill_yr <- c(opt_log_bbinom_ncpc$par[1], l_a_ci_ncpc$y, r_a_ci_ncpc$y) 
pn_ci_ncpc_hill_yr <- c(opt_log_bbinom_ncpc$par[2], l_pn_ci_ncpc$y, r_pn_ci_ncpc$y)
h_ci_ncpc_hill_yr <- c(opt_log_bbinom_ncpc$par[3], l_h_ci_ncpc$y, r_h_ci_ncpc$y)
c_ci_ncpc_hill_yr <- c(opt_log_bbinom_ncpc$par[4], l_c_ci_ncpc$y, r_c_ci_ncpc$y)
theta_ci_ncpc_hill_yr <- c(opt_log_bbinom_ncpc$par[5], l_theta_ci_ncpc$y, r_theta_ci_ncpc$y)

a_ci_ncpc_hill_yr
pn_ci_ncpc_hill_yr
h_ci_ncpc_hill_yr
c_ci_ncpc_hill_yr
theta_ci_ncpc_hill_yr


#######################################
#######################################
######                           ######
######           SEASONS         ######
######                           ######
#######################################
#######################################

###############################
##                           ##
## Data reading and morphing ##
##                           ##
###############################

# Set date factor correctly, drop na and incorrect lon/lat values
dat <- dat_base
dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter(lon != 0 & lat != 0)

dat <- dat %>% 
  mutate(season = ifelse((dat$year == 2013) | (dat$year == 2014 & dat$month %in% 1:3), 1, 
                         ifelse((dat$year == 2014 & dat$month %in% 4:12) | (dat$year == 2015 & dat$month %in% 1:3), 2, 
                                ifelse((dat$year == 2015 & dat$month %in% 4:12) | (dat$year == 2016 & dat$month %in% 1:3), 3, 
                                       ifelse((dat$year == 2016 & dat$month %in% 4:12) | (dat$year == 2017 & dat$month %in% 1:3), 4, 5)))))


# Set number of seasons, used for length of loops and such
seasons <- unique(dat$season)

# Assumed origin
gallipoli <- c(lon=17.992615,lat=40.055851)

# Calculate distance of data points to origin
dat <- dat %>%
  mutate(dist = round((distHaversine(tibble(lon = dat$lon, lat = dat$lat), 
                                     c(gallipoli[1], gallipoli[2])))/1000))

# Find the maximum distance measured for the distance circels and the maximum distance where a positive is found for visualization
max_dist <- max(dat$dist)
max_pos <- max(dat[dat$result == 1,]$dist)

# Split the data for each year
dat <- dat %>% 
  split(dat$season)


###################
##               ##
## Model fitting ##
##               ##
###################

#################
# Original data #
#################

#
# Data morphing and distance circle creation
#

# Set distance circles
dc <- 0:(max_dist - 1)
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))

dat_2 <- list()

for(i in 1:length(seasons)){
  for(j in dc){
    pos[[i]][j + 1] <- sum(dat[[i]]$result[dat[[i]]$dist >= j & dat[[i]]$dist < j+1])
    n[[i]][j + 1] <- length(dat[[i]]$result[dat[[i]]$dist >= j & dat[[i]]$dist < j+1])
  }
  dat_2[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n)
}

# Create dataset with seasons for breaks
dat_3 <- list()
for(i in 1:length(seasons)){
  dat_2[[i]] <- dat_2[[i]] %>% 
    mutate(season = (i))
  dat_3 <- rbind(dat_3, dat_2[[i]])
}

# Plot dat_3
ggplot() +
  geom_point(data = dat_3[dat_3$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  labs(
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C19. Data points of the original data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function
mlsp_fun_od <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- (1/(1 + exp(par[1] * (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3]))))) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
  return(nll)
}

# Setup parameter starting values
mlsp_par_od <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_od <- optim(par = mlsp_par_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, 
                     db1 = dat_3$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_3$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_3$season >= 4, 
                     db4 = dat_3$season >= 5, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_od <- mlsp_opt_od$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_3[dat_3$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 0*mlsp_coefs_od[3]))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 1*mlsp_coefs_od[3]))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 2*mlsp_coefs_od[3]))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 3*mlsp_coefs_od[3]))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_od[1]*(x - (mlsp_coefs_od[2] + 4*mlsp_coefs_od[3]))))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C20. Logistic function through the original data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_od <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_od <- function(par, x, z, n, x50, db1, db2, db3, db4){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_od <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_od <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4)*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_od = c(mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4]) # Starting parameters are the estimated parameters of the optimized function
avec_od = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_od = numeric(100) # Prepare profile of likelihood values

pars_x50_od = c(mlsp_coefs_od[1], mlsp_coefs_od[3], mlsp_coefs_od[4])
x50vec_od = seq(-120, 100, length = 100)
x50prof_od = numeric(100)

pars_c_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[4])
cvec_od = seq(0.1, 50, length = 100)
cprof_od = numeric(100)

pars_theta_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3])
thetavec_od = seq(0.01, 10, length = 100)
thetaprof_od = numeric(100)

# Optimize original function without fixed parameters
pars_2_od <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_od <- optim(par = pars_2_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n,
                           db1 = dat_3$season >= 2, 
                           db2 = dat_3$season >= 3, 
                           db3 = dat_3$season >= 4, 
                           db4 = dat_3$season >= 5, 
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_od <- list(rep(list(), times = 100))
x50_coefs_vec_od <- list(rep(list(), times = 100))
c_coefs_vec_od <- list(rep(list(), times = 100))
theta_coefs_vec_od <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_od <- optim(par = pars_a_od, fn = ci_a_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, a = avec_od[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_3$season >= 2, 
                    db2 = dat_3$season >= 3, 
                    db3 = dat_3$season >= 4, 
                    db4 = dat_3$season >= 5, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_od[[j]] = opt_a_od$par # Store parameters for control
  # print(opt_a$convergence)
  aprof_od[j] <- opt_a_od$value # Set likelihood value for profile
  
  opt_x50_od <- optim(par = pars_x50_od, fn = ci_x50_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, x50 = x50vec_od[j],
                      db1 = dat_3$season >= 2, 
                      db2 = dat_3$season >= 3, 
                      db3 = dat_3$season >= 4, 
                      db4 = dat_3$season >= 5, 
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_od[j] <- opt_x50_od$value
  x50_coefs_vec_od[[j]] = opt_x50_od$par
  
  opt_c_od <- optim(par = pars_c_od, fn = ci_c_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, c = cvec_od[j],
                    db1 = dat_3$season >= 2, 
                    db2 = dat_3$season >= 3, 
                    db3 = dat_3$season >= 4, 
                    db4 = dat_3$season >= 5, 
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_od[j] <- opt_c_od$value
  c_coefs_vec_od[[j]] = opt_c_od$par
  
  opt_theta_od <- optim(par = pars_theta_od, fn = ci_theta_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, theta = thetavec_od[j], 
                        db1 = dat_3$season >= 2, 
                        db2 = dat_3$season >= 3, 
                        db3 = dat_3$season >= 4, 
                        db4 = dat_3$season >= 5, 
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_od[j] <- opt_theta_od$value
  theta_coefs_vec_od[[j]] = opt_theta_od$par
}

# Set the 95% confidence limits
aprof_lower_od <- aprof_od[1:which.min(aprof_od)] # Likelihood values below the best estimate
avec_lower_od <- avec_od[1:which.min(aprof_od)] # Parameter values below the best estimate
aprof_higher_od <- aprof_od[which.min(aprof_od):length(aprof_od)] # Likelihood values above the best estimate
avec_higher_od <- avec_od[which.min(aprof_od):length(aprof_od)] # Parameter values above the best estimate
l_a_ci_od <- approx(aprof_lower_od, avec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_od <- approx(aprof_higher_od, avec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x50prof_lower_od <- x50prof_od[1:which.min(x50prof_od)]
x50vec_lower_od <- x50vec_od[1:which.min(x50prof_od)]
x50prof_higher_od <- x50prof_od[which.min(x50prof_od):length(x50prof_od)]
x50vec_higher_od <- x50vec_od[which.min(x50prof_od):length(x50prof_od)]
l_x50_ci_od <- approx(x50prof_lower_od, x50vec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_x50_ci_od <- approx(x50prof_higher_od, x50vec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

cprof_lower_od <- cprof_od[1:which.min(cprof_od)]
cvec_lower_od <- cvec_od[1:which.min(cprof_od)]
cprof_higher_od <- cprof_od[which.min(cprof_od):length(cprof_od)]
cvec_higher_od <- cvec_od[which.min(cprof_od):length(cprof_od)]
l_c_ci_od <- approx(cprof_lower_od, cvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_c_ci_od <- approx(cprof_higher_od, cvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

thetaprof_lower_od <- thetaprof_od[1:which.min(thetaprof_od)]
thetavec_lower_od <- thetavec_od[1:which.min(thetaprof_od)]
thetaprof_higher_od <- thetaprof_od[which.min(thetaprof_od):length(thetaprof_od)]
thetavec_higher_od <- thetavec_od[which.min(thetaprof_od):length(thetaprof_od)]
l_theta_ci_od <- approx(thetaprof_lower_od, thetavec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_theta_ci_od <- approx(thetaprof_higher_od, thetavec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_od_log_sn <- c(opt_log_bbinom_od$par[1], l_a_ci_od$y, r_a_ci_od$y) 
x50_ci_od_log_sn <- c(opt_log_bbinom_od$par[2], l_x50_ci_od$y, r_x50_ci_od$y)
c_ci_od_log_sn <- c(opt_log_bbinom_od$par[3], l_c_ci_od$y, r_c_ci_od$y)
theta_ci_od_log_sn <- c(opt_log_bbinom_od$par[4], l_theta_ci_od$y, r_theta_ci_od$y)

a_ci_od_log_sn
x50_ci_od_log_sn
c_ci_od_log_sn
theta_ci_od_log_sn


##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function
mlsp_fun_od <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_od <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_od <- optim(par = mlsp_par_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, 
                     db1 = dat_3$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_3$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_3$season >= 4, 
                     db4 = dat_3$season >= 5, 
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_od <- mlsp_opt_od$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_3[dat_3$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((0) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((0) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((1) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((1) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((2) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((2) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((3) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((3) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((4) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_od[1]*(mlsp_coefs_od[2] + ((4) * mlsp_coefs_od[3])))) * exp(-mlsp_coefs_od[1] * x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C21. CNE function through the original data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_od <- function(par, x, z, n, b, db1, db2, db3, db4){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_od <- function(par, x, z, n, x100, db1, db2, db3, db4){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_od <- function(par, x, z, n, c, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_od <- function(par, x, z, n, theta, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_od = c(mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_od = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_od = numeric(100) # Prepare profile of likelihood values

pars_x100_od = c(mlsp_coefs_od[1], mlsp_coefs_od[3], mlsp_coefs_od[4])
x100vec_od = seq(-120, 100, length = 100)
x100prof_od = numeric(100)

pars_c_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[4])
cvec_od = seq(0.1, 50, length = 100)
cprof_od = numeric(100)

pars_theta_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3])
thetavec_od = seq(0.01, 10, length = 100)
thetaprof_od = numeric(100)

# Optimize original function without fixed parameters
pars_2_od <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_od <- optim(par = pars_2_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n,
                           db1 = dat_3$season >= 2, 
                           db2 = dat_3$season >= 3, 
                           db3 = dat_3$season >= 4, 
                           db4 = dat_3$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_od <- list(rep(list(), times = 100))
x100_coefs_vec_od <- list(rep(list(), times = 100))
c_coefs_vec_od <- list(rep(list(), times = 100))
theta_coefs_vec_od <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_od <- optim(par = pars_b_od, fn = ci_b_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, b = bvec_od[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_3$season >= 2, 
                    db2 = dat_3$season >= 3, 
                    db3 = dat_3$season >= 4, 
                    db4 = dat_3$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_od[[j]] = opt_b_od$par # Store parameters for control
  # print(opt_a$convergence)
  bprof_od[j] <- opt_b_od$value # Set likelihood value for profile
  
  opt_x100_od <- optim(par = pars_x100_od, fn = ci_x100_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, x100 = x100vec_od[j],
                       db1 = dat_3$season >= 2, 
                       db2 = dat_3$season >= 3, 
                       db3 = dat_3$season >= 4, 
                       db4 = dat_3$season >= 5, 
                       
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_od[j] <- opt_x100_od$value
  x100_coefs_vec_od[[j]] = opt_x100_od$par
  
  opt_c_od <- optim(par = pars_c_od, fn = ci_c_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, c = cvec_od[j],
                    db1 = dat_3$season >= 2, 
                    db2 = dat_3$season >= 3, 
                    db3 = dat_3$season >= 4, 
                    db4 = dat_3$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_od[j] <- opt_c_od$value
  c_coefs_vec_od[[j]] = opt_c_od$par
  
  opt_theta_od <- optim(par = pars_theta_od, fn = ci_theta_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, theta = thetavec_od[j], 
                        db1 = dat_3$season >= 2, 
                        db2 = dat_3$season >= 3, 
                        db3 = dat_3$season >= 4, 
                        db4 = dat_3$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_od[j] <- opt_theta_od$value
  theta_coefs_vec_od[[j]] = opt_theta_od$par
}

# Set the 95% confidence limits
bprof_lower_od <- bprof_od[1:which.min(bprof_od)] # Likelihood values below the best estimate
bvec_lower_od <- bvec_od[1:which.min(bprof_od)] # Parameter values below the best estimate
bprof_higher_od <- bprof_od[which.min(bprof_od):length(bprof_od)] # Likelihood values above the best estimate
bvec_higher_od <- bvec_od[which.min(bprof_od):length(bprof_od)] # Parameter values above the best estimate
l_b_ci_od <- approx(bprof_lower_od, bvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_od <- approx(bprof_higher_od, bvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_od <- x100prof_od[1:which.min(x100prof_od)]
x100vec_lower_od <- x100vec_od[1:which.min(x100prof_od)]
x100prof_higher_od <- x100prof_od[which.min(x100prof_od):length(x100prof_od)]
x100vec_higher_od <- x100vec_od[which.min(x100prof_od):length(x100prof_od)]
l_x100_ci_od <- approx(x100prof_lower_od, x100vec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_x100_ci_od <- approx(x100prof_higher_od, x100vec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

cprof_lower_od <- cprof_od[1:which.min(cprof_od)]
cvec_lower_od <- cvec_od[1:which.min(cprof_od)]
cprof_higher_od <- cprof_od[which.min(cprof_od):length(cprof_od)]
cvec_higher_od <- cvec_od[which.min(cprof_od):length(cprof_od)]
l_c_ci_od <- approx(cprof_lower_od, cvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_c_ci_od <- approx(cprof_higher_od, cvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

thetaprof_lower_od <- thetaprof_od[1:which.min(thetaprof_od)]
thetavec_lower_od <- thetavec_od[1:which.min(thetaprof_od)]
thetaprof_higher_od <- thetaprof_od[which.min(thetaprof_od):length(thetaprof_od)]
thetavec_higher_od <- thetavec_od[which.min(thetaprof_od):length(thetaprof_od)]
l_theta_ci_od <- approx(thetaprof_lower_od, thetavec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)
r_theta_ci_od <- approx(thetaprof_higher_od, thetavec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_od_cne_sn <- c(opt_log_bbinom_od$par[1], l_b_ci_od$y, r_b_ci_od$y) 
x100_ci_od_cne_sn <- c(opt_log_bbinom_od$par[2], l_x100_ci_od$y, r_x100_ci_od$y)
c_ci_od_cne_sn <- c(opt_log_bbinom_od$par[3], l_c_ci_od$y, r_c_ci_od$y)
theta_ci_od_cne_sn <- c(opt_log_bbinom_od$par[4], l_theta_ci_od$y, r_theta_ci_od$y)

b_ci_od_cne_sn
x100_ci_od_cne_sn
c_ci_od_cne_sn
theta_ci_od_cne_sn


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function
mlsp_fun_od <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_od <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_od <- optim(par = mlsp_par_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, 
                     db1 = dat_3$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_3$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_3$season >= 4, 
                     db4 = dat_3$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_od <- mlsp_opt_od$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_3[dat_3$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_3[dat_3$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(0)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(0))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(0))^mlsp_coefs_od[2])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(1)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(1))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(1))^mlsp_coefs_od[2])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(2)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(2))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(2))^mlsp_coefs_od[2])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(3)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(3))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(3))^mlsp_coefs_od[2])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_od[4]*(4)), mlsp_coefs_od[1], 
                                        (mlsp_coefs_od[1]*(1 - ((x - mlsp_coefs_od[4]*(4))^mlsp_coefs_od[2])/
                                                   (mlsp_coefs_od[3]^mlsp_coefs_od[2] + (x - mlsp_coefs_od[4]*(4))^mlsp_coefs_od[2])))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C21. Hill function through the original data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_od <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_od <- function(par, x, z, n, pn, db1, db2, db3, db4){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_od <- function(par, x, z, n, h, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_od <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_od <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_od = c(mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4], mlsp_coefs_od[5]) # Starting parameters are the estimated parameters of the optimized function
avec_od = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_od = numeric(100) # Prepare profile of likelihood values

pars_pn_od = c(mlsp_coefs_od[1], mlsp_coefs_od[3], mlsp_coefs_od[4], mlsp_coefs_od[5])
pnvec_od = seq(0.1, 10, length = 100)
pnprof_od = numeric(100)

pars_h_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[4], mlsp_coefs_od[5])
hvec_od = seq(5, 100, length = 100)
hprof_od = numeric(100)

pars_c_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[5])
cvec_od = seq(0.1, 50, length = 100)
cprof_od = numeric(100)

pars_theta_od = c(mlsp_coefs_od[1], mlsp_coefs_od[2], mlsp_coefs_od[3], mlsp_coefs_od[4])
thetavec_od = seq(0.01, 10, length = 100)
thetaprof_od = numeric(100)

# Optimize original function without fixed parameters
pars_2_od <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_od <- optim(par = pars_2_od, fn = mlsp_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n,
                           db1 = dat_3$season >= 2, 
                           db2 = dat_3$season >= 3, 
                           db3 = dat_3$season >= 4, 
                           db4 = dat_3$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_od <- list(rep(list(), times = 100))
pn_coefs_vec_od <- list(rep(list(), times = 100))
h_coefs_vec_od <- list(rep(list(), times = 100))
c_coefs_vec_od <- list(rep(list(), times = 100))
theta_coefs_vec_od <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_od <- optim(par = pars_a_od, fn = ci_a_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, a = avec_od[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_3$season >= 2, 
                    db2 = dat_3$season >= 3, 
                    db3 = dat_3$season >= 4, 
                    db4 = dat_3$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_od[[j]] = opt_a_od$par # Store parameters for control
  aprof_od[j] <- opt_a_od$value # Set likelihood value for profile
  
  opt_pn_od <- optim(par = pars_pn_od, fn = ci_pn_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, pn = pnvec_od[j],
                     db1 = dat_3$season >= 2, 
                     db2 = dat_3$season >= 3, 
                     db3 = dat_3$season >= 4, 
                     db4 = dat_3$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_od[j] <- opt_pn_od$value
  pn_coefs_vec_od[[j]] = opt_pn_od$par
  
  opt_h_od <- optim(par = pars_h_od, fn = ci_h_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, h = hvec_od[j],
                    db1 = dat_3$season >= 2, 
                    db2 = dat_3$season >= 3, 
                    db3 = dat_3$season >= 4, 
                    db4 = dat_3$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_od[j] <- opt_h_od$value
  h_coefs_vec_od[[j]] = opt_h_od$par
  
  opt_c_od <- optim(par = pars_c_od, fn = ci_c_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, c = cvec_od[j],
                    db1 = dat_3$season >= 2, 
                    db2 = dat_3$season >= 3, 
                    db3 = dat_3$season >= 4, 
                    db4 = dat_3$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_od[j] <- opt_c_od$value
  c_coefs_vec_od[[j]] = opt_c_od$par
  
  opt_theta_od <- optim(par = pars_theta_od, fn = ci_theta_fun_od, x = dat_3$dist, z = dat_3$pos, n = dat_3$n, theta = thetavec_od[j], 
                        db1 = dat_3$season >= 2, 
                        db2 = dat_3$season >= 3, 
                        db3 = dat_3$season >= 4, 
                        db4 = dat_3$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_od[j] <- opt_theta_od$value
  theta_coefs_vec_od[[j]] = opt_theta_od$par
}

# Set the 95% confidence limits
aprof_lower_od <- aprof_od[1:which.min(aprof_od)] # Likelihood values below the best estimate
avec_lower_od <- avec_od[1:which.min(aprof_od)] # Parameter values below the best estimate
aprof_higher_od <- aprof_od[which.min(aprof_od):length(aprof_od)] # Likelihood values above the best estimate
avec_higher_od <- avec_od[which.min(aprof_od):length(aprof_od)] # Parameter values above the best estimate
l_a_ci_od <- approx(aprof_lower_od, avec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_od <- approx(aprof_higher_od, avec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_od <- pnprof_od[1:which.min(pnprof_od)]
pnvec_lower_od <- pnvec_od[1:which.min(pnprof_od)]
pnprof_higher_od <- pnprof_od[which.min(pnprof_od):length(pnprof_od)]
pnvec_higher_od <- pnvec_od[which.min(pnprof_od):length(pnprof_od)]
l_pn_ci_od <- approx(pnprof_lower_od, pnvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_pn_ci_od <- approx(pnprof_higher_od, pnvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

hprof_lower_od <- hprof_od[1:which.min(hprof_od)]
hvec_lower_od <- hvec_od[1:which.min(hprof_od)]
hprof_higher_od <- hprof_od[which.min(hprof_od):length(hprof_od)]
hvec_higher_od <- hvec_od[which.min(hprof_od):length(hprof_od)]
l_h_ci_od <- approx(hprof_lower_od, hvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_h_ci_od <- approx(hprof_higher_od, hvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

cprof_lower_od <- cprof_od[1:which.min(cprof_od)]
cvec_lower_od <- cvec_od[1:which.min(cprof_od)]
cprof_higher_od <- cprof_od[which.min(cprof_od):length(cprof_od)]
cvec_higher_od <- cvec_od[which.min(cprof_od):length(cprof_od)]
l_c_ci_od <- approx(cprof_lower_od, cvec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_c_ci_od <- approx(cprof_higher_od, cvec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

thetaprof_lower_od <- thetaprof_od[1:which.min(thetaprof_od)]
thetavec_lower_od <- thetavec_od[1:which.min(thetaprof_od)]
thetaprof_higher_od <- thetaprof_od[which.min(thetaprof_od):length(thetaprof_od)]
thetavec_higher_od <- thetavec_od[which.min(thetaprof_od):length(thetaprof_od)]
l_theta_ci_od <- approx(thetaprof_lower_od, thetavec_lower_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)
r_theta_ci_od <- approx(thetaprof_higher_od, thetavec_higher_od, xout = opt_log_bbinom_od$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_od_hill_sn <- c(opt_log_bbinom_od$par[1], l_a_ci_od$y, r_a_ci_od$y) 
pn_ci_od_hill_sn <- c(opt_log_bbinom_od$par[2], l_pn_ci_od$y, r_pn_ci_od$y)
h_ci_od_hill_sn <- c(opt_log_bbinom_od$par[3], l_h_ci_od$y, r_h_ci_od$y)
c_ci_od_hill_sn <- c(opt_log_bbinom_od$par[4], l_c_ci_od$y, r_c_ci_od$y)
theta_ci_od_hill_sn <- c(opt_log_bbinom_od$par[5], l_theta_ci_od$y, r_theta_ci_od$y)

a_ci_od_hill_sn
pn_ci_od_hill_sn
h_ci_od_hill_sn
c_ci_od_hill_sn
theta_ci_od_hill_sn


############################
# Positive carry-over data #
############################

#
# Data morphing and distance circle creation
#

dat_4 <- dat

# make positives carry over
for(i in 2:length(seasons)){ # For-loop through every year from 2014
  co <- which(dat_4[[i-1]]$result == 1) # Which lines of the previous year have a result of 1
  dat_4[[i]] <- rbind(dat_4[[i]], dat_4[[i-1]][co,]) # Include the lines of the previous year with a result of 1 with the current year
  co_2 <- which(dat_4[[i]]$lon == dat_4[[i-1]]$lon[co] & dat_4[[i]]$lat == dat_4[[i-1]]$lat[co] & dat_4[[i]]$season == (2012 + i)) # Which lines of this year have the same coordinates of the positives last year
  # print(co_2)
  if(length(co_2) >= 1){ # If this is not 0, then there has been a measurement at these coordinates this year as well as last year, or multiple times in a single year (with these data, only multiple measurements in 2013)
    dat_4[[i]] <- dat_4[[i]][-co_2,] # Remove all the duplicates. If the duplicate measurement from this year is a 0, it will still be a 1, because I assume that once a disease is there, it will not leave. I assume that 0 to be a false negative.
  }
}

# Set distance circles
dc <- 1:(max_dist) # Number of distance circles is the maximum distance a measurement is done.
pos <- rep(list(rep(NA, times = length(dc))), times = length(seasons)) # Prepare list for number of positives in each dc
n <- rep(list(rep(NA, times = length(dc))), times = length(seasons)) # Prepare list for total number of measurements in each dc

dat_5 <- list() # Prepare data list

for(i in 1:length(seasons)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j] <- sum((dat_4[[i]]$result[dat_4[[i]]$dist >= j & (dat_4[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j] <- length((dat_4[[i]]$result[dat_4[[i]]$dist >= j & (dat_4[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_5[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
}

# Combine seasons and create year column for year 
dat_6 <- list()
for(i in 1:length(seasons)){
  dat_5[[i]] <- dat_5[[i]] %>% 
    mutate(season = (i))
  dat_6 <- rbind(dat_6, dat_5[[i]])
}

# Plot the data
ggplot() +
  geom_point(data = dat_6[dat_6$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  labs(
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C22. Data points of the positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function
mlsp_fun_pc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- (1/(1 + exp(par[1] * (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3]))))) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
  return(nll)
}

# Setup parameter starting values
mlsp_par_pc <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_pc <- optim(par = mlsp_par_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, 
                     db1 = dat_6$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_6$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_6$season >= 4, 
                     db4 = dat_6$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_pc <- mlsp_opt_pc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_6[dat_6$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 0*mlsp_coefs_pc[3]))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 1*mlsp_coefs_pc[3]))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 2*mlsp_coefs_pc[3]))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 3*mlsp_coefs_pc[3]))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_pc[1]*(x - (mlsp_coefs_pc[2] + 4*mlsp_coefs_pc[3]))))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C23. Logistic function through the positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_pc <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_pc <- function(par, x, z, n, x50, db1, db2, db3, db4){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_pc <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_pc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_pc = c(mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4]) # Starting parameters are the estimated parameters of the optimized function
avec_pc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_pc = numeric(100) # Prepare profile of likelihood values

pars_x50_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[3], mlsp_coefs_pc[4])
x50vec_pc = seq(-120, 100, length = 100)
x50prof_pc = numeric(100)

pars_c_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[4])
cvec_pc = seq(0.1, 50, length = 100)
cprof_pc = numeric(100)

pars_theta_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3])
thetavec_pc = seq(0.01, 10, length = 100)
thetaprof_pc = numeric(100)

# Optimize original function without fixed parameters
pars_2_pc <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_pc <- optim(par = pars_2_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n,
                           db1 = dat_6$season >= 2, 
                           db2 = dat_6$season >= 3, 
                           db3 = dat_6$season >= 4, 
                           db4 = dat_6$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_pc <- list(rep(list(), times = 100))
x50_coefs_vec_pc <- list(rep(list(), times = 100))
c_coefs_vec_pc <- list(rep(list(), times = 100))
theta_coefs_vec_pc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_pc <- optim(par = pars_a_pc, fn = ci_a_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, a = avec_pc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_6$season >= 2, 
                    db2 = dat_6$season >= 3, 
                    db3 = dat_6$season >= 4, 
                    db4 = dat_6$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_pc[[j]] = opt_a_pc$par # Store parameters for control
  # print(opt_a$convergence)
  aprof_pc[j] <- opt_a_pc$value # Set likelihood value for profile
  
  opt_x50_pc <- optim(par = pars_x50_pc, fn = ci_x50_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, x50 = x50vec_pc[j],
                      db1 = dat_6$season >= 2, 
                      db2 = dat_6$season >= 3, 
                      db3 = dat_6$season >= 4, 
                      db4 = dat_6$season >= 5, 
                      
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_pc[j] <- opt_x50_pc$value
  x50_coefs_vec_pc[[j]] = opt_x50_pc$par
  
  opt_c_pc <- optim(par = pars_c_pc, fn = ci_c_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, c = cvec_pc[j],
                    db1 = dat_6$season >= 2, 
                    db2 = dat_6$season >= 3, 
                    db3 = dat_6$season >= 4, 
                    db4 = dat_6$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_pc[j] <- opt_c_pc$value
  c_coefs_vec_pc[[j]] = opt_c_pc$par
  
  opt_theta_pc <- optim(par = pars_theta_pc, fn = ci_theta_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, theta = thetavec_pc[j], 
                        db1 = dat_6$season >= 2, 
                        db2 = dat_6$season >= 3, 
                        db3 = dat_6$season >= 4, 
                        db4 = dat_6$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_pc[j] <- opt_theta_pc$value
  theta_coefs_vec_pc[[j]] = opt_theta_pc$par
}

# Set the 95% confidence limits
aprof_lower_pc <- aprof_pc[1:which.min(aprof_pc)] # Likelihood values below the best estimate
avec_lower_pc <- avec_pc[1:which.min(aprof_pc)] # Parameter values below the best estimate
aprof_higher_pc <- aprof_pc[which.min(aprof_pc):length(aprof_pc)] # Likelihood values above the best estimate
avec_higher_pc <- avec_pc[which.min(aprof_pc):length(aprof_pc)] # Parameter values above the best estimate
l_a_ci_pc <- approx(aprof_lower_pc, avec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_pc <- approx(aprof_higher_pc, avec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x50prof_lower_pc <- x50prof_pc[1:which.min(x50prof_pc)]
x50vec_lower_pc <- x50vec_pc[1:which.min(x50prof_pc)]
x50prof_higher_pc <- x50prof_pc[which.min(x50prof_pc):length(x50prof_pc)]
x50vec_higher_pc <- x50vec_pc[which.min(x50prof_pc):length(x50prof_pc)]
l_x50_ci_pc <- approx(x50prof_lower_pc, x50vec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_x50_ci_pc <- approx(x50prof_higher_pc, x50vec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

cprof_lower_pc <- cprof_pc[1:which.min(cprof_pc)]
cvec_lower_pc <- cvec_pc[1:which.min(cprof_pc)]
cprof_higher_pc <- cprof_pc[which.min(cprof_pc):length(cprof_pc)]
cvec_higher_pc <- cvec_pc[which.min(cprof_pc):length(cprof_pc)]
l_c_ci_pc <- approx(cprof_lower_pc, cvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_c_ci_pc <- approx(cprof_higher_pc, cvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

thetaprof_lower_pc <- thetaprof_pc[1:which.min(thetaprof_pc)]
thetavec_lower_pc <- thetavec_pc[1:which.min(thetaprof_pc)]
thetaprof_higher_pc <- thetaprof_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
thetavec_higher_pc <- thetavec_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
l_theta_ci_pc <- approx(thetaprof_lower_pc, thetavec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_theta_ci_pc <- approx(thetaprof_higher_pc, thetavec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_pc_log_sn <- c(opt_log_bbinom_pc$par[1], l_a_ci_pc$y, r_a_ci_pc$y) 
x50_ci_pc_log_sn <- c(opt_log_bbinom_pc$par[2], l_x50_ci_pc$y, r_x50_ci_pc$y)
c_ci_pc_log_sn <- c(opt_log_bbinom_pc$par[3], l_c_ci_pc$y, r_c_ci_pc$y)
theta_ci_pc_log_sn <- c(opt_log_bbinom_pc$par[4], l_theta_ci_pc$y, r_theta_ci_pc$y)

a_ci_pc_log_sn
x50_ci_pc_log_sn
c_ci_pc_log_sn
theta_ci_pc_log_sn


##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function
mlsp_fun_pc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_pc <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_pc <- optim(par = mlsp_par_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, 
                     db1 = dat_6$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_6$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_6$season >= 4, 
                     db4 = dat_6$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_pc <- mlsp_opt_pc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_6[dat_6$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((0) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((0) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((1) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((1) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((2) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((2) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((3) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((3) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((4) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_pc[1]*(mlsp_coefs_pc[2] + ((4) * mlsp_coefs_pc[3])))) * exp(-mlsp_coefs_pc[1] * x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C24. CNE function through the positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_pc <- function(par, x, z, n, b, db1, db2, db3, db4){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_pc <- function(par, x, z, n, x100, db1, db2, db3, db4){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_pc <- function(par, x, z, n, c, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_pc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_pc = c(mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_pc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_pc = numeric(100) # Prepare profile of likelihood values

pars_x100_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[3], mlsp_coefs_pc[4])
x100vec_pc = seq(-120, 100, length = 100)
x100prof_pc = numeric(100)

pars_c_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[4])
cvec_pc = seq(0.1, 50, length = 100)
cprof_pc = numeric(100)

pars_theta_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3])
thetavec_pc = seq(0.01, 10, length = 100)
thetaprof_pc = numeric(100)

# Optimize original function without fixed parameters
pars_2_pc <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_pc <- optim(par = pars_2_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n,
                           db1 = dat_6$season >= 2, 
                           db2 = dat_6$season >= 3, 
                           db3 = dat_6$season >= 4, 
                           db4 = dat_6$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_pc <- list(rep(list(), times = 100))
x100_coefs_vec_pc <- list(rep(list(), times = 100))
c_coefs_vec_pc <- list(rep(list(), times = 100))
theta_coefs_vec_pc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_pc <- optim(par = pars_b_pc, fn = ci_b_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, b = bvec_pc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_6$season >= 2, 
                    db2 = dat_6$season >= 3, 
                    db3 = dat_6$season >= 4, 
                    db4 = dat_6$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_pc[[j]] = opt_b_pc$par # Store parameters for control
  # print(opt_a$convergence)
  bprof_pc[j] <- opt_b_pc$value # Set likelihood value for profile
  
  opt_x100_pc <- optim(par = pars_x100_pc, fn = ci_x100_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, x100 = x100vec_pc[j],
                       db1 = dat_6$season >= 2, 
                       db2 = dat_6$season >= 3, 
                       db3 = dat_6$season >= 4, 
                       db4 = dat_6$season >= 5, 
                       
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_pc[j] <- opt_x100_pc$value
  x100_coefs_vec_pc[[j]] = opt_x100_pc$par
  
  opt_c_pc <- optim(par = pars_c_pc, fn = ci_c_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, c = cvec_pc[j],
                    db1 = dat_6$season >= 2, 
                    db2 = dat_6$season >= 3, 
                    db3 = dat_6$season >= 4, 
                    db4 = dat_6$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_pc[j] <- opt_c_pc$value
  c_coefs_vec_pc[[j]] = opt_c_pc$par
  
  opt_theta_pc <- optim(par = pars_theta_pc, fn = ci_theta_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, theta = thetavec_pc[j], 
                        db1 = dat_6$season >= 2, 
                        db2 = dat_6$season >= 3, 
                        db3 = dat_6$season >= 4, 
                        db4 = dat_6$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_pc[j] <- opt_theta_pc$value
  theta_coefs_vec_pc[[j]] = opt_theta_pc$par
}

# Set the 95% confidence limits
bprof_lower_pc <- bprof_pc[1:which.min(bprof_pc)] # Likelihood values below the best estimate
bvec_lower_pc <- bvec_pc[1:which.min(bprof_pc)] # Parameter values below the best estimate
bprof_higher_pc <- bprof_pc[which.min(bprof_pc):length(bprof_pc)] # Likelihood values above the best estimate
bvec_higher_pc <- bvec_pc[which.min(bprof_pc):length(bprof_pc)] # Parameter values above the best estimate
l_b_ci_pc <- approx(bprof_lower_pc, bvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_pc <- approx(bprof_higher_pc, bvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_pc <- x100prof_pc[1:which.min(x100prof_pc)]
x100vec_lower_pc <- x100vec_pc[1:which.min(x100prof_pc)]
x100prof_higher_pc <- x100prof_pc[which.min(x100prof_pc):length(x100prof_pc)]
x100vec_higher_pc <- x100vec_pc[which.min(x100prof_pc):length(x100prof_pc)]
l_x100_ci_pc <- approx(x100prof_lower_pc, x100vec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_x100_ci_pc <- approx(x100prof_higher_pc, x100vec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

cprof_lower_pc <- cprof_pc[1:which.min(cprof_pc)]
cvec_lower_pc <- cvec_pc[1:which.min(cprof_pc)]
cprof_higher_pc <- cprof_pc[which.min(cprof_pc):length(cprof_pc)]
cvec_higher_pc <- cvec_pc[which.min(cprof_pc):length(cprof_pc)]
l_c_ci_pc <- approx(cprof_lower_pc, cvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_c_ci_pc <- approx(cprof_higher_pc, cvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

thetaprof_lower_pc <- thetaprof_pc[1:which.min(thetaprof_pc)]
thetavec_lower_pc <- thetavec_pc[1:which.min(thetaprof_pc)]
thetaprof_higher_pc <- thetaprof_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
thetavec_higher_pc <- thetavec_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
l_theta_ci_pc <- approx(thetaprof_lower_pc, thetavec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)
r_theta_ci_pc <- approx(thetaprof_higher_pc, thetavec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_pc_cne_sn <- c(opt_log_bbinom_pc$par[1], l_b_ci_pc$y, r_b_ci_pc$y) 
x100_ci_pc_cne_sn <- c(opt_log_bbinom_pc$par[2], l_x100_ci_pc$y, r_x100_ci_pc$y)
c_ci_pc_cne_sn <- c(opt_log_bbinom_pc$par[3], l_c_ci_pc$y, r_c_ci_pc$y)
theta_ci_pc_cne_sn <- c(opt_log_bbinom_pc$par[4], l_theta_ci_pc$y, r_theta_ci_pc$y)

b_ci_pc_cne_sn
x100_ci_pc_cne_sn
c_ci_pc_cne_sn
theta_ci_pc_cne_sn


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function
mlsp_fun_pc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_pc <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_pc <- optim(par = mlsp_par_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, 
                     db1 = dat_6$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_6$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_6$season >= 4, 
                     db4 = dat_6$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_pc <- mlsp_opt_pc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_6[dat_6$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_6[dat_6$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(0)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(0))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(0))^mlsp_coefs_pc[2])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(1)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(1))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(1))^mlsp_coefs_pc[2])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(2)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(2))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(2))^mlsp_coefs_pc[2])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(3)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(3))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(3))^mlsp_coefs_pc[2])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_pc[4]*(4)), mlsp_coefs_pc[1], 
                                        (mlsp_coefs_pc[1]*(1 - ((x - mlsp_coefs_pc[4]*(4))^mlsp_coefs_pc[2])/
                                                   (mlsp_coefs_pc[3]^mlsp_coefs_pc[2] + (x - mlsp_coefs_pc[4]*(4))^mlsp_coefs_pc[2])))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C25. Hill function through the positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_pc <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_pc <- function(par, x, z, n, pn, db1, db2, db3, db4){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_pc <- function(par, x, z, n, h, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_pc <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_pc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_pc = c(mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4], mlsp_coefs_pc[5]) # Starting parameters are the estimated parameters of the optimized function
avec_pc = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_pc = numeric(100) # Prepare profile of likelihood values

pars_pn_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[3], mlsp_coefs_pc[4], mlsp_coefs_pc[5])
pnvec_pc = seq(0.1, 10, length = 100)
pnprof_pc = numeric(100)

pars_h_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[4], mlsp_coefs_pc[5])
hvec_pc = seq(2, 100, length = 100)
hprof_pc = numeric(100)

pars_c_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[5])
cvec_pc = seq(0.1, 50, length = 100)
cprof_pc = numeric(100)

pars_theta_pc = c(mlsp_coefs_pc[1], mlsp_coefs_pc[2], mlsp_coefs_pc[3], mlsp_coefs_pc[4])
thetavec_pc = seq(0.01, 10, length = 100)
thetaprof_pc = numeric(100)

# Optimize original function without fixed parameters
pars_2_pc <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_pc <- optim(par = pars_2_pc, fn = mlsp_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n,
                           db1 = dat_6$season >= 2, 
                           db2 = dat_6$season >= 3, 
                           db3 = dat_6$season >= 4, 
                           db4 = dat_6$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_pc <- list(rep(list(), times = 100))
pn_coefs_vec_pc <- list(rep(list(), times = 100))
h_coefs_vec_pc <- list(rep(list(), times = 100))
c_coefs_vec_pc <- list(rep(list(), times = 100))
theta_coefs_vec_pc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_pc <- optim(par = pars_a_pc, fn = ci_a_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, a = avec_pc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_6$season >= 2, 
                    db2 = dat_6$season >= 3, 
                    db3 = dat_6$season >= 4, 
                    db4 = dat_6$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_pc[[j]] = opt_a_pc$par # Store parameters for control
  aprof_pc[j] <- opt_a_pc$value # Set likelihood value for profile
  
  opt_pn_pc <- optim(par = pars_pn_pc, fn = ci_pn_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, pn = pnvec_pc[j],
                     db1 = dat_6$season >= 2, 
                     db2 = dat_6$season >= 3, 
                     db3 = dat_6$season >= 4, 
                     db4 = dat_6$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_pc[j] <- opt_pn_pc$value
  pn_coefs_vec_pc[[j]] = opt_pn_pc$par
  
  opt_h_pc <- optim(par = pars_h_pc, fn = ci_h_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, h = hvec_pc[j],
                    db1 = dat_6$season >= 2, 
                    db2 = dat_6$season >= 3, 
                    db3 = dat_6$season >= 4, 
                    db4 = dat_6$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_pc[j] <- opt_h_pc$value
  h_coefs_vec_pc[[j]] = opt_h_pc$par
  
  opt_c_pc <- optim(par = pars_c_pc, fn = ci_c_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, c = cvec_pc[j],
                    db1 = dat_6$season >= 2, 
                    db2 = dat_6$season >= 3, 
                    db3 = dat_6$season >= 4, 
                    db4 = dat_6$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_pc[j] <- opt_c_pc$value
  c_coefs_vec_pc[[j]] = opt_c_pc$par
  
  opt_theta_pc <- optim(par = pars_theta_pc, fn = ci_theta_fun_pc, x = dat_6$dist, z = dat_6$pos, n = dat_6$n, theta = thetavec_pc[j], 
                        db1 = dat_6$season >= 2, 
                        db2 = dat_6$season >= 3, 
                        db3 = dat_6$season >= 4, 
                        db4 = dat_6$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_pc[j] <- opt_theta_pc$value
  theta_coefs_vec_pc[[j]] = opt_theta_pc$par
}

# Set the 95% confidence limits
aprof_lower_pc <- aprof_pc[1:which.min(aprof_pc)] # Likelihood values below the best estimate
avec_lower_pc <- avec_pc[1:which.min(aprof_pc)] # Parameter values below the best estimate
aprof_higher_pc <- aprof_pc[which.min(aprof_pc):length(aprof_pc)] # Likelihood values above the best estimate
avec_higher_pc <- avec_pc[which.min(aprof_pc):length(aprof_pc)] # Parameter values above the best estimate
l_a_ci_pc <- approx(aprof_lower_pc, avec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_pc <- approx(aprof_higher_pc, avec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_pc <- pnprof_pc[1:which.min(pnprof_pc)]
pnvec_lower_pc <- pnvec_pc[1:which.min(pnprof_pc)]
pnprof_higher_pc <- pnprof_pc[which.min(pnprof_pc):length(pnprof_pc)]
pnvec_higher_pc <- pnvec_pc[which.min(pnprof_pc):length(pnprof_pc)]
l_pn_ci_pc <- approx(pnprof_lower_pc, pnvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_pn_ci_pc <- approx(pnprof_higher_pc, pnvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

hprof_lower_pc <- hprof_pc[1:which.min(hprof_pc)]
hvec_lower_pc <- hvec_pc[1:which.min(hprof_pc)]
hprof_higher_pc <- hprof_pc[which.min(hprof_pc):length(hprof_pc)]
hvec_higher_pc <- hvec_pc[which.min(hprof_pc):length(hprof_pc)]
l_h_ci_pc <- approx(hprof_lower_pc, hvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_h_ci_pc <- approx(hprof_higher_pc, hvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

cprof_lower_pc <- cprof_pc[1:which.min(cprof_pc)]
cvec_lower_pc <- cvec_pc[1:which.min(cprof_pc)]
cprof_higher_pc <- cprof_pc[which.min(cprof_pc):length(cprof_pc)]
cvec_higher_pc <- cvec_pc[which.min(cprof_pc):length(cprof_pc)]
l_c_ci_pc <- approx(cprof_lower_pc, cvec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_c_ci_pc <- approx(cprof_higher_pc, cvec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

thetaprof_lower_pc <- thetaprof_pc[1:which.min(thetaprof_pc)]
thetavec_lower_pc <- thetavec_pc[1:which.min(thetaprof_pc)]
thetaprof_higher_pc <- thetaprof_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
thetavec_higher_pc <- thetavec_pc[which.min(thetaprof_pc):length(thetaprof_pc)]
l_theta_ci_pc <- approx(thetaprof_lower_pc, thetavec_lower_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)
r_theta_ci_pc <- approx(thetaprof_higher_pc, thetavec_higher_pc, xout = opt_log_bbinom_pc$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_pc_hill_sn <- c(opt_log_bbinom_pc$par[1], l_a_ci_pc$y, r_a_ci_pc$y) 
pn_ci_pc_hill_sn <- c(opt_log_bbinom_pc$par[2], l_pn_ci_pc$y, r_pn_ci_pc$y)
h_ci_pc_hill_sn <- c(opt_log_bbinom_pc$par[3], l_h_ci_pc$y, r_h_ci_pc$y)
c_ci_pc_hill_sn <- c(opt_log_bbinom_pc$par[4], l_c_ci_pc$y, r_c_ci_pc$y)
theta_ci_pc_hill_sn <- c(opt_log_bbinom_pc$par[5], l_theta_ci_pc$y, r_theta_ci_pc$y)

a_ci_pc_hill_sn
pn_ci_pc_hill_sn
h_ci_pc_hill_sn
c_ci_pc_hill_sn
theta_ci_pc_hill_sn


##############################
# Negative clairvoyance data #
##############################

#
# Data morphing and distance circle creation
#

dat_7 <- dat

# Make data clairvoyant
for(i in length(seasons):2){ # Work backwards, from 2018 to 2014
  co <- which(dat_7[[i]]$result == 0) # Which of the data have a result of 0
  dat_7[[i-1]] <- rbind(dat_7[[i-1]], dat_7[[i]][co,]) # Combine the previous year with the negative results of the later year
  co_2 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & dat_7[[i-1]]$season == (2012 + i)) # Check for multiple samples at the same location
  if(length(co_2) >= 1){ # If there are multiple samples at one location
    for(j in 1:length(co_2)){ # For every duplicate sample
      if(dat_7[[i-1]]$result[co_2[j]] == 1){ # If the duplicate from the year before has result 1, the result at this location should be a 1
        co_3 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & dat_7[[i-1]]$season == (2012 + i))
        dat_7[[i-1]] <- dat_7[[i-1]][-co_3[j],]
      }
      if(dat_7[[i-1]]$result[co_2[j]] == 0){ # If the result is a 0, it will stay a 0.
        dat_7[[i-1]] <- dat_7[[i-1]][-co_2[j],]
      }
    }
  } # Note that this should never be the case, since a 1 would have turned into a 0. This part of code is a safeguard to be used in future data or simulated data.
}

# Set distance circles
dc <- 1:(max_dist) # Number of distance circles is the maximum distance a measurement is done.
pos <- rep(list(rep(NA, times = length(dc))), times = length(seasons)) # Prepare list for number of positives in each dc
n <- rep(list(rep(NA, times = length(dc))), times = length(seasons)) # Prepare list for total number of measurements in each dc

dat_8 <- list() # Prepare data list

for(i in 1:length(seasons)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j] <- sum((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j] <- length((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_8[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
}

# Combine seasons and create year column for year 
dat_9 <- list()
for(i in 1:length(seasons)){
  dat_8[[i]] <- dat_8[[i]] %>% 
    mutate(season = (i))
  dat_9 <- rbind(dat_9, dat_8[[i]])
}

# Plot the data
ggplot() +
  geom_point(data = dat_9[dat_9$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  labs(
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C26. Data points of the negative clairvoyance data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function
mlsp_fun_nc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- (1/(1 + exp(par[1] * (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3]))))) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) 
  return(nll)
}

# Setup parameter starting values
mlsp_par_nc <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_nc <- optim(par = mlsp_par_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, 
                     db1 = dat_9$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_9$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_9$season >= 4, 
                     db4 = dat_9$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_nc <- mlsp_opt_nc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_9[dat_9$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 0*mlsp_coefs_nc[3]))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 1*mlsp_coefs_nc[3]))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 2*mlsp_coefs_nc[3]))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 3*mlsp_coefs_nc[3]))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_nc[1]*(x - (mlsp_coefs_nc[2] + 4*mlsp_coefs_nc[3]))))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C27. Logistic function through the negative clairvoyance data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_nc <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_nc <- function(par, x, z, n, x50, db1, db2, db3, db4){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_nc <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_nc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_nc = c(mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4]) # Starting parameters are the estimated parameters of the optimized function
avec_nc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_nc = numeric(100) # Prepare profile of likelihood values

pars_x50_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[3], mlsp_coefs_nc[4])
x50vec_nc = seq(-120, 100, length = 100)
x50prof_nc = numeric(100)

pars_c_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[4])
cvec_nc = seq(0.1, 50, length = 100)
cprof_nc = numeric(100)

pars_theta_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3])
thetavec_nc = seq(0.01, 10, length = 100)
thetaprof_nc = numeric(100)

# Optimize original function without fixed parameters
pars_2_nc <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_nc <- optim(par = pars_2_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n,
                           db1 = dat_9$season >= 2, 
                           db2 = dat_9$season >= 3, 
                           db3 = dat_9$season >= 4, 
                           db4 = dat_9$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_nc <- list(rep(list(), times = 100))
x50_coefs_vec_nc <- list(rep(list(), times = 100))
c_coefs_vec_nc <- list(rep(list(), times = 100))
theta_coefs_vec_nc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_nc <- optim(par = pars_a_nc, fn = ci_a_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, a = avec_nc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_9$season >= 2, 
                    db2 = dat_9$season >= 3, 
                    db3 = dat_9$season >= 4, 
                    db4 = dat_9$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_nc[[j]] = opt_a_nc$par # Store parameters for control
  # print(opt_a$convergence)
  aprof_nc[j] <- opt_a_nc$value # Set likelihood value for profile
  
  opt_x50_nc <- optim(par = pars_x50_nc, fn = ci_x50_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, x50 = x50vec_nc[j],
                      db1 = dat_9$season >= 2, 
                      db2 = dat_9$season >= 3, 
                      db3 = dat_9$season >= 4, 
                      db4 = dat_9$season >= 5, 
                      
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_nc[j] <- opt_x50_nc$value
  x50_coefs_vec_nc[[j]] = opt_x50_nc$par
  
  opt_c_nc <- optim(par = pars_c_nc, fn = ci_c_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, c = cvec_nc[j],
                    db1 = dat_9$season >= 2, 
                    db2 = dat_9$season >= 3, 
                    db3 = dat_9$season >= 4, 
                    db4 = dat_9$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_nc[j] <- opt_c_nc$value
  c_coefs_vec_nc[[j]] = opt_c_nc$par
  
  opt_theta_nc <- optim(par = pars_theta_nc, fn = ci_theta_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, theta = thetavec_nc[j], 
                        db1 = dat_9$season >= 2, 
                        db2 = dat_9$season >= 3, 
                        db3 = dat_9$season >= 4, 
                        db4 = dat_9$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_nc[j] <- opt_theta_nc$value
  theta_coefs_vec_nc[[j]] = opt_theta_nc$par
}

# Set the 95% confidence limits
aprof_lower_nc <- aprof_nc[1:which.min(aprof_nc)] # Likelihood values below the best estimate
avec_lower_nc <- avec_nc[1:which.min(aprof_nc)] # Parameter values below the best estimate
aprof_higher_nc <- aprof_nc[which.min(aprof_nc):length(aprof_nc)] # Likelihood values above the best estimate
avec_higher_nc <- avec_nc[which.min(aprof_nc):length(aprof_nc)] # Parameter values above the best estimate
l_a_ci_nc <- approx(aprof_lower_nc, avec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_nc <- approx(aprof_higher_nc, avec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x50prof_lower_nc <- x50prof_nc[1:which.min(x50prof_nc)]
x50vec_lower_nc <- x50vec_nc[1:which.min(x50prof_nc)]
x50prof_higher_nc <- x50prof_nc[which.min(x50prof_nc):length(x50prof_nc)]
x50vec_higher_nc <- x50vec_nc[which.min(x50prof_nc):length(x50prof_nc)]
l_x50_ci_nc <- approx(x50prof_lower_nc, x50vec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_x50_ci_nc <- approx(x50prof_higher_nc, x50vec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

cprof_lower_nc <- cprof_nc[1:which.min(cprof_nc)]
cvec_lower_nc <- cvec_nc[1:which.min(cprof_nc)]
cprof_higher_nc <- cprof_nc[which.min(cprof_nc):length(cprof_nc)]
cvec_higher_nc <- cvec_nc[which.min(cprof_nc):length(cprof_nc)]
l_c_ci_nc <- approx(cprof_lower_nc, cvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_c_ci_nc <- approx(cprof_higher_nc, cvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

thetaprof_lower_nc <- thetaprof_nc[1:which.min(thetaprof_nc)]
thetavec_lower_nc <- thetavec_nc[1:which.min(thetaprof_nc)]
thetaprof_higher_nc <- thetaprof_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
thetavec_higher_nc <- thetavec_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
l_theta_ci_nc <- approx(thetaprof_lower_nc, thetavec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_theta_ci_nc <- approx(thetaprof_higher_nc, thetavec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_nc_log_sn <- c(opt_log_bbinom_nc$par[1], l_a_ci_nc$y, r_a_ci_nc$y) 
x50_ci_nc_log_sn <- c(opt_log_bbinom_nc$par[2], l_x50_ci_nc$y, r_x50_ci_nc$y)
c_ci_nc_log_sn <- c(opt_log_bbinom_nc$par[3], l_c_ci_nc$y, r_c_ci_nc$y)
theta_ci_nc_log_sn <- c(opt_log_bbinom_nc$par[4], l_theta_ci_nc$y, r_theta_ci_nc$y)

a_ci_nc_log_sn
x50_ci_nc_log_sn
c_ci_nc_log_sn
theta_ci_nc_log_sn


##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function
mlsp_fun_nc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_nc <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_nc <- optim(par = mlsp_par_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, 
                     db1 = dat_9$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_9$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_9$season >= 4, 
                     db4 = dat_9$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_nc <- mlsp_opt_nc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_9[dat_9$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((0) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((0) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((1) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((1) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((2) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((2) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((3) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((3) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((4) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_nc[1]*(mlsp_coefs_nc[2] + ((4) * mlsp_coefs_nc[3])))) * exp(-mlsp_coefs_nc[1] * x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C28. CNE function through the negative clairvoyance data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_nc <- function(par, x, z, n, b, db1, db2, db3, db4){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_nc <- function(par, x, z, n, x100, db1, db2, db3, db4){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_nc <- function(par, x, z, n, c, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_nc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_nc = c(mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_nc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_nc = numeric(100) # Prepare profile of likelihood values

pars_x100_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[3], mlsp_coefs_nc[4])
x100vec_nc = seq(-120, 100, length = 100)
x100prof_nc = numeric(100)

pars_c_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[4])
cvec_nc = seq(0.1, 50, length = 100)
cprof_nc = numeric(100)

pars_theta_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3])
thetavec_nc = seq(0.01, 10, length = 100)
thetaprof_nc = numeric(100)

# Optimize original function without fixed parameters
pars_2_nc <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_nc <- optim(par = pars_2_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n,
                           db1 = dat_9$season >= 2, 
                           db2 = dat_9$season >= 3, 
                           db3 = dat_9$season >= 4, 
                           db4 = dat_9$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_nc <- list(rep(list(), times = 100))
x100_coefs_vec_nc <- list(rep(list(), times = 100))
c_coefs_vec_nc <- list(rep(list(), times = 100))
theta_coefs_vec_nc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_nc <- optim(par = pars_b_nc, fn = ci_b_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, b = bvec_nc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_9$season >= 2, 
                    db2 = dat_9$season >= 3, 
                    db3 = dat_9$season >= 4, 
                    db4 = dat_9$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_nc[[j]] = opt_b_nc$par # Store parameters for control
  bprof_nc[j] <- opt_b_nc$value # Set likelihood value for profile
  
  opt_x100_nc <- optim(par = pars_x100_nc, fn = ci_x100_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, x100 = x100vec_nc[j],
                       db1 = dat_9$season >= 2, 
                       db2 = dat_9$season >= 3, 
                       db3 = dat_9$season >= 4, 
                       db4 = dat_9$season >= 5, 
                       
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_nc[j] <- opt_x100_nc$value
  x100_coefs_vec_nc[[j]] = opt_x100_nc$par
  
  opt_c_nc <- optim(par = pars_c_nc, fn = ci_c_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, c = cvec_nc[j],
                    db1 = dat_9$season >= 2, 
                    db2 = dat_9$season >= 3, 
                    db3 = dat_9$season >= 4, 
                    db4 = dat_9$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_nc[j] <- opt_c_nc$value
  c_coefs_vec_nc[[j]] = opt_c_nc$par
  
  opt_theta_nc <- optim(par = pars_theta_nc, fn = ci_theta_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, theta = thetavec_nc[j], 
                        db1 = dat_9$season >= 2, 
                        db2 = dat_9$season >= 3, 
                        db3 = dat_9$season >= 4, 
                        db4 = dat_9$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_nc[j] <- opt_theta_nc$value
  theta_coefs_vec_nc[[j]] = opt_theta_nc$par
}

# Set the 95% confidence limits
bprof_lower_nc <- bprof_nc[1:which.min(bprof_nc)] # Likelihood values below the best estimate
bvec_lower_nc <- bvec_nc[1:which.min(bprof_nc)] # Parameter values below the best estimate
bprof_higher_nc <- bprof_nc[which.min(bprof_nc):length(bprof_nc)] # Likelihood values above the best estimate
bvec_higher_nc <- bvec_nc[which.min(bprof_nc):length(bprof_nc)] # Parameter values above the best estimate
l_b_ci_nc <- approx(bprof_lower_nc, bvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_nc <- approx(bprof_higher_nc, bvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_nc <- x100prof_nc[1:which.min(x100prof_nc)]
x100vec_lower_nc <- x100vec_nc[1:which.min(x100prof_nc)]
x100prof_higher_nc <- x100prof_nc[which.min(x100prof_nc):length(x100prof_nc)]
x100vec_higher_nc <- x100vec_nc[which.min(x100prof_nc):length(x100prof_nc)]
l_x100_ci_nc <- approx(x100prof_lower_nc, x100vec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_x100_ci_nc <- approx(x100prof_higher_nc, x100vec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

cprof_lower_nc <- cprof_nc[1:which.min(cprof_nc)]
cvec_lower_nc <- cvec_nc[1:which.min(cprof_nc)]
cprof_higher_nc <- cprof_nc[which.min(cprof_nc):length(cprof_nc)]
cvec_higher_nc <- cvec_nc[which.min(cprof_nc):length(cprof_nc)]
l_c_ci_nc <- approx(cprof_lower_nc, cvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_c_ci_nc <- approx(cprof_higher_nc, cvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

thetaprof_lower_nc <- thetaprof_nc[1:which.min(thetaprof_nc)]
thetavec_lower_nc <- thetavec_nc[1:which.min(thetaprof_nc)]
thetaprof_higher_nc <- thetaprof_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
thetavec_higher_nc <- thetavec_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
l_theta_ci_nc <- approx(thetaprof_lower_nc, thetavec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)
r_theta_ci_nc <- approx(thetaprof_higher_nc, thetavec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_nc_cne_sn <- c(opt_log_bbinom_nc$par[1], l_b_ci_nc$y, r_b_ci_nc$y) 
x100_ci_nc_cne_sn <- c(opt_log_bbinom_nc$par[2], l_x100_ci_nc$y, r_x100_ci_nc$y)
c_ci_nc_cne_sn <- c(opt_log_bbinom_nc$par[3], l_c_ci_nc$y, r_c_ci_nc$y)
theta_ci_nc_cne_sn <- c(opt_log_bbinom_nc$par[4], l_theta_ci_nc$y, r_theta_ci_nc$y)

b_ci_nc_cne_sn
x100_ci_nc_cne_sn
c_ci_nc_cne_sn
theta_ci_nc_cne_sn


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function
mlsp_fun_nc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_nc <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_nc <- optim(par = mlsp_par_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, 
                     db1 = dat_9$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_9$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_9$season >= 4, 
                     db4 = dat_9$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_nc <- mlsp_opt_nc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_9[dat_9$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_9[dat_9$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(0)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(0))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(0))^mlsp_coefs_nc[2])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(1)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(1))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(1))^mlsp_coefs_nc[2])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(2)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(2))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(2))^mlsp_coefs_nc[2])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(3)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(3))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(3))^mlsp_coefs_nc[2])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_nc[4]*(4)), mlsp_coefs_nc[1], 
                                        (mlsp_coefs_nc[1]*(1 - ((x - mlsp_coefs_nc[4]*(4))^mlsp_coefs_nc[2])/
                                                   (mlsp_coefs_nc[3]^mlsp_coefs_nc[2] + (x - mlsp_coefs_nc[4]*(4))^mlsp_coefs_nc[2])))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C29. Hill function through the negative clairvoyance data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_nc <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_nc <- function(par, x, z, n, pn, db1, db2, db3, db4){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_nc <- function(par, x, z, n, h, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_nc <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_nc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_nc = c(mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4], mlsp_coefs_nc[5]) # Starting parameters are the estimated parameters of the optimized function
avec_nc = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_nc = numeric(100) # Prepare profile of likelihood values

pars_pn_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[3], mlsp_coefs_nc[4], mlsp_coefs_nc[5])
pnvec_nc = seq(0.1, 10, length = 100)
pnprof_nc = numeric(100)

pars_h_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[4], mlsp_coefs_nc[5])
hvec_nc = seq(4, 100, length = 100)
hprof_nc = numeric(100)

pars_c_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[5])
cvec_nc = seq(0.1, 50, length = 100)
cprof_nc = numeric(100)

pars_theta_nc = c(mlsp_coefs_nc[1], mlsp_coefs_nc[2], mlsp_coefs_nc[3], mlsp_coefs_nc[4])
thetavec_nc = seq(0.01, 10, length = 100)
thetaprof_nc = numeric(100)

# Optimize original function without fixed parameters
pars_2_nc <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_nc <- optim(par = pars_2_nc, fn = mlsp_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n,
                           db1 = dat_9$season >= 2, 
                           db2 = dat_9$season >= 3, 
                           db3 = dat_9$season >= 4, 
                           db4 = dat_9$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_nc <- list(rep(list(), times = 100))
pn_coefs_vec_nc <- list(rep(list(), times = 100))
h_coefs_vec_nc <- list(rep(list(), times = 100))
c_coefs_vec_nc <- list(rep(list(), times = 100))
theta_coefs_vec_nc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_nc <- optim(par = pars_a_nc, fn = ci_a_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, a = avec_nc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_9$season >= 2, 
                    db2 = dat_9$season >= 3, 
                    db3 = dat_9$season >= 4, 
                    db4 = dat_9$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_nc[[j]] = opt_a_nc$par # Store parameters for control
  aprof_nc[j] <- opt_a_nc$value # Set likelihood value for profile
  
  opt_pn_nc <- optim(par = pars_pn_nc, fn = ci_pn_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, pn = pnvec_nc[j],
                     db1 = dat_9$season >= 2, 
                     db2 = dat_9$season >= 3, 
                     db3 = dat_9$season >= 4, 
                     db4 = dat_9$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_nc[j] <- opt_pn_nc$value
  pn_coefs_vec_nc[[j]] = opt_pn_nc$par
  
  opt_h_nc <- optim(par = pars_h_nc, fn = ci_h_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, h = hvec_nc[j],
                    db1 = dat_9$season >= 2, 
                    db2 = dat_9$season >= 3, 
                    db3 = dat_9$season >= 4, 
                    db4 = dat_9$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_nc[j] <- opt_h_nc$value
  h_coefs_vec_nc[[j]] = opt_h_nc$par
  
  opt_c_nc <- optim(par = pars_c_nc, fn = ci_c_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, c = cvec_nc[j],
                    db1 = dat_9$season >= 2, 
                    db2 = dat_9$season >= 3, 
                    db3 = dat_9$season >= 4, 
                    db4 = dat_9$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_nc[j] <- opt_c_nc$value
  c_coefs_vec_nc[[j]] = opt_c_nc$par
  
  opt_theta_nc <- optim(par = pars_theta_nc, fn = ci_theta_fun_nc, x = dat_9$dist, z = dat_9$pos, n = dat_9$n, theta = thetavec_nc[j], 
                        db1 = dat_9$season >= 2, 
                        db2 = dat_9$season >= 3, 
                        db3 = dat_9$season >= 4, 
                        db4 = dat_9$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_nc[j] <- opt_theta_nc$value
  theta_coefs_vec_nc[[j]] = opt_theta_nc$par
}

# Set the 95% confidence limits
aprof_lower_nc <- aprof_nc[1:which.min(aprof_nc)] # Likelihood values below the best estimate
avec_lower_nc <- avec_nc[1:which.min(aprof_nc)] # Parameter values below the best estimate
aprof_higher_nc <- aprof_nc[which.min(aprof_nc):length(aprof_nc)] # Likelihood values above the best estimate
avec_higher_nc <- avec_nc[which.min(aprof_nc):length(aprof_nc)] # Parameter values above the best estimate
l_a_ci_nc <- approx(aprof_lower_nc, avec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_nc <- approx(aprof_higher_nc, avec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_nc <- pnprof_nc[1:which.min(pnprof_nc)]
pnvec_lower_nc <- pnvec_nc[1:which.min(pnprof_nc)]
pnprof_higher_nc <- pnprof_nc[which.min(pnprof_nc):length(pnprof_nc)]
pnvec_higher_nc <- pnvec_nc[which.min(pnprof_nc):length(pnprof_nc)]
l_pn_ci_nc <- approx(pnprof_lower_nc, pnvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_pn_ci_nc <- approx(pnprof_higher_nc, pnvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

hprof_lower_nc <- hprof_nc[1:which.min(hprof_nc)]
hvec_lower_nc <- hvec_nc[1:which.min(hprof_nc)]
hprof_higher_nc <- hprof_nc[which.min(hprof_nc):length(hprof_nc)]
hvec_higher_nc <- hvec_nc[which.min(hprof_nc):length(hprof_nc)]
l_h_ci_nc <- approx(hprof_lower_nc, hvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_h_ci_nc <- approx(hprof_higher_nc, hvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

cprof_lower_nc <- cprof_nc[1:which.min(cprof_nc)]
cvec_lower_nc <- cvec_nc[1:which.min(cprof_nc)]
cprof_higher_nc <- cprof_nc[which.min(cprof_nc):length(cprof_nc)]
cvec_higher_nc <- cvec_nc[which.min(cprof_nc):length(cprof_nc)]
l_c_ci_nc <- approx(cprof_lower_nc, cvec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_c_ci_nc <- approx(cprof_higher_nc, cvec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

thetaprof_lower_nc <- thetaprof_nc[1:which.min(thetaprof_nc)]
thetavec_lower_nc <- thetavec_nc[1:which.min(thetaprof_nc)]
thetaprof_higher_nc <- thetaprof_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
thetavec_higher_nc <- thetavec_nc[which.min(thetaprof_nc):length(thetaprof_nc)]
l_theta_ci_nc <- approx(thetaprof_lower_nc, thetavec_lower_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)
r_theta_ci_nc <- approx(thetaprof_higher_nc, thetavec_higher_nc, xout = opt_log_bbinom_nc$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_nc_hill_sn <- c(opt_log_bbinom_nc$par[1], l_a_ci_nc$y, r_a_ci_nc$y) 
pn_ci_nc_hill_sn <- c(opt_log_bbinom_nc$par[2], l_pn_ci_nc$y, r_pn_ci_nc$y)
h_ci_nc_hill_sn <- c(opt_log_bbinom_nc$par[3], l_h_ci_nc$y, r_h_ci_nc$y)
c_ci_nc_hill_sn <- c(opt_log_bbinom_nc$par[4], l_c_ci_nc$y, r_c_ci_nc$y)
theta_ci_nc_hill_sn <- c(opt_log_bbinom_nc$par[5], l_theta_ci_nc$y, r_theta_ci_nc$y)

a_ci_nc_hill_sn
pn_ci_nc_hill_sn
h_ci_nc_hill_sn
c_ci_nc_hill_sn
theta_ci_nc_hill_sn


####################################################
# Negative clairvoyance & positive carry-over data #
####################################################

#
# Data morphing and distance circle creation
#

dat_10 <- dat

# Make data clairvoyant
for(i in length(seasons):2){ # Work backwards, from 2018 to 2014
  co <- which(dat_10[[i]]$result == 0) # Which of the data have a result of 0
  dat_10[[i-1]] <- rbind(dat_10[[i-1]], dat_10[[i]][co,]) # Combine the previous year with the negative results of the later year
  co_2 <- which(dat_10[[i-1]]$lon == dat_10[[i]]$lon[co] & dat_10[[i-1]]$lat == dat_10[[i]]$lat[co] & dat_10[[i-1]]$season == (2012 + i)) # Check for multiple samples at the same location
  if(length(co_2) >= 1){ # If there are multiple samples at one location
    for(j in 1:length(co_2)){ # For every duplicate sample
      if(dat_10[[i-1]]$result[co_2[j]] == 1){ # If the duplicate from the year before has result 1, the result at this location should be a 1
        co_3 <- which(dat_10[[i-1]]$lon == dat_10[[i]]$lon[co] & dat_10[[i-1]]$lat == dat_10[[i]]$lat[co] & dat_10[[i-1]]$season == (2012 + i))
        dat_10[[i-1]] <- dat_10[[i-1]][-co_3[j],]
      }
      if(dat_10[[i-1]]$result[co_2[j]] == 0){ # If the result is a 0, it will stay a 0.
        dat_10[[i-1]] <- dat_10[[i-1]][-co_2[j],]
      }
    }
  } # Note that this should never be the case, since a 1 would have turned into a 0. This part of code is a safeguard to be used in future data or simulated data.
}

# make positives carry over
for(i in 2:length(seasons)){ # For-loop through every year from 2014
  co <- which(dat_10[[i-1]]$result == 1) # Which lines of the previous year have a result of 1
  dat_10[[i]] <- rbind(dat_10[[i]], dat_10[[i-1]][co,]) # Include the lines of the previous year with a result of 1 with the current year
  co_2 <- which(dat_10[[i]]$lon == dat_10[[i-1]]$lon[co] & dat_10[[i]]$lat == dat_10[[i-1]]$lat[co] & dat_10[[i]]$season == (2012 + i)) # Which lines of this year have the same coordinates of the positives last year
  # print(co_2)
  if(length(co_2) >= 1){ # If this is not 0, then there has been a measurement at these coordinates this year as well as last year, or multiple times in a single year (with these data, only multiple measurements in 2013)
    dat_10[[i]] <- dat_10[[i]][-co_2,] # Remove all the duplicates. If the duplicate measurement from this year is a 0, it will still be a 1, because I assume that once a disease is there, it will not leave. I assume that 0 to be a false negative.
  }
}

# Set distance circles
dc <- 1:(max_dist) # Number of distance circles is the maximum distance a measurement is done.
pos <- rep(list(rep(NA, times = length(dc))), times = length(seasons)) # Prepare list for number of positives in each dc
n <- rep(list(rep(NA, times = length(dc))), times = length(seasons)) # Prepare list for total number of measurements in each dc

dat_11 <- list() # Prepare data list

for(i in 1:length(seasons)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j] <- sum((dat_10[[i]]$result[dat_10[[i]]$dist >= j & (dat_10[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j] <- length((dat_10[[i]]$result[dat_10[[i]]$dist >= j & (dat_10[[i]]$dist < j + 1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_11[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
}

# Combine seasons and create year column for year 
dat_12 <- list()
for(i in 1:length(seasons)){
  dat_11[[i]] <- dat_11[[i]] %>% 
    mutate(season = (i))
  dat_12 <- rbind(dat_12, dat_11[[i]])
}

# Plot the data
ggplot() +
  geom_point(data = dat_12[dat_12$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  labs(
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure C30. Data points of the negative clairvoyance & positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##                   ##
## Logistic function ##
##                   ##

#
# Model setup
#

# Setup logistic function
mlsp_fun_ncpc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- (1/(1 + exp(par[1] * (x - (par[2] + (1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3]))))) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_ncpc <- c(a = 0.05, x50 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_ncpc <- optim(par = mlsp_par_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, 
                     db1 = dat_12$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_12$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_12$season >= 4, 
                     db4 = dat_12$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_ncpc <- mlsp_opt_ncpc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_12[dat_12$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 0*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 1*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 2*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 3*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1/(1+exp(mlsp_coefs_ncpc[1]*(x - (mlsp_coefs_ncpc[2] + 4*mlsp_coefs_ncpc[3]))))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive", 
       caption = "Figure C31. Logistic function through the negative clairvoyance & positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_ncpc <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  x50 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x50_fun_ncpc <- function(par, x, z, n, x50, db1, db2, db3, db4){
  a = par[1]
  x50 = x50
  c = par[2]
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_ncpc <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = c
  theta = par[3]
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_ncpc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  x50 = par[2]
  c = par[3]
  theta = theta
  mu <- (1/(1 + exp(a *
                      (x - (x50 + (1*db1 + 1*db2 + 1*db3 + 1*db4 )*c)))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_ncpc = c(mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4]) # Starting parameters are the estimated parameters of the optimized function
avec_ncpc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
aprof_ncpc = numeric(100) # Prepare profile of likelihood values

pars_x50_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4])
x50vec_ncpc = seq(-120, 100, length = 100)
x50prof_ncpc = numeric(100)

pars_c_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[4])
cvec_ncpc = seq(0.1, 50, length = 100)
cprof_ncpc = numeric(100)

pars_theta_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3])
thetavec_ncpc = seq(0.01, 10, length = 100)
thetaprof_ncpc = numeric(100)

# Optimize original function without fixed parameters
pars_2_ncpc <- c(a = 0.05, x50 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_ncpc <- optim(par = pars_2_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n,
                           db1 = dat_12$season >= 2, 
                           db2 = dat_12$season >= 3, 
                           db3 = dat_12$season >= 4, 
                           db4 = dat_12$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_ncpc <- list(rep(list(), times = 100))
x50_coefs_vec_ncpc <- list(rep(list(), times = 100))
c_coefs_vec_ncpc <- list(rep(list(), times = 100))
theta_coefs_vec_ncpc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_ncpc <- optim(par = pars_a_ncpc, fn = ci_a_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, a = avec_ncpc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_12$season >= 2, 
                    db2 = dat_12$season >= 3, 
                    db3 = dat_12$season >= 4, 
                    db4 = dat_12$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_ncpc[[j]] = opt_a_ncpc$par # Store parameters for control
  # print(opt_a$convergence)
  aprof_ncpc[j] <- opt_a_ncpc$value # Set likelihood value for profile
  
  opt_x50_ncpc <- optim(par = pars_x50_ncpc, fn = ci_x50_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, x50 = x50vec_ncpc[j],
                      db1 = dat_12$season >= 2, 
                      db2 = dat_12$season >= 3, 
                      db3 = dat_12$season >= 4, 
                      db4 = dat_12$season >= 5, 
                      
                      control = list(maxit = 10000), method = "Nelder-Mead")
  x50prof_ncpc[j] <- opt_x50_ncpc$value
  x50_coefs_vec_ncpc[[j]] = opt_x50_ncpc$par
  
  opt_c_ncpc <- optim(par = pars_c_ncpc, fn = ci_c_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, c = cvec_ncpc[j],
                    db1 = dat_12$season >= 2, 
                    db2 = dat_12$season >= 3, 
                    db3 = dat_12$season >= 4, 
                    db4 = dat_12$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_ncpc[j] <- opt_c_ncpc$value
  c_coefs_vec_ncpc[[j]] = opt_c_ncpc$par
  
  opt_theta_ncpc <- optim(par = pars_theta_ncpc, fn = ci_theta_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, theta = thetavec_ncpc[j], 
                        db1 = dat_12$season >= 2, 
                        db2 = dat_12$season >= 3, 
                        db3 = dat_12$season >= 4, 
                        db4 = dat_12$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_ncpc[j] <- opt_theta_ncpc$value
  theta_coefs_vec_ncpc[[j]] = opt_theta_ncpc$par
}

# Set the 95% confidence limits
aprof_lower_ncpc <- aprof_ncpc[1:which.min(aprof_ncpc)] # Likelihood values below the best estimate
avec_lower_ncpc <- avec_ncpc[1:which.min(aprof_ncpc)] # Parameter values below the best estimate
aprof_higher_ncpc <- aprof_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Likelihood values above the best estimate
avec_higher_ncpc <- avec_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Parameter values above the best estimate
l_a_ci_ncpc <- approx(aprof_lower_ncpc, avec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_ncpc <- approx(aprof_higher_ncpc, avec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x50prof_lower_ncpc <- x50prof_ncpc[1:which.min(x50prof_ncpc)]
x50vec_lower_ncpc <- x50vec_ncpc[1:which.min(x50prof_ncpc)]
x50prof_higher_ncpc <- x50prof_ncpc[which.min(x50prof_ncpc):length(x50prof_ncpc)]
x50vec_higher_ncpc <- x50vec_ncpc[which.min(x50prof_ncpc):length(x50prof_ncpc)]
l_x50_ci_ncpc <- approx(x50prof_lower_ncpc, x50vec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_x50_ci_ncpc <- approx(x50prof_higher_ncpc, x50vec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

cprof_lower_ncpc <- cprof_ncpc[1:which.min(cprof_ncpc)]
cvec_lower_ncpc <- cvec_ncpc[1:which.min(cprof_ncpc)]
cprof_higher_ncpc <- cprof_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
cvec_higher_ncpc <- cvec_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
l_c_ci_ncpc <- approx(cprof_lower_ncpc, cvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_c_ci_ncpc <- approx(cprof_higher_ncpc, cvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

thetaprof_lower_ncpc <- thetaprof_ncpc[1:which.min(thetaprof_ncpc)]
thetavec_lower_ncpc <- thetavec_ncpc[1:which.min(thetaprof_ncpc)]
thetaprof_higher_ncpc <- thetaprof_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
thetavec_higher_ncpc <- thetavec_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
l_theta_ci_ncpc <- approx(thetaprof_lower_ncpc, thetavec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_theta_ci_ncpc <- approx(thetaprof_higher_ncpc, thetavec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_ncpc_log_sn <- c(opt_log_bbinom_ncpc$par[1], l_a_ci_ncpc$y, r_a_ci_ncpc$y) 
x50_ci_ncpc_log_sn <- c(opt_log_bbinom_ncpc$par[2], l_x50_ci_ncpc$y, r_x50_ci_ncpc$y)
c_ci_ncpc_log_sn <- c(opt_log_bbinom_ncpc$par[3], l_c_ci_ncpc$y, r_c_ci_ncpc$y)
theta_ci_ncpc_log_sn <- c(opt_log_bbinom_ncpc$par[4], l_theta_ci_ncpc$y, r_theta_ci_ncpc$y)

a_ci_ncpc_log_sn
x50_ci_ncpc_log_sn
c_ci_ncpc_log_sn
theta_ci_ncpc_log_sn


##              ##
## CNE function ##
##              ##

#
# Model setup
#

# Setup CNE function
mlsp_fun_ncpc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse((1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3])))) * exp(-par[1] * x) < 0.999,
               (1/exp(-par[1]*(par[2] + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * par[3])))) * exp(-par[1] * x), 0.999) 
  
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[4], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_ncpc <- c(b = 0.1, x100 = -10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_ncpc <- optim(par = mlsp_par_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, 
                     db1 = dat_12$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_12$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_12$season >= 4, 
                     db4 = dat_12$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_ncpc <- mlsp_opt_ncpc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_12[dat_12$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((0) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((0) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((1) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((1) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((2) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((2) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((3) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((3) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse((1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((4) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x) < 0.999,
                                        (1/exp(-mlsp_coefs_ncpc[1]*(mlsp_coefs_ncpc[2] + ((4) * mlsp_coefs_ncpc[3])))) * exp(-mlsp_coefs_ncpc[1] * x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C32. CNE function through the negative clairvoyance & positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_b_fun_ncpc <- function(par, x, z, n, b, db1, db2, db3, db4){
  b = b # This parameter is fixed
  x100 = par[1] # These parameters will be optimized
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999) 
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_x100_fun_ncpc <- function(par, x, z, n, x100, db1, db2, db3, db4){
  b = par[1]
  x100 = x100
  c = par[2]
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_ncpc <- function(par, x, z, n, c, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = c
  theta = par[3]
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_ncpc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  b = par[1]
  x100 = par[2]
  c = par[3]
  theta = theta
  mu <- ifelse((1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x) < 0.999,
               (1/exp(-b*(x100 + ((1*db1 + 1*db2 + 1*db3 + 1*db4 ) * c)))) * exp(-b * x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_b_ncpc = c(mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4]) # Starting parameters are the estimated parameters of the optimized function
bvec_ncpc = seq(0.0001, 0.30, length = 100) # Set vector for the fixed parameter
bprof_ncpc = numeric(100) # Prepare profile of likelihood values

pars_x100_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4])
x100vec_ncpc = seq(-120, 100, length = 100)
x100prof_ncpc = numeric(100)

pars_c_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[4])
cvec_ncpc = seq(0.1, 50, length = 100)
cprof_ncpc = numeric(100)

pars_theta_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3])
thetavec_ncpc = seq(0.01, 10, length = 100)
thetaprof_ncpc = numeric(100)

# Optimize original function without fixed parameters
pars_2_ncpc <- c(b = 0.1, x100 = -10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_ncpc <- optim(par = pars_2_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n,
                           db1 = dat_12$season >= 2, 
                           db2 = dat_12$season >= 3, 
                           db3 = dat_12$season >= 4, 
                           db4 = dat_12$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
b_coefs_vec_ncpc <- list(rep(list(), times = 100))
x100_coefs_vec_ncpc <- list(rep(list(), times = 100))
c_coefs_vec_ncpc <- list(rep(list(), times = 100))
theta_coefs_vec_ncpc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_b_ncpc <- optim(par = pars_b_ncpc, fn = ci_b_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, b = bvec_ncpc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_12$season >= 2, 
                    db2 = dat_12$season >= 3, 
                    db3 = dat_12$season >= 4, 
                    db4 = dat_12$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  b_coefs_vec_ncpc[[j]] = opt_b_ncpc$par # Store parameters for control
  # print(opt_a$convergence)
  bprof_ncpc[j] <- opt_b_ncpc$value # Set likelihood value for profile
  
  opt_x100_ncpc <- optim(par = pars_x100_ncpc, fn = ci_x100_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, x100 = x100vec_ncpc[j],
                       db1 = dat_12$season >= 2, 
                       db2 = dat_12$season >= 3, 
                       db3 = dat_12$season >= 4, 
                       db4 = dat_12$season >= 5, 
                       
                       control = list(maxit = 10000), method = "Nelder-Mead")
  x100prof_ncpc[j] <- opt_x100_ncpc$value
  x100_coefs_vec_ncpc[[j]] = opt_x100_ncpc$par
  
  opt_c_ncpc <- optim(par = pars_c_ncpc, fn = ci_c_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, c = cvec_ncpc[j],
                    db1 = dat_12$season >= 2, 
                    db2 = dat_12$season >= 3, 
                    db3 = dat_12$season >= 4, 
                    db4 = dat_12$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_ncpc[j] <- opt_c_ncpc$value
  c_coefs_vec_ncpc[[j]] = opt_c_ncpc$par
  
  opt_theta_ncpc <- optim(par = pars_theta_ncpc, fn = ci_theta_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, theta = thetavec_ncpc[j], 
                        db1 = dat_12$season >= 2, 
                        db2 = dat_12$season >= 3, 
                        db3 = dat_12$season >= 4, 
                        db4 = dat_12$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_ncpc[j] <- opt_theta_ncpc$value
  theta_coefs_vec_ncpc[[j]] = opt_theta_ncpc$par
}

# Set the 95% confidence limits
bprof_lower_ncpc <- bprof_ncpc[1:which.min(bprof_ncpc)] # Likelihood values below the best estimate
bvec_lower_ncpc <- bvec_ncpc[1:which.min(bprof_ncpc)] # Parameter values below the best estimate
bprof_higher_ncpc <- bprof_ncpc[which.min(bprof_ncpc):length(bprof_ncpc)] # Likelihood values above the best estimate
bvec_higher_ncpc <- bvec_ncpc[which.min(bprof_ncpc):length(bprof_ncpc)] # Parameter values above the best estimate
l_b_ci_ncpc <- approx(bprof_lower_ncpc, bvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_b_ci_ncpc <- approx(bprof_higher_ncpc, bvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

x100prof_lower_ncpc <- x100prof_ncpc[1:which.min(x100prof_ncpc)]
x100vec_lower_ncpc <- x100vec_ncpc[1:which.min(x100prof_ncpc)]
x100prof_higher_ncpc <- x100prof_ncpc[which.min(x100prof_ncpc):length(x100prof_ncpc)]
x100vec_higher_ncpc <- x100vec_ncpc[which.min(x100prof_ncpc):length(x100prof_ncpc)]
l_x100_ci_ncpc <- approx(x100prof_lower_ncpc, x100vec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_x100_ci_ncpc <- approx(x100prof_higher_ncpc, x100vec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

cprof_lower_ncpc <- cprof_ncpc[1:which.min(cprof_ncpc)]
cvec_lower_ncpc <- cvec_ncpc[1:which.min(cprof_ncpc)]
cprof_higher_ncpc <- cprof_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
cvec_higher_ncpc <- cvec_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
l_c_ci_ncpc <- approx(cprof_lower_ncpc, cvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_c_ci_ncpc <- approx(cprof_higher_ncpc, cvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

thetaprof_lower_ncpc <- thetaprof_ncpc[1:which.min(thetaprof_ncpc)]
thetavec_lower_ncpc <- thetavec_ncpc[1:which.min(thetaprof_ncpc)]
thetaprof_higher_ncpc <- thetaprof_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
thetavec_higher_ncpc <- thetavec_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
l_theta_ci_ncpc <- approx(thetaprof_lower_ncpc, thetavec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)
r_theta_ci_ncpc <- approx(thetaprof_higher_ncpc, thetavec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 3)/2)

# Set the parameter estimate with low and high 95% CI's
b_ci_ncpc_cne_sn <- c(opt_log_bbinom_ncpc$par[1], l_b_ci_ncpc$y, r_b_ci_ncpc$y) 
x100_ci_ncpc_cne_sn <- c(opt_log_bbinom_ncpc$par[2], l_x100_ci_ncpc$y, r_x100_ci_ncpc$y)
c_ci_ncpc_cne_sn <- c(opt_log_bbinom_ncpc$par[3], l_c_ci_ncpc$y, r_c_ci_ncpc$y)
theta_ci_ncpc_cne_sn <- c(opt_log_bbinom_ncpc$par[4], l_theta_ci_ncpc$y, r_theta_ci_ncpc$y)

b_ci_ncpc_cne_sn
x100_ci_ncpc_cne_sn
c_ci_ncpc_cne_sn
theta_ci_ncpc_cne_sn


##               ##
## Hill function ##
##               ##

#
# Model setup
#

# Setup Hill function
mlsp_fun_ncpc <- function(par, x, z, n, db1, db2, db3, db4){ 
  mu <- ifelse(x < (par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4)), par[1], 
               (par[1]*(1 - ((x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2])/
                          (par[3]^par[2] + (x - par[4]*(1*db1 + 1*db2 + 1*db3 + 1*db4))^par[2]))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par[5], log = TRUE)) # nll calculation
  return(nll)
}

# Setup parameter starting values
mlsp_par_ncpc <- c(a = 0.99, n = 2, h = 10, c = 10, theta = 1)


#
# Model fitting
#

# Model fit
mlsp_opt_ncpc <- optim(par = mlsp_par_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, 
                     db1 = dat_12$season >= 2, # db1 is TRUE for 2014 and above, otherwise FALSE
                     db2 = dat_12$season >= 3, # db2 is TRUE for 2015 and above, otherwise FALSE
                     db3 = dat_12$season >= 4, 
                     db4 = dat_12$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")

mlsp_coefs_ncpc <- mlsp_opt_ncpc$par # Set parameters

# Plot the lines
ggplot() +
  geom_point(data = dat_12[dat_12$season == 1,], aes(x = dist, y = prop, color = "2013/2014"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 2,], aes(x = dist, y = prop, color = "2014/2015"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 3,], aes(x = dist, y = prop, color = "2015/2016"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 4,], aes(x = dist, y = prop, color = "2016/2017"), shape = 1) +
  geom_point(data = dat_12[dat_12$season == 5,], aes(x = dist, y = prop, color = "2017/2018"), shape = 1) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(0)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(0))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(0))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(1)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(1))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(1))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(2)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(2))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(2))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(3)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(3))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(3))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)ifelse(x < (mlsp_coefs_ncpc[4]*(4)), mlsp_coefs_ncpc[1], 
                                        (mlsp_coefs_ncpc[1]*(1 - ((x - mlsp_coefs_ncpc[4]*(4))^mlsp_coefs_ncpc[2])/
                                                   (mlsp_coefs_ncpc[3]^mlsp_coefs_ncpc[2] + (x - mlsp_coefs_ncpc[4]*(4))^mlsp_coefs_ncpc[2])))), 
                aes(color = "2017/2018", size = "s")) +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure C33. Hill function through the negative clairvoyance & positive carry-over data with seasons.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


#
# 95% CI calculation
#

# Create NLL functions for CI calculations
ci_a_fun_ncpc <- function(par, x, z, n, a, db1, db2, db3, db4){
  a = a # This parameter is fixed
  pn = par[1] # pn for parameter n, since there already is an n for the number of samples
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_pn_fun_ncpc <- function(par, x, z, n, pn, db1, db2, db3, db4){
  a = par[1] # This parameter is fixed
  pn = pn
  h = par[2] # These parameters will be optimized
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_h_fun_ncpc <- function(par, x, z, n, h, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = h
  c = par[3]
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_c_fun_ncpc <- function(par, x, z, n, c, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = c
  theta = par[4]
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

ci_theta_fun_ncpc <- function(par, x, z, n, theta, db1, db2, db3, db4){
  a = par[1]
  pn = par[2]
  h = par[3]
  c = par[4]
  theta = theta
  mu <- ifelse(x < (c*(1*db1 + 1*db2 + 1*db3 + 1*db4)), a, 
               (a*(1 - ((x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn)/
                     (h^pn + (x - c*(1*db1 + 1*db2 + 1*db3 + 1*db4))^pn))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = theta, log = TRUE))
  return(nll)
}

# Set the starting parameters for the optimization of the NLL functions for CI calculations
pars_a_ncpc = c(mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4], mlsp_coefs_ncpc[5]) # Starting parameters are the estimated parameters of the optimized function
avec_ncpc = seq(0.01, 0.99, length = 100) # Set vector for the fixed parameter
aprof_ncpc = numeric(100) # Prepare profile of likelihood values

pars_pn_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4], mlsp_coefs_ncpc[5])
pnvec_ncpc = seq(0.1, 10, length = 100)
pnprof_ncpc = numeric(100)

pars_h_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[4], mlsp_coefs_ncpc[5])
hvec_ncpc = seq(2, 100, length = 100)
hprof_ncpc = numeric(100)

pars_c_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[5])
cvec_ncpc = seq(0.1, 50, length = 100)
cprof_ncpc = numeric(100)

pars_theta_ncpc = c(mlsp_coefs_ncpc[1], mlsp_coefs_ncpc[2], mlsp_coefs_ncpc[3], mlsp_coefs_ncpc[4])
thetavec_ncpc = seq(0.01, 10, length = 100)
thetaprof_ncpc = numeric(100)

# Optimize original function without fixed parameters
pars_2_ncpc <- c(a = 0.99, pn = 2, h = 10, c = 10, theta = 1) # Set starting parameters for optimization of function without fixed parameters
opt_log_bbinom_ncpc <- optim(par = pars_2_ncpc, fn = mlsp_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n,
                           db1 = dat_12$season >= 2, 
                           db2 = dat_12$season >= 3, 
                           db3 = dat_12$season >= 4, 
                           db4 = dat_12$season >= 5, 
                           
                           control = list(maxit = 10000), method = "Nelder-Mead")

# Prepare lists to store parameter values for each fixed value, for control
a_coefs_vec_ncpc <- list(rep(list(), times = 100))
pn_coefs_vec_ncpc <- list(rep(list(), times = 100))
h_coefs_vec_ncpc <- list(rep(list(), times = 100))
c_coefs_vec_ncpc <- list(rep(list(), times = 100))
theta_coefs_vec_ncpc <- list(rep(list(), times = 100))

# Run optimization for each fixed value
for(j in 1:100){ # 100 fixed values for each parameter
  opt_a_ncpc <- optim(par = pars_a_ncpc, fn = ci_a_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, a = avec_ncpc[j], # Optimize all parameters except for the fixed parameter. Run for every fixed value.
                    db1 = dat_12$season >= 2, 
                    db2 = dat_12$season >= 3, 
                    db3 = dat_12$season >= 4, 
                    db4 = dat_12$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  a_coefs_vec_ncpc[[j]] = opt_a_ncpc$par # Store parameters for control
  aprof_ncpc[j] <- opt_a_ncpc$value # Set likelihood value for profile
  
  opt_pn_ncpc <- optim(par = pars_pn_ncpc, fn = ci_pn_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, pn = pnvec_ncpc[j],
                     db1 = dat_12$season >= 2, 
                     db2 = dat_12$season >= 3, 
                     db3 = dat_12$season >= 4, 
                     db4 = dat_12$season >= 5, 
                     
                     control = list(maxit = 10000), method = "Nelder-Mead")
  pnprof_ncpc[j] <- opt_pn_ncpc$value
  pn_coefs_vec_ncpc[[j]] = opt_pn_ncpc$par
  
  opt_h_ncpc <- optim(par = pars_h_ncpc, fn = ci_h_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, h = hvec_ncpc[j],
                    db1 = dat_12$season >= 2, 
                    db2 = dat_12$season >= 3, 
                    db3 = dat_12$season >= 4, 
                    db4 = dat_12$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  hprof_ncpc[j] <- opt_h_ncpc$value
  h_coefs_vec_ncpc[[j]] = opt_h_ncpc$par
  
  opt_c_ncpc <- optim(par = pars_c_ncpc, fn = ci_c_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, c = cvec_ncpc[j],
                    db1 = dat_12$season >= 2, 
                    db2 = dat_12$season >= 3, 
                    db3 = dat_12$season >= 4, 
                    db4 = dat_12$season >= 5, 
                    
                    control = list(maxit = 10000), method = "Nelder-Mead")
  cprof_ncpc[j] <- opt_c_ncpc$value
  c_coefs_vec_ncpc[[j]] = opt_c_ncpc$par
  
  opt_theta_ncpc <- optim(par = pars_theta_ncpc, fn = ci_theta_fun_ncpc, x = dat_12$dist, z = dat_12$pos, n = dat_12$n, theta = thetavec_ncpc[j], 
                        db1 = dat_12$season >= 2, 
                        db2 = dat_12$season >= 3, 
                        db3 = dat_12$season >= 4, 
                        db4 = dat_12$season >= 5, 
                        
                        control = list(maxit = 10000), method = "Nelder-Mead")
  thetaprof_ncpc[j] <- opt_theta_ncpc$value
  theta_coefs_vec_ncpc[[j]] = opt_theta_ncpc$par
}

# Set the 95% confidence limits
aprof_lower_ncpc <- aprof_ncpc[1:which.min(aprof_ncpc)] # Likelihood values below the best estimate
avec_lower_ncpc <- avec_ncpc[1:which.min(aprof_ncpc)] # Parameter values below the best estimate
aprof_higher_ncpc <- aprof_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Likelihood values above the best estimate
avec_higher_ncpc <- avec_ncpc[which.min(aprof_ncpc):length(aprof_ncpc)] # Parameter values above the best estimate
l_a_ci_ncpc <- approx(aprof_lower_ncpc, avec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2) # Calculate the 95% lower confidence limit with a chisquare distribution with 3 df
r_a_ci_ncpc <- approx(aprof_higher_ncpc, avec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2) # Calculate the 95% upper confidence limit with a chisquare dsitribution with 3 df

pnprof_lower_ncpc <- pnprof_ncpc[1:which.min(pnprof_ncpc)]
pnvec_lower_ncpc <- pnvec_ncpc[1:which.min(pnprof_ncpc)]
pnprof_higher_ncpc <- pnprof_ncpc[which.min(pnprof_ncpc):length(pnprof_ncpc)]
pnvec_higher_ncpc <- pnvec_ncpc[which.min(pnprof_ncpc):length(pnprof_ncpc)]
l_pn_ci_ncpc <- approx(pnprof_lower_ncpc, pnvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_pn_ci_ncpc <- approx(pnprof_higher_ncpc, pnvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

hprof_lower_ncpc <- hprof_ncpc[1:which.min(hprof_ncpc)]
hvec_lower_ncpc <- hvec_ncpc[1:which.min(hprof_ncpc)]
hprof_higher_ncpc <- hprof_ncpc[which.min(hprof_ncpc):length(hprof_ncpc)]
hvec_higher_ncpc <- hvec_ncpc[which.min(hprof_ncpc):length(hprof_ncpc)]
l_h_ci_ncpc <- approx(hprof_lower_ncpc, hvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_h_ci_ncpc <- approx(hprof_higher_ncpc, hvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

cprof_lower_ncpc <- cprof_ncpc[1:which.min(cprof_ncpc)]
cvec_lower_ncpc <- cvec_ncpc[1:which.min(cprof_ncpc)]
cprof_higher_ncpc <- cprof_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
cvec_higher_ncpc <- cvec_ncpc[which.min(cprof_ncpc):length(cprof_ncpc)]
l_c_ci_ncpc <- approx(cprof_lower_ncpc, cvec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_c_ci_ncpc <- approx(cprof_higher_ncpc, cvec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

thetaprof_lower_ncpc <- thetaprof_ncpc[1:which.min(thetaprof_ncpc)]
thetavec_lower_ncpc <- thetavec_ncpc[1:which.min(thetaprof_ncpc)]
thetaprof_higher_ncpc <- thetaprof_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
thetavec_higher_ncpc <- thetavec_ncpc[which.min(thetaprof_ncpc):length(thetaprof_ncpc)]
l_theta_ci_ncpc <- approx(thetaprof_lower_ncpc, thetavec_lower_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)
r_theta_ci_ncpc <- approx(thetaprof_higher_ncpc, thetavec_higher_ncpc, xout = opt_log_bbinom_ncpc$value + qchisq(0.95, 4)/2)

# Set the parameter estimate with low and high 95% CI's
a_ci_ncpc_hill_sn <- c(opt_log_bbinom_ncpc$par[1], l_a_ci_ncpc$y, r_a_ci_ncpc$y) 
pn_ci_ncpc_hill_sn <- c(opt_log_bbinom_ncpc$par[2], l_pn_ci_ncpc$y, r_pn_ci_ncpc$y)
h_ci_ncpc_hill_sn <- c(opt_log_bbinom_ncpc$par[3], l_h_ci_ncpc$y, r_h_ci_ncpc$y)
c_ci_ncpc_hill_sn <- c(opt_log_bbinom_ncpc$par[4], l_c_ci_ncpc$y, r_c_ci_ncpc$y)
theta_ci_ncpc_hill_sn <- c(opt_log_bbinom_ncpc$par[5], l_theta_ci_ncpc$y, r_theta_ci_ncpc$y)

a_ci_ncpc_hill_sn
pn_ci_ncpc_hill_sn
h_ci_ncpc_hill_sn
c_ci_ncpc_hill_sn
theta_ci_ncpc_hill_sn


##             ##
## All results ##
##             ##

a_ci_od_log_yr
x50_ci_od_log_yr
c_ci_od_log_yr
theta_ci_od_log_yr

b_ci_od_cne_yr
x100_ci_od_cne_yr
c_ci_od_cne_yr
theta_ci_od_cne_yr

a_ci_od_hill_yr
pn_ci_od_hill_yr
h_ci_od_hill_yr
c_ci_od_hill_yr
theta_ci_od_hill_yr


a_ci_pc_log_yr
x50_ci_pc_log_yr
c_ci_pc_log_yr
theta_ci_pc_log_yr

b_ci_pc_cne_yr
x100_ci_pc_cne_yr
c_ci_pc_cne_yr
theta_ci_pc_cne_yr

a_ci_pc_hill_yr
pn_ci_pc_hill_yr
h_ci_pc_hill_yr
c_ci_pc_hill_yr
theta_ci_pc_hill_yr


a_ci_nc_log_yr
x50_ci_nc_log_yr
c_ci_nc_log_yr
theta_ci_nc_log_yr

b_ci_nc_cne_yr
x100_ci_nc_cne_yr
c_ci_nc_cne_yr
theta_ci_nc_cne_yr

a_ci_nc_hill_yr
pn_ci_nc_hill_yr
h_ci_nc_hill_yr
c_ci_nc_hill_yr
theta_ci_nc_hill_yr


a_ci_ncpc_log_yr
x50_ci_ncpc_log_yr
c_ci_ncpc_log_yr
theta_ci_ncpc_log_yr

b_ci_ncpc_cne_yr
x100_ci_ncpc_cne_yr
c_ci_ncpc_cne_yr
theta_ci_ncpc_cne_yr

a_ci_ncpc_hill_yr
pn_ci_ncpc_hill_yr
h_ci_ncpc_hill_yr
c_ci_ncpc_hill_yr
theta_ci_ncpc_hill_yr


a_ci_od_log_sn
x50_ci_od_log_sn
c_ci_od_log_sn
theta_ci_od_log_sn

b_ci_od_cne_sn
x100_ci_od_cne_sn
c_ci_od_cne_sn
theta_ci_od_cne_sn

a_ci_od_hill_sn
pn_ci_od_hill_sn
h_ci_od_hill_sn
c_ci_od_hill_sn
theta_ci_od_hill_sn


a_ci_pc_log_sn
x50_ci_pc_log_sn
c_ci_pc_log_sn
theta_ci_pc_log_sn

b_ci_pc_cne_sn
x100_ci_pc_cne_sn
c_ci_pc_cne_sn
theta_ci_pc_cne_sn

a_ci_pc_hill_sn
pn_ci_pc_hill_sn
h_ci_pc_hill_sn
c_ci_pc_hill_sn
theta_ci_pc_hill_sn


a_ci_nc_log_sn
x50_ci_nc_log_sn
c_ci_nc_log_sn
theta_ci_nc_log_sn

b_ci_nc_cne_sn
x100_ci_nc_cne_sn
c_ci_nc_cne_sn
theta_ci_nc_cne_sn

a_ci_nc_hill_sn
pn_ci_nc_hill_sn
h_ci_nc_hill_sn
c_ci_nc_hill_sn
theta_ci_nc_hill_sn


a_ci_ncpc_log_sn
x50_ci_ncpc_log_sn
c_ci_ncpc_log_sn
theta_ci_ncpc_log_sn

b_ci_ncpc_cne_sn
x100_ci_ncpc_cne_sn
c_ci_ncpc_cne_sn
theta_ci_ncpc_cne_sn

a_ci_ncpc_hill_sn
pn_ci_ncpc_hill_sn
h_ci_ncpc_hill_sn
c_ci_ncpc_hill_sn
theta_ci_ncpc_hill_sn