###########
##       ##
## Setup ##
##       ##
###########

rm(list = ls())
setwd("C:/Users/dbkot/OneDrive/Documents/Wageningen University/Msc Biology WUR/Thesis/R")

library(readxl)
library(reshape2)
library(geosphere)
library(ggplot2)
library(emdbook)
library(tidyverse)
library(dplyr)

# Set ggplot graph theme
theme_set(theme_bw())


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
    x = "Distance (km)", 
    y = "Proportion positives") +
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
  labs(x = "Distance (km)", 
       y = "Proportion positives") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.05) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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
  labs(x = "Distance (km)", 
       y = "Proportion positives") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.01) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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