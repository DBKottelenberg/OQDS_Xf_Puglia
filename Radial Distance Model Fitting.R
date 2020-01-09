###########
##       ##
## Setup ##
##       ##
###########

rm(list = ls())
setwd("C:/Users/dbkot/OneDrive - WageningenUR/Msc Biology WUR/Thesis/R")

library(readxl)
library(geosphere)
library(ggplot2)
library(emdbook)
library(tidyverse)
library(dplyr)
library(segmented)

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

# Set date factor correctly, drop NA and incorrect lon/lat values
dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter((lon != 0) | (lat != 0))

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
dc <- 0:(max_dist - 1) # Number of dc's is the maximum distance a measurement is done.

# Prepare list for number of positives and total number of measurements in each dc.
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(years)) 
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(years)) 

dat_2 <- list() # Prepare data list

for(i in 1:length(years)){ # For-loop running for every year
  for(j in dc){ # For-loop running for every distance circle
    pos[[i]][j + 1] <- sum((dat[[i]]$result[dat[[i]]$dist >= j & (dat[[i]]$dist < j+1)])) # For every year, for every dc, calculate the total number of positives
    n[[i]][j + 1] <- length((dat[[i]]$result[dat[[i]]$dist >= j & (dat[[i]]$dist < j+1)])) # For every year, for every dc, calculate the total number of measurements
  }
  dat_2[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% # Create list of data frames for every year, which have a row for every distance circle, set proportion of positives.
    mutate(prop = pos/n) 
}

# Plot of the data
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, y = prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure A1. A graph with the datapoints of the original 'year' data aggregated with 
       distance circles.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(plot.caption = element_text(hjust = 0))


#
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_nexp_binom <- list(c_2013 <- c(NA, NA),
                         c_2014 <- c(NA, NA),
                         c_2015 <- c(NA, NA),
                         c_2016 <- c(NA, NA),
                         c_2017 <- c(NA, NA),
                         c_2018 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_2[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(a = 0.1, b = 0.1) # Starting parameters
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n) # Optimization
  c_opt <- opt_nexp_binom$par # Store optimized parameters
  coefs_nexp_binom[[i]] <- c_opt 
  comb_AIC_nexp_binom[i] <- 2*2 + 2*as.numeric(opt_nexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[6]]['a']*exp(-coefs_nexp_binom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution;
       Type of data: Original data", 
    caption = "Figure A2. A graph with the datapoints of the original data and lines as fitted with a 
       negative exponential function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                          c_2014 <- c(NA, NA, NA),
                          c_2015 <- c(NA, NA, NA),
                          c_2016 <- c(NA, NA, NA),
                          c_2017 <- c(NA, NA, NA),
                          c_2018 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(years))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(years)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*3 + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[6]]['a']*exp(-coefs_nexp_bbinom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution; 
       Type of data: Original data", 
    caption = "Figure A3. A graph with the datapoints of the original data and lines as fitted with a negative 
       exponential function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_2013 <- c(NA, NA),
                        c_2014 <- c(NA, NA),
                        c_2015 <- c(NA, NA),
                        c_2016 <- c(NA, NA),
                        c_2017 <- c(NA, NA),
                        c_2018 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(years))

fun_log_binom <- function(par, x, z, n){
  mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 0.1, b = 0.1)

for(i in 1:length(years)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*2 + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_binom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(exp(coefs_log_binom[[1]]['a'] + coefs_log_binom[[1]]['b']*x))/(1 + exp(coefs_log_binom[[1]]['a'] + coefs_log_binom[[1]]['b']*x)), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[2]]['a'] + coefs_log_binom[[2]]['b']*x))/(1 + exp(coefs_log_binom[[2]]['a'] + coefs_log_binom[[2]]['b']*x)), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[3]]['a'] + coefs_log_binom[[3]]['b']*x))/(1 + exp(coefs_log_binom[[3]]['a'] + coefs_log_binom[[3]]['b']*x)), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[4]]['a'] + coefs_log_binom[[4]]['b']*x))/(1 + exp(coefs_log_binom[[4]]['a'] + coefs_log_binom[[4]]['b']*x)), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[5]]['a'] + coefs_log_binom[[5]]['b']*x))/(1 + exp(coefs_log_binom[[5]]['a'] + coefs_log_binom[[5]]['b']*x)), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[6]]['a'] + coefs_log_binom[[6]]['b']*x))/(1 + exp(coefs_log_binom[[6]]['a'] + coefs_log_binom[[6]]['b']*x)), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution; 
       Type of data: Original data",
    caption = "Figure A4. A graph with the datapoints of the original data and lines as fitted with a logistic 
       function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(years))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 0.1, b = 0.1, theta = 1)

for(i in 1:length(years)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n)
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[1]]['a'] + coefs_log_bbinom[[1]]['b']*x))/(1 + exp(coefs_log_bbinom[[1]]['a'] + coefs_log_bbinom[[1]]['b']*x)), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[2]]['a'] + coefs_log_bbinom[[2]]['b']*x))/(1 + exp(coefs_log_bbinom[[2]]['a'] + coefs_log_bbinom[[2]]['b']*x)), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[3]]['a'] + coefs_log_bbinom[[3]]['b']*x))/(1 + exp(coefs_log_bbinom[[3]]['a'] + coefs_log_bbinom[[3]]['b']*x)), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[4]]['a'] + coefs_log_bbinom[[4]]['b']*x))/(1 + exp(coefs_log_bbinom[[4]]['a'] + coefs_log_bbinom[[4]]['b']*x)), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[5]]['a'] + coefs_log_bbinom[[5]]['b']*x))/(1 + exp(coefs_log_bbinom[[5]]['a'] + coefs_log_bbinom[[5]]['b']*x)), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[6]]['a'] + coefs_log_bbinom[[6]]['b']*x))/(1 + exp(coefs_log_bbinom[[6]]['a'] + coefs_log_bbinom[[6]]['b']*x)), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution; 
       Type of data: Original data", 
    caption = "Figure A5. A graph with the datapoints of the original data and lines as fitted with a logistic 
       function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
##  CNE, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_binom <- list(c_2013 <- c(NA, NA),
                           c_2014 <- c(NA, NA),
                           c_2015 <- c(NA, NA),
                           c_2016 <- c(NA, NA),
                           c_2017 <- c(NA, NA),
                           c_2018 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_2[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.1, x100 = -10) # Starting parameters
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_binom$par # Store optimized parameters
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution; 
       Type of data: Original data", 
    caption = "Figure A6. A graph with the datapoints of the original data and lines as fitted with a CNE function 
       with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                            c_2014 <- c(NA, NA, NA),
                            c_2015 <- c(NA, NA, NA),
                            c_2016 <- c(NA, NA, NA),
                            c_2017 <- c(NA, NA, NA),
                            c_2018 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_2[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.1, x100 = -10, theta = 1) # Starting parameters
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_bbinom$par # Store optimized parameters
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

# Create the plot
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution; 
       Type of data: Original data", 
    caption = "Figure A7. A graph with the datapoints of the original data and lines as fitted with a CNE function 
       with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function/binomial distribution
##
coefs_hill_binom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(years))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(years)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 5000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[6]]['a']*(1 - (x^(coefs_hill_binom[[6]]['n'])/(coefs_hill_binom[[6]]['h']^(coefs_hill_binom[[6]]['n']) + x^(coefs_hill_binom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Binomial distribution; 
       Type of data: Original data", 
    caption = "Figure A8. A graph with the datapoints of the original data and lines as fitted with a Hill function 
       with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function/beta-binomial distribution
##
coefs_hill_bbinom <- list(c_2013 <- c(NA, NA, NA, NA),
                          c_2014 <- c(NA, NA, NA, NA),
                          c_2015 <- c(NA, NA, NA, NA),
                          c_2016 <- c(NA, NA, NA, NA),
                          c_2017 <- c(NA, NA, NA, NA),
                          c_2018 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(years))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(years)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[6]]['a']*(1 - (x^(coefs_hill_bbinom[[6]]['n'])/(coefs_hill_bbinom[[6]]['h']^(coefs_hill_bbinom[[6]]['n']) + x^(coefs_hill_bbinom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution; 
       Type of data: Original data", 
    caption = "Figure A9. A graph with the datapoints of the original data and lines as fitted with a Hill function 
       with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


############################
# Positive carry-over data #
############################

#
# Data morphing and distance circle creation
#

dat_3 <- dat

# make positives carry over
for(i in 2:length(years)){ # For-loop through every year from 2014
  co <- which(dat_3[[i-1]]$result == 1) # Which lines of the previous year have a result of 1
  
  dat_3[[i]] <- rbind(dat_3[[i]], dat_3[[i-1]][co,]) # Include the lines of the previous year with a result of 1 with the current year
  co_2 <- which(dat_3[[i]]$lon == dat_3[[i-1]]$lon[co] & 
                  dat_3[[i]]$lat == dat_3[[i-1]]$lat[co] & 
                  dat_3[[i]]$year == (2012 + i)) # Which lines of this year have the same coordinates of the positives last year
  
  if(length(co_2) >= 1){ # If this is not 0, then there has been a measurement at these coordinates this year as well as last year, or multiple times 
                         # in a single year (with the Apulia data, only multiple measurements in 2013)
    dat_3[[i]] <- dat_3[[i]][-co_2,] # Remove all the duplicates. If the duplicate measurement from this year is a 0, it will still be a 1,
                                     # because I assume that once the disease is there, it will not leave. I assume that 0 to be a false negative.
  }
}

# Set distance circles
dc <- 0:(max_dist - 1)
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(years))
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(years))

dat_4 <- list()

for(i in 1:length(years)){
  for(j in dc){
    pos[[i]][j + 1] <- sum((dat_3[[i]]$result[dat_3[[i]]$dist >= j & (dat_3[[i]]$dist < j+1)]))
    n[[i]][j + 1] <- length((dat_3[[i]]$result[dat_3[[i]]$dist >= j & (dat_3[[i]]$dist < j+1)]))
  }
  dat_4[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n)
}

# Plot of the data
ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, y = prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  labs(subtitle = "Positive carry-over data",
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A10. A graph with the datapoints of the positive carry-over data aggregated with 
       distance circles.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(plot.caption = element_text(hjust = 0))


#
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames.
coefs_nexp_binom <- list(c_2013 <- c(NA, NA),
                         c_2014 <- c(NA, NA),
                         c_2015 <- c(NA, NA),
                         c_2016 <- c(NA, NA),
                         c_2017 <- c(NA, NA),
                         c_2018 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1)
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n)
  c_opt <- opt_nexp_binom$par
  coefs_nexp_binom[[i]] <- c_opt
  comb_AIC_nexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_nexp_binom$value)
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[6]]['a']*exp(-coefs_nexp_binom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A11. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       negative exponential function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                          c_2014 <- c(NA, NA, NA),
                          c_2015 <- c(NA, NA, NA),
                          c_2016 <- c(NA, NA, NA),
                          c_2017 <- c(NA, NA, NA),
                          c_2018 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(years))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(years)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[6]]['a']*exp(-coefs_nexp_bbinom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A12. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       negative exponential function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_2013 <- c(NA, NA),
                        c_2014 <- c(NA, NA),
                        c_2015 <- c(NA, NA),
                        c_2016 <- c(NA, NA),
                        c_2017 <- c(NA, NA),
                        c_2018 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(years))

fun_log_binom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['b'])))))
  # mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, b = -30)

for(i in 1:length(years)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*length(pars) + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_binom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[1]]['a'] * (x - coefs_log_binom[[1]]['b'])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[2]]['a'] * (x - coefs_log_binom[[2]]['b'])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[3]]['a'] * (x - coefs_log_binom[[3]]['b'])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[4]]['a'] * (x - coefs_log_binom[[4]]['b'])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[5]]['a'] * (x - coefs_log_binom[[5]]['b'])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[6]]['a'] * (x - coefs_log_binom[[6]]['b'])))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A13. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       logistic function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(years))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['c'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, c = 30, theta = 5)

for(i in 1:length(years)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[1]]['a'] * (x - coefs_log_bbinom[[1]]['c'])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[2]]['a'] * (x - coefs_log_bbinom[[2]]['c'])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[3]]['a'] * (x - coefs_log_bbinom[[3]]['c'])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[4]]['a'] * (x - coefs_log_bbinom[[4]]['c'])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[5]]['a'] * (x - coefs_log_bbinom[[5]]['c'])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[6]]['a'] * (x - coefs_log_bbinom[[6]]['c'])))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A14. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       logistic function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_binom <- list(c_2013 <- c(NA, NA),
                           c_2014 <- c(NA, NA),
                           c_2015 <- c(NA, NA),
                           c_2016 <- c(NA, NA),
                           c_2017 <- c(NA, NA),
                           c_2018 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_4[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.1, x100 = -10) # Starting parameters
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_binom$par # Store optimized parameters
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A15. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       CNE function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                            c_2014 <- c(NA, NA, NA),
                            c_2015 <- c(NA, NA, NA),
                            c_2016 <- c(NA, NA, NA),
                            c_2017 <- c(NA, NA, NA),
                            c_2018 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_4[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.05, x100 = -10, theta = 1) # Starting parameters
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_bbinom$par # Store optimized parameters
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

# Create the plot
ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A16. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       CNE function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function/binomial distribution
##
coefs_hill_binom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(years))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(years)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[6]]['a']*(1 - (x^(coefs_hill_binom[[6]]['n'])/(coefs_hill_binom[[6]]['h']^(coefs_hill_binom[[6]]['n']) + x^(coefs_hill_binom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A17. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       Hill function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))

##
## Hill function/beta-binomial distribution
##
coefs_hill_bbinom <- list(c_2013 <- c(NA, NA, NA, NA),
                          c_2014 <- c(NA, NA, NA, NA),
                          c_2015 <- c(NA, NA, NA, NA),
                          c_2016 <- c(NA, NA, NA, NA),
                          c_2017 <- c(NA, NA, NA, NA),
                          c_2018 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(years))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(years)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[6]]['a']*(1 - (x^(coefs_hill_bbinom[[6]]['n'])/(coefs_hill_bbinom[[6]]['h']^(coefs_hill_bbinom[[6]]['n']) + x^(coefs_hill_bbinom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution; 
       Type of data: positive carry-over data", 
    caption = "Figure A18. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       Hill function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##############################
# Negative clairvoyance data #
##############################

#
# Data morphing and distance circle creation
#

dat_5 <- dat

# Make data clairvoyant
for(i in length(years):2){ # Work backwards, from 2018 to 2014
  co <- which(dat_5[[i]]$result == 0) # Which of the data have a result of 0
  dat_5[[i-1]] <- rbind(dat_5[[i-1]], dat_5[[i]][co,]) # Combine the previous year with the negative results of the later year
  co_2 <- which(dat_5[[i-1]]$lon == dat_5[[i]]$lon[co] &     # Check for multiple samples at the same location
            dat_5[[i-1]]$lat == dat_5[[i]]$lat[co] &  
            dat_5[[i-1]]$year == (2012 + i)) 
  if(length(co_2) >= 1){ # If there are multiple samples at one location
    for(j in 1:length(co_2)){ # For every duplicate sample
      if(dat_5[[i-1]]$result[co_2[j]] == 1){ # If the duplicate from the year before has result = 1, the result at this location should be a 1
        co_3 <- which(dat_5[[i-1]]$lon == dat_5[[i]]$lon[co] & 
                        dat_5[[i-1]]$lat == dat_5[[i]]$lat[co] & 
                        dat_5[[i-1]]$year == (2012 + i))
        dat_5[[i-1]] <- dat_5[[i-1]][-co_3[j],]
      }
      if(dat_5[[i-1]]$result[co_2[j]] == 0){    # If the duplicate result is a 0, it will stay a 0.
        dat_5[[i-1]] <- dat_5[[i-1]][-co_2[j],] 
      }
    }
  } # Note that this should never be the case, since a 1 would have turned into a 0. This is indeed never the case with the Apulia data.
    # This part of code is a safeguard to be used in future data or simulated data.
}

# Set distance circles
dc <- 0:(max_dist - 1)
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(years))
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(years))

dat_6 <- list()

for(i in 1:length(years)){
  for(j in dc){
    pos[[i]][j + 1] <- sum((dat_5[[i]]$result[dat_5[[i]]$dist >= j & (dat_5[[i]]$dist < j+1)]))
    n[[i]][j + 1] <- length((dat_5[[i]]$result[dat_5[[i]]$dist >= j & (dat_5[[i]]$dist < j+1)]))
  }
  dat_6[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n)
}

# Plot of the data
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, y = prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  labs(subtitle = "Negative clairvoyance data",
    caption = "Figure A19. A graph with the datapoints of the negative clairvoyance data aggregated 
       with distance circles.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(plot.caption = element_text(hjust = 0))

# 
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames.
coefs_nexp_binom <- list(c_2013 <- c(NA, NA),
                         c_2014 <- c(NA, NA),
                         c_2015 <- c(NA, NA),
                         c_2016 <- c(NA, NA),
                         c_2017 <- c(NA, NA),
                         c_2018 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1)
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n)
  c_opt <- opt_nexp_binom$par
  coefs_nexp_binom[[i]] <- c_opt
  comb_AIC_nexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_nexp_binom$value)
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[6]]['a']*exp(-coefs_nexp_binom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(title = "Proportion of positives with increasing distance", 
       subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution; 
       Type of data: negative clairvoyance data", 
       caption = "Figure A20. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a negative exponential function with binomial distribution.",
       x = "Distance", 
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                          c_2014 <- c(NA, NA, NA),
                          c_2015 <- c(NA, NA, NA),
                          c_2016 <- c(NA, NA, NA),
                          c_2017 <- c(NA, NA, NA),
                          c_2018 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(years))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(years)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[6]]['a']*exp(-coefs_nexp_bbinom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution; 
       Type of data: negative clairvoyance data", 
    caption = "Figure A21. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a negative exponential function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_2013 <- c(NA, NA),
                        c_2014 <- c(NA, NA),
                        c_2015 <- c(NA, NA),
                        c_2016 <- c(NA, NA),
                        c_2017 <- c(NA, NA),
                        c_2018 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(years))

fun_log_binom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['b'])))))
  # mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, b = -30)

for(i in 1:length(years)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*length(pars) + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_binom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[1]]['a'] * (x - coefs_log_binom[[1]]['b'])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[2]]['a'] * (x - coefs_log_binom[[2]]['b'])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[3]]['a'] * (x - coefs_log_binom[[3]]['b'])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[4]]['a'] * (x - coefs_log_binom[[4]]['b'])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[5]]['a'] * (x - coefs_log_binom[[5]]['b'])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[6]]['a'] * (x - coefs_log_binom[[6]]['b'])))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution; 
       Type of data: negative clairvoyance data", 
    caption = "Figure A22. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a logistic function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(years))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['c'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, c = 30, theta = 5)

for(i in 1:length(years)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[1]]['a'] * (x - coefs_log_bbinom[[1]]['c'])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[2]]['a'] * (x - coefs_log_bbinom[[2]]['c'])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[3]]['a'] * (x - coefs_log_bbinom[[3]]['c'])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[4]]['a'] * (x - coefs_log_bbinom[[4]]['c'])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[5]]['a'] * (x - coefs_log_bbinom[[5]]['c'])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[6]]['a'] * (x - coefs_log_bbinom[[6]]['c'])))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution; 
       Type of data: negative clairvoyance data",
    caption = "Figure A23. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a logistic function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_binom <- list(c_2013 <- c(NA, NA),
                           c_2014 <- c(NA, NA),
                           c_2015 <- c(NA, NA),
                           c_2016 <- c(NA, NA),
                           c_2017 <- c(NA, NA),
                           c_2018 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_6[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.1, x100 = -10) # Starting parameters
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_binom$par # Store optimized parameters
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution 
       Type of data: negative clairvoyance data", 
    caption = "Figure A24. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a CNE function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                            c_2014 <- c(NA, NA, NA),
                            c_2015 <- c(NA, NA, NA),
                            c_2016 <- c(NA, NA, NA),
                            c_2017 <- c(NA, NA, NA),
                            c_2018 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_6[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.05, x100 = -10, theta = 1) # Starting parameters
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_bbinom$par # Store optimized parameters
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

# Create the plot
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution 
       Type of data: negative clairvoyance data", 
    caption = "Figure A25. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a CNE function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function/binomial distribution
##
coefs_hill_binom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(years))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(years)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[6]]['a']*(1 - (x^(coefs_hill_binom[[6]]['n'])/(coefs_hill_binom[[6]]['h']^(coefs_hill_binom[[6]]['n']) + x^(coefs_hill_binom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution 
       Type of data: negative clairvoyance data", 
    caption = "Figure A26. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a Hill function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))

##
## Hill function/beta-binomial distribution
##
coefs_hill_bbinom <- list(c_2013 <- c(NA, NA, NA, NA),
                          c_2014 <- c(NA, NA, NA, NA),
                          c_2015 <- c(NA, NA, NA, NA),
                          c_2016 <- c(NA, NA, NA, NA),
                          c_2017 <- c(NA, NA, NA, NA),
                          c_2018 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(years))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(years)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[6]]['a']*(1 - (x^(coefs_hill_bbinom[[6]]['n'])/(coefs_hill_bbinom[[6]]['h']^(coefs_hill_bbinom[[6]]['n']) + x^(coefs_hill_bbinom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution; 
       Type of data: negative clairvoyance data", 
    caption = "Figure A27. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a Hill function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


####################################################
# Negative clairvoyance & positive carry-over data #
####################################################

#
# Data morphing and distance circle creation
#

dat_7 <- dat

# Make data clairvoyant
for(i in length(years):2){ 
  co <- which(dat_7[[i]]$result == 0)
  dat_7[[i-1]] <- rbind(dat_7[[i-1]], dat_7[[i]][co,])
  co_2 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & 
                  dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & 
                  dat_7[[i-1]]$year == (2012 + i))
  if(length(co_2) >= 1){
    for(j in 1:length(co_2)){
      if(dat_7[[i-1]]$result[co_2[j]] == 1){ 
        co_3 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & 
                        dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & 
                        dat_7[[i-1]]$year == (2012 + i))
        dat_7[[i-1]] <- dat_7[[i-1]][-co_3[j],]
      }
      if(dat_7[[i-1]]$result[co_2[j]] == 0){ 
        dat_7[[i-1]] <- dat_7[[i-1]][-co_2[j],]
      }
    }
  } 
}

# make positives carry over
for(i in 2:length(years)){
  co <- which(dat_7[[i-1]]$result == 1)
  dat_7[[i]] <- rbind(dat_7[[i]], dat_7[[i-1]][co,])
  co_2 <- which(dat_7[[i]]$lon == dat_7[[i-1]]$lon[co] & 
                  dat_7[[i]]$lat == dat_7[[i-1]]$lat[co] & 
                  dat_7[[i]]$year == (2012 + i))
  
  if(length(co_2) >= 1){ 
    dat_7[[i]] <- dat_7[[i]][-co_2,] 
  }
}

# Set distance circles
dc <- 0:(max_dist - 1)
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(years))
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(years))

dat_8 <- list()

for(i in 1:length(years)){
  for(j in dc){
    pos[[i]][j + 1] <- sum((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j+1)]))
    n[[i]][j + 1] <- length((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j+1)]))
  }
  dat_8[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n)
}

# Plot of the data
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, y = prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  labs(subtitle = "Negative clairvoyance & positive carry-over data",
    caption = "Figure A28. A graph with the datapoints of the negative clairvoyance & positive 
       carry-over data aggregated with distance circles.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(plot.caption = element_text(hjust = 0))


#
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames.
coefs_nexp_binom <- list(c_2013 <- c(NA, NA),
                         c_2014 <- c(NA, NA),
                         c_2015 <- c(NA, NA),
                         c_2016 <- c(NA, NA),
                         c_2017 <- c(NA, NA),
                         c_2018 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1)
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n)
  c_opt <- opt_nexp_binom$par
  coefs_nexp_binom[[i]] <- c_opt
  comb_AIC_nexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_nexp_binom$value)
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[6]]['a']*exp(-coefs_nexp_binom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution; 
       Type of data: negative clairvoyance & positive carry-over data",
    caption = 
      "Figure A29. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a negative exponential function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", 
                                               "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                          c_2014 <- c(NA, NA, NA),
                          c_2015 <- c(NA, NA, NA),
                          c_2016 <- c(NA, NA, NA),
                          c_2017 <- c(NA, NA, NA),
                          c_2018 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(years))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(years)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[6]]['a']*exp(-coefs_nexp_bbinom[[6]]['b']*x), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution; 
       Type of data: negative clairvoyance & positive carry-over data", 
    caption = "Figure A30. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a negative exponential function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_2013 <- c(NA, NA),
                        c_2014 <- c(NA, NA),
                        c_2015 <- c(NA, NA),
                        c_2016 <- c(NA, NA),
                        c_2017 <- c(NA, NA),
                        c_2018 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(years))

fun_log_binom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['b'])))))
  # mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, b = -30)

for(i in 1:length(years)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*length(pars) + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_binom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[1]]['a'] * (x - coefs_log_binom[[1]]['b'])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[2]]['a'] * (x - coefs_log_binom[[2]]['b'])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[3]]['a'] * (x - coefs_log_binom[[3]]['b'])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[4]]['a'] * (x - coefs_log_binom[[4]]['b'])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[5]]['a'] * (x - coefs_log_binom[[5]]['b'])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[6]]['a'] * (x - coefs_log_binom[[6]]['b'])))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution; 
       Type of data: negative clairvoyance & positive carry-over data", 
    caption = "Figure A31. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a logistic function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(years))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['c'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, c = 30, theta = 5)

for(i in 1:length(years)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = 2013:2018, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = c(0.8, 0.75, 0.7, 0.65, 0.6, 0.55))
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[1]]['a'] * (x - coefs_log_bbinom[[1]]['c'])))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[2]]['a'] * (x - coefs_log_bbinom[[2]]['c'])))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[3]]['a'] * (x - coefs_log_bbinom[[3]]['c'])))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[4]]['a'] * (x - coefs_log_bbinom[[4]]['c'])))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[5]]['a'] * (x - coefs_log_bbinom[[5]]['c'])))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[6]]['a'] * (x - coefs_log_bbinom[[6]]['c'])))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution; 
       Type of data: negative clairvoyance & positive carry-over data",
    caption = "Figure A32. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a logistic function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_binom <- list(c_2013 <- c(NA, NA),
                           c_2014 <- c(NA, NA),
                           c_2015 <- c(NA, NA),
                           c_2016 <- c(NA, NA),
                           c_2017 <- c(NA, NA),
                           c_2018 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_8[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.1, x100 = -10) # Starting parameters
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_binom$par # Store optimized parameters
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[6]]['b']*coefs_conexp_binom[[6]]['x100']))*exp(-coefs_conexp_binom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution; 
       Type of data: negative clairvoyance & positive carry-over data", 
    caption = "Figure A33. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a CNE function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_bbinom <- list(c_2013 <- c(NA, NA, NA),
                            c_2014 <- c(NA, NA, NA),
                            c_2015 <- c(NA, NA, NA),
                            c_2016 <- c(NA, NA, NA),
                            c_2017 <- c(NA, NA, NA),
                            c_2018 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(years))

## Setup model function
fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/exp(-par['b']*par['x100']))*exp(-par['b']*x) < 0.999, (1/exp(-par['b']*par['x100']))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(years)){ # For-loop through every year
  df <- dat_8[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.05, x100 = -10, theta = 1) # Starting parameters
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_bbinom$par # Store optimized parameters
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(years))
for(i in 1:length(years)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = years, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

# Create the plot
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[6]]['b']*coefs_conexp_bbinom[[6]]['x100']))*exp(-coefs_conexp_bbinom[[6]]['b']*x), 0.999), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution; 
       Type of data: negative clairvoyance & positive carry-over data", 
    caption = "Figure A34. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a CNE function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function/binomial distribution
##
coefs_hill_binom <- list(c_2013 <- c(NA, NA, NA),
                         c_2014 <- c(NA, NA, NA),
                         c_2015 <- c(NA, NA, NA),
                         c_2016 <- c(NA, NA, NA),
                         c_2017 <- c(NA, NA, NA),
                         c_2018 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(years))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(years)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[6]]['a']*(1 - (x^(coefs_hill_binom[[6]]['n'])/(coefs_hill_binom[[6]]['h']^(coefs_hill_binom[[6]]['n']) + x^(coefs_hill_binom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Binomial distribution; 
       Type of data: negative clairvoyance data & positive carry-over",
    caption = "Figure A35. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a Hill function with binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))

##
## Hill function/beta-binomial distribution
##
coefs_hill_bbinom <- list(c_2013 <- c(NA, NA, NA, NA),
                          c_2014 <- c(NA, NA, NA, NA),
                          c_2015 <- c(NA, NA, NA, NA),
                          c_2016 <- c(NA, NA, NA, NA),
                          c_2017 <- c(NA, NA, NA, NA),
                          c_2018 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(years))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(years)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = years, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[6]], aes(x = dist, prop, color = "2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[6]]['a']*(1 - (x^(coefs_hill_bbinom[[6]]['n'])/(coefs_hill_bbinom[[6]]['h']^(coefs_hill_bbinom[[6]]['n']) + x^(coefs_hill_bbinom[[6]]['n']))))), 
                aes(color = "2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution; 
       Type of data: negative clairvoyance & positive carry-over data", 
    caption = "Figure A36. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a Hill function with beta-binomial distribution.",
    x = "Distance", 
    y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


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

# Read in data
dat <- dat_base

# Adjust for date, NA's, and missing coordinates
dat <- dat %>%
  mutate(date = as.Date(date)) %>%
  drop_na() %>%
  filter((lon != 0) | (lat != 0))

# Set seasons. Month 3 is March, month 4 is April.
dat <- dat %>% 
  mutate(season = 
           ifelse((dat$year == 2013) | (dat$year == 2014 & dat$month %in% 1:3), 1, 
                  ifelse((dat$year == 2014 & dat$month %in% 4:12) | 
                           (dat$year == 2015 & dat$month %in% 1:3), 2, 
                         ifelse((dat$year == 2015 & dat$month %in% 4:12) | 
                                  (dat$year == 2016 & dat$month %in% 1:3), 3, 
                                ifelse((dat$year == 2016 & dat$month %in% 4:12) | 
                                         (dat$year == 2017 & dat$month %in% 1:3), 4, 5)))))


# Set number of seasons, used for length of loops and such
seasons <- unique(dat$season)

# Calculate distance of data points to origin
dat <- dat %>%
  mutate(dist = round((distHaversine(tibble(lon = dat$lon, lat = dat$lat), 
                                     c(gallipoli[1], gallipoli[2])))/1000))

# Find the maximum distance measured for the distance circels and the maximum distance 
# where a positive is found for visualization
max_dist <- max(dat$dist)
max_pos <- max(dat[dat$result == 1,]$dist)

# Split the data for each season
dat <- dat %>% 
  split(dat$season)

names(dat) <- c("2013/2014", "2014/2015", "2015/2016", "2016/2017", "2017/2018")

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
    pos[[i]][j + 1] <- sum((dat[[i]]$result[dat[[i]]$dist >= j & (dat[[i]]$dist < j+1)])) 
    n[[i]][j + 1] <- length((dat[[i]]$result[dat[[i]]$dist >= j & (dat[[i]]$dist < j+1)])) 
  }
  dat_2[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n) 
}

# Plot of the data
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, y = prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  labs(x = "Distance", 
       y = "Proportion positive",
       caption = "Figure A37. A graph with the datapoints of the original 'season' data aggregated with distance circles.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange")) +
  theme(plot.caption = element_text(hjust = 0))


#
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_nexp_binom <- list(c_1314 <- c(NA, NA),
                         c_1415 <- c(NA, NA),
                         c_1516 <- c(NA, NA),
                         c_1617 <- c(NA, NA),
                         c_1718 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(seasons))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){ # For-loop through every year
  df <- dat_2[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(a = 0.1, b = 0.1) # Starting parameters
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n) # Optimization
  c_opt <- opt_nexp_binom$par # Store optimized parameters
  coefs_nexp_binom[[i]] <- c_opt 
  comb_AIC_nexp_binom[i] <- 2*2 + 2*as.numeric(opt_nexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = 
      "Figure A38. A graph with the datapoints of the original data and lines as fitted with a negative 
       exponential function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", 
                                               "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                          c_1415 <- c(NA, NA, NA),
                          c_1516 <- c(NA, NA, NA),
                          c_1617 <- c(NA, NA, NA),
                          c_1718 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(seasons))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(seasons)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*3 + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = tjust)
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A39. A graph with the datapoints of the original data and lines as fitted with a negative 
       exponential function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_1314 <- c(NA, NA),
                        c_1415 <- c(NA, NA),
                        c_1516 <- c(NA, NA),
                        c_1617 <- c(NA, NA),
                        c_1718 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(seasons))

fun_log_binom <- function(par, x, z, n){
  mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 0.1, b = 0.1)

for(i in 1:length(seasons)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*2 + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_binom, digits = 0), tjust = tjust)
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(exp(coefs_log_binom[[1]]['a'] + coefs_log_binom[[1]]['b']*x))/(1 + exp(coefs_log_binom[[1]]['a'] + coefs_log_binom[[1]]['b']*x)), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[2]]['a'] + coefs_log_binom[[2]]['b']*x))/(1 + exp(coefs_log_binom[[2]]['a'] + coefs_log_binom[[2]]['b']*x)), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[3]]['a'] + coefs_log_binom[[3]]['b']*x))/(1 + exp(coefs_log_binom[[3]]['a'] + coefs_log_binom[[3]]['b']*x)), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[4]]['a'] + coefs_log_binom[[4]]['b']*x))/(1 + exp(coefs_log_binom[[4]]['a'] + coefs_log_binom[[4]]['b']*x)), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_binom[[5]]['a'] + coefs_log_binom[[5]]['b']*x))/(1 + exp(coefs_log_binom[[5]]['a'] + coefs_log_binom[[5]]['b']*x)), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A40. A graph with the datapoints of the original data and lines as fitted with a logistic 
       function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(seasons))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 0.1, b = 0.1, theta = 1)

for(i in 1:length(seasons)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n)
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*3 + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = tjust)
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[1]]['a'] + coefs_log_bbinom[[1]]['b']*x))/(1 + exp(coefs_log_bbinom[[1]]['a'] + coefs_log_bbinom[[1]]['b']*x)), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[2]]['a'] + coefs_log_bbinom[[2]]['b']*x))/(1 + exp(coefs_log_bbinom[[2]]['a'] + coefs_log_bbinom[[2]]['b']*x)), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[3]]['a'] + coefs_log_bbinom[[3]]['b']*x))/(1 + exp(coefs_log_bbinom[[3]]['a'] + coefs_log_bbinom[[3]]['b']*x)), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[4]]['a'] + coefs_log_bbinom[[4]]['b']*x))/(1 + exp(coefs_log_bbinom[[4]]['a'] + coefs_log_bbinom[[4]]['b']*x)), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(exp(coefs_log_bbinom[[5]]['a'] + coefs_log_bbinom[[5]]['b']*x))/(1 + exp(coefs_log_bbinom[[5]]['a'] + coefs_log_bbinom[[5]]['b']*x)), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A41. A graph with the datapoints of the original data and lines as fitted with a logistic 
       function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_binom <- list(c_1314 <- c(NA, NA),
                           c_1415 <- c(NA, NA),
                           c_1516 <- c(NA, NA),
                           c_1617 <- c(NA, NA),
                           c_1718 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(seasons))

## Setup model function
fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){ # For-loop through every year
  df <- dat_2[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.1, x100 = -10) # Starting parameters
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_binom$par # Store optimized parameters
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A42. A graph with the datapoints of the original data and lines as fitted with a CNE 
       function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                            c_1415 <- c(NA, NA, NA),
                            c_1516 <- c(NA, NA, NA),
                            c_1617 <- c(NA, NA, NA),
                            c_1718 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(seasons))

## Setup model function
fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){ # For-loop through every year
  df <- dat_2[[i]] # Data frame of the current year with distance, positives, and number of measurements of every dc that year
  x <- df$dist # Distance
  z <- df$pos # Positives
  n <- df$n # Number of measurements
  
  pars <- c(b = 0.1, x100 = -10, theta = 1) # Starting parameters
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_bbinom$par # Store optimized parameters
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

# Create the plot
ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A43. A graph with the datapoints of the original data and lines as fitted with a CNE 
       function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, binomial distribution
##

coefs_hill_binom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(seasons))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(seasons)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A44. A graph with the datapoints of the original data and lines as fitted with a Hill 
       function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, beta-binomial distribution
##

coefs_hill_bbinom <- list(c_1314 <- c(NA, NA, NA, NA),
                          c_1415 <- c(NA, NA, NA, NA),
                          c_1516 <- c(NA, NA, NA, NA),
                          c_1617 <- c(NA, NA, NA, NA),
                          c_1718 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(seasons))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(seasons)){
  t <- dat_2[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_2[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_2[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution;
       Type of data: original data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A45. A graph with the datapoints of the original data and lines as fitted with a Hill 
       function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


############################
# Positive carry-over data #
############################

#
# Data morphing and distance circle creation
#

dat_3 <- dat

# make positives carry over
for(i in 2:length(seasons)){ 
  co <- which(dat_3[[i-1]]$result == 1) 
  dat_3[[i]] <- rbind(dat_3[[i]], dat_3[[i-1]][co,]) 
  co_2 <- which(dat_3[[i]]$lon == dat_3[[i-1]]$lon[co] & dat_3[[i]]$lat == dat_3[[i-1]]$lat[co] & dat_3[[i]]$year == (2012 + i)) 
  if(length(co_2) >= 1){ 
    dat_3[[i]] <- dat_3[[i]][-co_2,] 
  }
}

# Set distance circles
dc <- 0:(max_dist - 1)
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))

dat_4 <- list()

for(i in 1:length(seasons)){
  for(j in dc){
    pos[[i]][j + 1] <- sum((dat_3[[i]]$result[dat_3[[i]]$dist >= j & (dat_3[[i]]$dist < j+1)]))
    n[[i]][j + 1] <- length((dat_3[[i]]$result[dat_3[[i]]$dist >= j & (dat_3[[i]]$dist < j+1)]))
  }
  dat_4[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n)
}

# Plot of the data
ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, y = prop - 0.005, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop - 0.01, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop - 0.015, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop - 0.02, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop - 0.025, color = "2017/2018"), shape = 1, position = "jitter") +
  labs(x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A46. A graph with the datapoints of the positive carry-over 'season' data aggregated with 
       distance circles.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(plot.caption = element_text(hjust = 0))

#
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames.
coefs_nexp_binom <- list(c_1314 <- c(NA, NA),
                         c_1415 <- c(NA, NA),
                         c_1516 <- c(NA, NA),
                         c_1617 <- c(NA, NA),
                         c_1718 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(seasons))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1)
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n)
  c_opt <- opt_nexp_binom$par
  coefs_nexp_binom[[i]] <- c_opt
  comb_AIC_nexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_nexp_binom$value)
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A47. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       negative exponential function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                          c_1415 <- c(NA, NA, NA),
                          c_1516 <- c(NA, NA, NA),
                          c_1617 <- c(NA, NA, NA),
                          c_1718 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(seasons))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(seasons)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*3 + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = tjust)
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A48. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       negative exponential function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_1314 <- c(NA, NA),
                        c_1415 <- c(NA, NA),
                        c_1516 <- c(NA, NA),
                        c_1617 <- c(NA, NA),
                        c_1718 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(seasons))

fun_log_binom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['b'])))))
  # mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, b = -30)

for(i in 1:length(seasons)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*2 + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_binom, digits = 0), tjust = tjust)
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[1]]['a'] * (x - coefs_log_binom[[1]]['b'])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[2]]['a'] * (x - coefs_log_binom[[2]]['b'])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[3]]['a'] * (x - coefs_log_binom[[3]]['b'])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[4]]['a'] * (x - coefs_log_binom[[4]]['b'])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[5]]['a'] * (x - coefs_log_binom[[5]]['b'])))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A49. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       logistic function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(seasons))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['c'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, c = 30, theta = 5)

for(i in 1:length(seasons)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*3 + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = tjust)
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[1]]['a'] * (x - coefs_log_bbinom[[1]]['c'])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[2]]['a'] * (x - coefs_log_bbinom[[2]]['c'])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[3]]['a'] * (x - coefs_log_bbinom[[3]]['c'])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[4]]['a'] * (x - coefs_log_bbinom[[4]]['c'])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[5]]['a'] * (x - coefs_log_bbinom[[5]]['c'])))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A50. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       logistic function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, binomial distribution
##

coefs_conexp_binom <- list(c_1314 <- c(NA, NA),
                           c_1415 <- c(NA, NA),
                           c_1516 <- c(NA, NA),
                           c_1617 <- c(NA, NA),
                           c_1718 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(seasons))

fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

for(i in 1:length(seasons)){ 
  df <- dat_4[[i]] 
  x <- df$dist 
  z <- df$pos 
  n <- df$n 
  
  pars <- c(b = 0.1, x100 = -10) 
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) 
  c_opt <- opt_conexp_binom$par 
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) 
}

tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A51. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       CNE function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

coefs_conexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                            c_1415 <- c(NA, NA, NA),
                            c_1516 <- c(NA, NA, NA),
                            c_1617 <- c(NA, NA, NA),
                            c_1718 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(seasons))

fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(seasons)){ 
  df <- dat_4[[i]] 
  x <- df$dist 
  z <- df$pos 
  n <- df$n 
  
  pars <- c(b = 0.1, x100 = -10, theta = 1) 
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) 
  c_opt <- opt_conexp_bbinom$par 
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value) 
}

tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A52. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       CNE function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, binomial distribution
##

coefs_hill_binom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(seasons))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(seasons)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A53. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       Hill function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, beta-binomial distribution
##

coefs_hill_bbinom <- list(c_1314 <- c(NA, NA, NA, NA),
                          c_1415 <- c(NA, NA, NA, NA),
                          c_1516 <- c(NA, NA, NA, NA),
                          c_1617 <- c(NA, NA, NA, NA),
                          c_1718 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(seasons))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(seasons)){
  t <- dat_4[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_4[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_4[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs(
    subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution;
       Type of data: positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A54. A graph with the datapoints of the positive carry-over data and lines as fitted with a 
       Hill function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##############################
# Negative clairvoyance data #
##############################

#
# Data morphing and distance circle creation
#

dat_5 <- dat

# Make data clairvoyant
for(i in length(seasons):2){ 
  co <- which(dat_5[[i]]$result == 0) 
  dat_5[[i-1]] <- rbind(dat_5[[i-1]], dat_5[[i]][co,]) 
  co_2 <- which(dat_5[[i-1]]$lon == dat_5[[i]]$lon[co] & dat_5[[i-1]]$lat == dat_5[[i]]$lat[co] & dat_5[[i-1]]$year == (2012 + i)) # Check for multiple samples at the same location
  if(length(co_2) >= 1){ 
    for(j in 1:length(co_2)){ 
      if(dat_5[[i-1]]$result[co_2[j]] == 1){ 
        co_3 <- which(dat_5[[i-1]]$lon == dat_5[[i]]$lon[co] & dat_5[[i-1]]$lat == dat_5[[i]]$lat[co] & dat_5[[i-1]]$year == (2012 + i))
        dat_5[[i-1]] <- dat_5[[i-1]][-co_3[j],]
      }
      if(dat_5[[i-1]]$result[co_2[j]] == 0){ 
        dat_5[[i-1]] <- dat_5[[i-1]][-co_2[j],]
      }
    }
  } 
}

# Set distance circles
dc <- 0:(max_dist - 1)
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))

dat_6 <- list()

for(i in 1:length(seasons)){
  for(j in dc){
    pos[[i]][j + 1] <- sum((dat_5[[i]]$result[dat_5[[i]]$dist >= j & (dat_5[[i]]$dist < j+1)]))
    n[[i]][j + 1] <- length((dat_5[[i]]$result[dat_5[[i]]$dist >= j & (dat_5[[i]]$dist < j+1)]))
  }
  dat_6[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n)
}

# Plot of the data
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, y = prop - 0.005, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop - 0.01, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop - 0.015, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop - 0.02, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop - 0.025, color = "2017/2018"), shape = 1, position = "jitter") +
  labs(x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A55. A graph with the datapoints of the negative clairvoyance 'season' data aggregated 
       with distance circles.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(plot.caption = element_text(hjust = 0))


#
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames.
coefs_nexp_binom <- list(c_1314 <- c(NA, NA),
                         c_1415 <- c(NA, NA),
                         c_1516 <- c(NA, NA),
                         c_1617 <- c(NA, NA),
                         c_1718 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(seasons))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1)
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n)
  c_opt <- opt_nexp_binom$par
  coefs_nexp_binom[[i]] <- c_opt
  comb_AIC_nexp_binom[i] <- 2*2 + 2*as.numeric(opt_nexp_binom$value)
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A56. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a negative exponential function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                          c_1415 <- c(NA, NA, NA),
                          c_1516 <- c(NA, NA, NA),
                          c_1617 <- c(NA, NA, NA),
                          c_1718 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(seasons))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(seasons)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*3 + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = tjust)
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A57. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a negative exponential function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_1314 <- c(NA, NA),
                        c_1415 <- c(NA, NA),
                        c_1516 <- c(NA, NA),
                        c_1617 <- c(NA, NA),
                        c_1718 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(seasons))

fun_log_binom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['b'])))))
  # mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, b = -30)

for(i in 1:length(seasons)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*2 + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_binom, digits = 0), tjust = tjust)
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[1]]['a'] * (x - coefs_log_binom[[1]]['b'])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[2]]['a'] * (x - coefs_log_binom[[2]]['b'])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[3]]['a'] * (x - coefs_log_binom[[3]]['b'])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[4]]['a'] * (x - coefs_log_binom[[4]]['b'])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[5]]['a'] * (x - coefs_log_binom[[5]]['b'])))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A58. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a logistic function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(seasons))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['c'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, c = 30, theta = 5)

for(i in 1:length(seasons)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = tjust)
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[1]]['a'] * (x - coefs_log_bbinom[[1]]['c'])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[2]]['a'] * (x - coefs_log_bbinom[[2]]['c'])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[3]]['a'] * (x - coefs_log_bbinom[[3]]['c'])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[4]]['a'] * (x - coefs_log_bbinom[[4]]['c'])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[5]]['a'] * (x - coefs_log_bbinom[[5]]['c'])))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(
    subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A59. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a logistic function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_binom <- list(c_1314 <- c(NA, NA),
                           c_1415 <- c(NA, NA),
                           c_1516 <- c(NA, NA),
                           c_1617 <- c(NA, NA),
                           c_1718 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(seasons))

## Setup model function
fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){ 
  df <- dat_6[[i]] 
  x <- df$dist 
  z <- df$pos 
  n <- df$n 
  
  pars <- c(b = 0.1, x100 = -10)
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) 
  c_opt <- opt_conexp_binom$par 
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) 
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A60. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a negative CNE with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                            c_1415 <- c(NA, NA, NA),
                            c_1516 <- c(NA, NA, NA),
                            c_1617 <- c(NA, NA, NA),
                            c_1718 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(seasons))

## Setup model function
fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){ 
  df <- dat_6[[i]] 
  x <- df$dist 
  z <- df$pos 
  n <- df$n 
  
  pars <- c(b = 0.1, x100 = -10, theta = 1)
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) 
  c_opt <- opt_conexp_bbinom$par 
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value) 
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

# Create the plot
ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A61. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a CNE function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, binomial distribution
##

coefs_hill_binom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(seasons))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(seasons)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A62. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a Hill function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, beta-binomial distribution
##

coefs_hill_bbinom <- list(c_1314 <- c(NA, NA, NA, NA),
                          c_1415 <- c(NA, NA, NA, NA),
                          c_1516 <- c(NA, NA, NA, NA),
                          c_1617 <- c(NA, NA, NA, NA),
                          c_1718 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(seasons))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(seasons)){
  t <- dat_6[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_6[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_6[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A63. A graph with the datapoints of the negative clairvoyance data and lines as fitted with 
       a Hill function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


######################################################
# Negative clairvoyance and positive carry-over data #
######################################################

#
# Data morphing and distance circle creation
#

dat_7 <- dat

# Make data clairvoyant
for(i in length(seasons):2){ 
  co <- which(dat_7[[i]]$result == 0) 
  dat_7[[i-1]] <- rbind(dat_7[[i-1]], dat_7[[i]][co,]) 
  co_2 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & dat_7[[i-1]]$year == (2012 + i)) 
  if(length(co_2) >= 1){ 
    for(j in 1:length(co_2)){ 
      if(dat_7[[i-1]]$result[co_2[j]] == 1){ 
        co_3 <- which(dat_7[[i-1]]$lon == dat_7[[i]]$lon[co] & dat_7[[i-1]]$lat == dat_7[[i]]$lat[co] & dat_7[[i-1]]$year == (2012 + i))
        dat_7[[i-1]] <- dat_7[[i-1]][-co_3[j],]
      }
      if(dat_7[[i-1]]$result[co_2[j]] == 0){ 
        dat_7[[i-1]] <- dat_7[[i-1]][-co_2[j],]
      }
    }
  } 
}

# make positives carry over
for(i in 2:length(seasons)){ 
  co <- which(dat_7[[i-1]]$result == 1) 
  dat_7[[i]] <- rbind(dat_7[[i]], dat_7[[i-1]][co,])
  co_2 <- which(dat_7[[i]]$lon == dat_7[[i-1]]$lon[co] & dat_7[[i]]$lat == dat_7[[i-1]]$lat[co] & dat_7[[i]]$year == (2012 + i)) 
  # print(co_2)
  if(length(co_2) >= 1){ 
    dat_7[[i]] <- dat_7[[i]][-co_2,] 
  }
}

# Set distance circles
dc <- 0:(max_dist - 1)
pos <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))
n <- rep(list(rep(NA, times = length(dc + 1))), times = length(seasons))

dat_8 <- list()

for(i in 1:length(seasons)){
  for(j in dc){
    pos[[i]][j + 1] <- sum((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j+1)]))
    n[[i]][j + 1] <- length((dat_7[[i]]$result[dat_7[[i]]$dist >= j & (dat_7[[i]]$dist < j+1)]))
  }
  dat_8[[i]] <- tibble(dist = 1:max_dist, pos = pos[[i]], n = n[[i]]) %>% 
    mutate(prop = pos/n)
}

# Plot of the data
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, y = prop - 0.005, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop - 0.01, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop - 0.015, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop - 0.02, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop - 0.025, color = "2017/2018"), shape = 1, position = "jitter") +
  labs(x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A64. A graph with the datapoints of the negative clairvoyance & positive carry-over 
       'season' data aggregated with distance circles.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(plot.caption = element_text(hjust = 0))

#
# Model fitting
#

##
## Negative exponential, binomial distribution
##

## Setup parameter value and AIC data frames.
coefs_nexp_binom <- list(c_1314 <- c(NA, NA),
                         c_1415 <- c(NA, NA),
                         c_1516 <- c(NA, NA),
                         c_1617 <- c(NA, NA),
                         c_1718 <- c(NA, NA))

comb_AIC_nexp_binom <- rep(NA, times = length(seasons))

## Setup model function
fun_nexp_binom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1)
  opt_nexp_binom <- optim(par = pars, fn = fun_nexp_binom, x = x, z = z, n = n)
  c_opt <- opt_nexp_binom$par
  coefs_nexp_binom[[i]] <- c_opt
  comb_AIC_nexp_binom[i] <- 2*2 + 2*as.numeric(opt_nexp_binom$value)
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_nexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_nexp_binom, digits = 0), 
             tjust)
AIC_nexp_binom_text$label <- 
  paste(AIC_nexp_binom_text$year, ": ", AIC_nexp_binom_text$AIC)

total_AIC_nexp_binom <- sum(comb_AIC_nexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[1]]['a']*exp(-coefs_nexp_binom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[2]]['a']*exp(-coefs_nexp_binom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[3]]['a']*exp(-coefs_nexp_binom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[4]]['a']*exp(-coefs_nexp_binom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)coefs_nexp_binom[[5]]['a']*exp(-coefs_nexp_binom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A65. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a negative exponential function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Negative exponential, beta-binomial distribution
##

coefs_nexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                          c_1415 <- c(NA, NA, NA),
                          c_1516 <- c(NA, NA, NA),
                          c_1617 <- c(NA, NA, NA),
                          c_1718 <- c(NA, NA, NA))

comb_AIC_nexp_bbinom <- rep(NA, times = length(seasons))

fun_nexp_bbinom <- function(par, x, z, n){
  mu <- par['a']*exp(-par['b']*x)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

for(i in 1:length(seasons)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  pars <- c(a = 0.1, b = 0.1, theta = 1)
  opt_nexp_bbinom <- optim(par = pars, fn = fun_nexp_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_nexp_bbinom$par
  coefs_nexp_bbinom[[i]] <- c_opt
  comb_AIC_nexp_bbinom[i] <- 2*3 + 2*as.numeric(opt_nexp_bbinom$value)
}

AIC_nexp_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_nexp_bbinom, digits = 0), tjust = tjust)
AIC_nexp_bbinom_text$label <- paste(AIC_nexp_bbinom_text$year, ": ", AIC_nexp_bbinom_text$AIC)

total_AIC_nexp_bbinom <- sum(comb_AIC_nexp_bbinom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)coefs_nexp_bbinom[[1]]['a']*exp(-coefs_nexp_bbinom[[1]]['b']*x), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[2]]['a']*exp(-coefs_nexp_bbinom[[2]]['b']*x), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[3]]['a']*exp(-coefs_nexp_bbinom[[3]]['b']*x), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[4]]['a']*exp(-coefs_nexp_bbinom[[4]]['b']*x), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)coefs_nexp_bbinom[[5]]['a']*exp(-coefs_nexp_bbinom[[5]]['b']*x), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A66. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a negative exponential function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, binomial distribution
##

coefs_log_binom <- list(c_1314 <- c(NA, NA),
                        c_1415 <- c(NA, NA),
                        c_1516 <- c(NA, NA),
                        c_1617 <- c(NA, NA),
                        c_1718 <- c(NA, NA))

comb_AIC_log_binom <- rep(NA, times = length(seasons))

fun_log_binom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['b'])))))
  # mu <- (exp(par['a'] + par['b']*x))/(1 + exp(par['a'] + par['b']*x))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, b = -30)

for(i in 1:length(seasons)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_binom <- optim(par = pars, fn = fun_log_binom, x = x, z = z, n = n)
  c_opt <- opt_log_binom$par
  coefs_log_binom[[i]] <- c_opt
  comb_AIC_log_binom[i] <- 2*2 + 2*as.numeric(opt_log_binom$value)
}

AIC_log_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_binom, digits = 0), tjust = tjust)
AIC_log_binom_text$label <- paste(AIC_log_binom_text$year, ": ", AIC_log_binom_text$AIC)

total_AIC_log_binom <- sum(comb_AIC_log_binom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[1]]['a'] * (x - coefs_log_binom[[1]]['b'])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[2]]['a'] * (x - coefs_log_binom[[2]]['b'])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[3]]['a'] * (x - coefs_log_binom[[3]]['b'])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[4]]['a'] * (x - coefs_log_binom[[4]]['b'])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_binom[[5]]['a'] * (x - coefs_log_binom[[5]]['b'])))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A67. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a logistic function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Logistic function, beta-binomial distribution
##

coefs_log_bbinom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_log_bbinom <- rep(NA, times = length(seasons))

fun_log_bbinom <- function(par, x, z, n){
  mu <- (1 / (1 + exp(-(par['a'] * (x - par['c'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = -0.1, c = 30, theta = 5)

for(i in 1:length(seasons)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_log_bbinom <- optim(par = pars, fn = fun_log_bbinom, x = x, z = z, n = n, control = list(maxit = 10000))
  c_opt <- opt_log_bbinom$par
  coefs_log_bbinom[[i]] <- c_opt
  comb_AIC_log_bbinom[i] <- 2*3 + 2*as.numeric(opt_log_bbinom$value)
}

AIC_log_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_log_bbinom, digits = 0), tjust = tjust)
AIC_log_bbinom_text$label <- paste(AIC_log_bbinom_text$year, ": ", AIC_log_bbinom_text$AIC)

total_AIC_log_bbinom <- sum(comb_AIC_log_bbinom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[1]]['a'] * (x - coefs_log_bbinom[[1]]['c'])))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[2]]['a'] * (x - coefs_log_bbinom[[2]]['c'])))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[3]]['a'] * (x - coefs_log_bbinom[[3]]['c'])))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[4]]['a'] * (x - coefs_log_bbinom[[4]]['c'])))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(1 / (1 + exp(-coefs_log_bbinom[[5]]['a'] * (x - coefs_log_bbinom[[5]]['c'])))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A68. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a logistic function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_binom <- list(c_1314 <- c(NA, NA),
                           c_1415 <- c(NA, NA),
                           c_1516 <- c(NA, NA),
                           c_1617 <- c(NA, NA),
                           c_1718 <- c(NA, NA))

comb_AIC_conexp_binom <- rep(NA, times = length(seasons))

## Setup model function
fun_conexp_binom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){ 
  df <- dat_8[[i]] 
  x <- df$dist 
  z <- df$pos
  n <- df$n 
  
  pars <- c(b = 0.1, x100 = -10) # Starting parameters
  opt_conexp_binom <- optim(par = pars, fn = fun_conexp_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) # Optimization
  c_opt <- opt_conexp_binom$par # Store optimized parameters
  coefs_conexp_binom[[i]] <- c_opt 
  comb_AIC_conexp_binom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_binom$value) # Calculate and store AIC
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_binom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_binom, digits = 0), 
             tjust)
AIC_conexp_binom_text$label <- 
  paste(AIC_conexp_binom_text$year, ": ", AIC_conexp_binom_text$AIC)

total_AIC_conexp_binom <- sum(comb_AIC_conexp_binom)

# Create the plot
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[1]]['b']*coefs_conexp_binom[[1]]['x100']))*exp(-coefs_conexp_binom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[2]]['b']*coefs_conexp_binom[[2]]['x100']))*exp(-coefs_conexp_binom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[3]]['b']*coefs_conexp_binom[[3]]['x100']))*exp(-coefs_conexp_binom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[4]]['b']*coefs_conexp_binom[[4]]['x100']))*exp(-coefs_conexp_binom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_binom[[5]]['b']*coefs_conexp_binom[[5]]['x100']))*exp(-coefs_conexp_binom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A69. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a CNE function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## CNE, Beta-binomial distribution
##

## Setup parameter value and AIC data frames to view later.
coefs_conexp_bbinom <- list(c_1314 <- c(NA, NA, NA),
                            c_1415 <- c(NA, NA, NA),
                            c_1516 <- c(NA, NA, NA),
                            c_1617 <- c(NA, NA, NA),
                            c_1718 <- c(NA, NA, NA))

comb_AIC_conexp_bbinom <- rep(NA, times = length(seasons))

## Setup model function
fun_conexp_bbinom <- function(par, x, z, n){
  mu <- ifelse((1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x) < 0.999, (1/(exp(-par['b']*par['x100'])))*exp(-par['b']*x), 0.999)
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}

## Run model fitting for every year
for(i in 1:length(seasons)){ 
  df <- dat_8[[i]] 
  x <- df$dist 
  z <- df$pos 
  n <- df$n 
  
  pars <- c(b = 0.1, x100 = -10, theta = 1)
  opt_conexp_bbinom <- optim(par = pars, fn = fun_conexp_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 2000)) 
  c_opt <- opt_conexp_bbinom$par
  coefs_conexp_bbinom[[i]] <- c_opt 
  comb_AIC_conexp_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_conexp_bbinom$value)
}

## Prepare AIC text for graph
tjust <- rep(NA, times = length(seasons))
for(i in 1:length(seasons)){
  tjust[i] = 0.8 - (0.05*(i - 1))
}
AIC_conexp_bbinom_text <- 
  data.frame(year = seasons, 
             AIC = round(comb_AIC_conexp_bbinom, digits = 0), 
             tjust)
AIC_conexp_bbinom_text$label <- 
  paste(AIC_conexp_bbinom_text$year, ": ", AIC_conexp_bbinom_text$AIC)

total_AIC_conexp_bbinom <- sum(comb_AIC_conexp_bbinom)

# Create the plot
ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[1]]['b']*coefs_conexp_bbinom[[1]]['x100']))*exp(-coefs_conexp_bbinom[[1]]['b']*x), 0.999), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[2]]['b']*coefs_conexp_bbinom[[2]]['x100']))*exp(-coefs_conexp_bbinom[[2]]['b']*x), 0.999), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun =
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[3]]['b']*coefs_conexp_bbinom[[3]]['x100']))*exp(-coefs_conexp_bbinom[[3]]['b']*x), 0.999), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[4]]['b']*coefs_conexp_bbinom[[4]]['x100']))*exp(-coefs_conexp_bbinom[[4]]['b']*x), 0.999), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = 
                  function(x)ifelse((1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x) < 0.999, (1/exp(-coefs_conexp_bbinom[[5]]['b']*coefs_conexp_bbinom[[5]]['x100']))*exp(-coefs_conexp_bbinom[[5]]['b']*x), 0.999), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A70. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a CNE function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, binomial distribution
##

coefs_hill_binom <- list(c_1314 <- c(NA, NA, NA),
                         c_1415 <- c(NA, NA, NA),
                         c_1516 <- c(NA, NA, NA),
                         c_1617 <- c(NA, NA, NA),
                         c_1718 <- c(NA, NA, NA))

comb_AIC_hill_binom <- rep(NA, times = length(seasons))

fun_hill_binom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30)

for(i in 1:length(seasons)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_binom <- optim(par = pars, fn = fun_hill_binom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_binom$par
  coefs_hill_binom[[i]] <- c_opt
  comb_AIC_hill_binom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_binom$value)
}

AIC_hill_binom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_binom, digits = 0), tjust = tjust)
AIC_hill_binom_text$label <- paste(AIC_hill_binom_text$year, ": ", AIC_hill_binom_text$AIC)

total_AIC_hill_binom <- sum(comb_AIC_hill_binom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_binom[[1]]['a']*(1 - (x^(coefs_hill_binom[[1]]['n'])/(coefs_hill_binom[[1]]['h']^(coefs_hill_binom[[1]]['n']) + x^(coefs_hill_binom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[2]]['a']*(1 - (x^(coefs_hill_binom[[2]]['n'])/(coefs_hill_binom[[2]]['h']^(coefs_hill_binom[[2]]['n']) + x^(coefs_hill_binom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[3]]['a']*(1 - (x^(coefs_hill_binom[[3]]['n'])/(coefs_hill_binom[[3]]['h']^(coefs_hill_binom[[3]]['n']) + x^(coefs_hill_binom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[4]]['a']*(1 - (x^(coefs_hill_binom[[4]]['n'])/(coefs_hill_binom[[4]]['h']^(coefs_hill_binom[[4]]['n']) + x^(coefs_hill_binom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_binom[[5]]['a']*(1 - (x^(coefs_hill_binom[[5]]['n'])/(coefs_hill_binom[[5]]['h']^(coefs_hill_binom[[5]]['n']) + x^(coefs_hill_binom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_binom_text, 
            aes(x = 100, y = AIC_hill_binom_text$tjust, label = AIC_hill_binom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_binom)),
                hjust = 0)) +
  labs(subtitle = "Deterministic: Hill function; Stochastic: Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A71. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a Hill function with binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


##
## Hill function, beta-binomial distribution
##

coefs_hill_bbinom <- list(c_1314 <- c(NA, NA, NA, NA),
                          c_1415 <- c(NA, NA, NA, NA),
                          c_1516 <- c(NA, NA, NA, NA),
                          c_1617 <- c(NA, NA, NA, NA),
                          c_1718 <- c(NA, NA, NA, NA))

comb_AIC_hill_bbinom <- rep(NA, times = length(seasons))

fun_hill_bbinom <- function(par, x, z, n){
  mu <- (par['a']*(1 - (x^(par['n'])/(par['h']^(par['n']) + x^(par['n'])))))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par['theta'], log = TRUE))
  return(nll)
}
pars <- c(a = 1, n = 2, h = 30, theta = 1)

for(i in 1:length(seasons)){
  t <- dat_8[[i]]
  x <- t$dist
  z <- t$pos
  n <- t$n
  
  opt_hill_bbinom <- optim(par = pars, fn = fun_hill_bbinom, x = x, z = z, n = n, method = "Nelder-Mead", control = list(maxit = 1000))
  c_opt <- opt_hill_bbinom$par
  coefs_hill_bbinom[[i]] <- c_opt
  comb_AIC_hill_bbinom[i] <- 2*length(pars) + 2*as.numeric(opt_hill_bbinom$value)
}

AIC_hill_bbinom_text <- data.frame(year = seasons, AIC = round(comb_AIC_hill_bbinom, digits = 0), tjust = tjust)
AIC_hill_bbinom_text$label <- paste(AIC_hill_bbinom_text$year, ": ", AIC_hill_bbinom_text$AIC)

total_AIC_hill_bbinom <- sum(comb_AIC_hill_bbinom)

ggplot() +
  geom_point(data = dat_8[[1]], aes(x = dist, prop, color = "2013/2014"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[2]], aes(x = dist, prop, color = "2014/2015"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[3]], aes(x = dist, prop, color = "2015/2016"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[4]], aes(x = dist, prop, color = "2016/2017"), shape = 1, position = "jitter") +
  geom_point(data = dat_8[[5]], aes(x = dist, prop, color = "2017/2018"), shape = 1, position = "jitter") +
  stat_function(fun = function(x)(coefs_hill_bbinom[[1]]['a']*(1 - (x^(coefs_hill_bbinom[[1]]['n'])/(coefs_hill_bbinom[[1]]['h']^(coefs_hill_bbinom[[1]]['n']) + x^(coefs_hill_bbinom[[1]]['n']))))), 
                aes(color = "2013/2014", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[2]]['a']*(1 - (x^(coefs_hill_bbinom[[2]]['n'])/(coefs_hill_bbinom[[2]]['h']^(coefs_hill_bbinom[[2]]['n']) + x^(coefs_hill_bbinom[[2]]['n']))))), 
                aes(color = "2014/2015", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[3]]['a']*(1 - (x^(coefs_hill_bbinom[[3]]['n'])/(coefs_hill_bbinom[[3]]['h']^(coefs_hill_bbinom[[3]]['n']) + x^(coefs_hill_bbinom[[3]]['n']))))), 
                aes(color = "2015/2016", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[4]]['a']*(1 - (x^(coefs_hill_bbinom[[4]]['n'])/(coefs_hill_bbinom[[4]]['h']^(coefs_hill_bbinom[[4]]['n']) + x^(coefs_hill_bbinom[[4]]['n']))))), 
                aes(color = "2016/2017", size = "s")) +
  stat_function(fun = function(x)(coefs_hill_bbinom[[5]]['a']*(1 - (x^(coefs_hill_bbinom[[5]]['n'])/(coefs_hill_bbinom[[5]]['h']^(coefs_hill_bbinom[[5]]['n']) + x^(coefs_hill_bbinom[[5]]['n']))))), 
                aes(color = "2017/2018", size = "s")) +
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:") +
  geom_text(data = AIC_hill_bbinom_text, 
            aes(x = 100, y = AIC_hill_bbinom_text$tjust, label = AIC_hill_bbinom_text$label, hjust = 0)) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_hill_bbinom)),
                hjust = 0)) +
  labs( 
    subtitle = "Deterministic: Hill function; Stochastic: Beta-Binomial distribution;
       Type of data: negative clairvoyance & positive carry-over data", 
    x = "Distance", 
    y = "Proportion positive",
    caption = "Figure A72. A graph with the datapoints of the negative clairvoyance & positive carry-over data 
       and lines as fitted with a Hill function with beta-binomial distribution.") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Season", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(plot.caption = element_text(hjust = 0))