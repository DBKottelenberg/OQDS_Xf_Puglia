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
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:", size = 6) +
  geom_text(data = AIC_nexp_binom_text, 
            aes(x = 100, 
                y = AIC_nexp_binom_text$tjust, 
                label = AIC_nexp_binom_text$label, hjust = 0), size = 6) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_binom)),
                hjust = 0),
            size = 6) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Binomial distribution", 
       x = "Distance", 
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:", size = 6) +
  geom_text(data = AIC_nexp_bbinom_text, 
            aes(x = 100, y = AIC_nexp_bbinom_text$tjust, label = AIC_nexp_bbinom_text$label, hjust = 0), size = 6) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_nexp_bbinom)),
                hjust = 0),
            size = 6) +
  labs(subtitle = "Deterministic: Negative exponential; Stochastic: Beta-Binomial distribution", 
       x = "Distance", 
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:", size = 6) +
  geom_text(data = AIC_log_binom_text, 
            aes(x = 100, y = AIC_log_binom_text$tjust, label = AIC_log_binom_text$label, hjust = 0), size = 6) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_binom)),
                hjust = 0),
            size = 6) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Binomial distribution",
       x = "Distance", 
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:", size = 6) +
  geom_text(data = AIC_log_bbinom_text, 
            aes(x = 100, y = AIC_log_bbinom_text$tjust, label = AIC_log_bbinom_text$label, hjust = 0), size = 6) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_log_bbinom)),
                hjust = 0),
            size = 6) +
  labs(subtitle = "Deterministic: Logistic function; Stochastic: Beta-Binomial distribution", 
       x = "Distance", 
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0,  (max_pos + 60))) +
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:", size = 6) +
  geom_text(data = AIC_conexp_binom_text, 
            aes(x = 100, 
                y = AIC_conexp_binom_text$tjust, 
                label = AIC_conexp_binom_text$label, hjust = 0),
            size = 6) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_binom)),
                hjust = 0),
            size = 6) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Binomial distribution", 
       x = "Distance", 
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))


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
  annotate(geom = "text", x = 112.5, y = 0.85, label = "AIC:", size = 6) +
  geom_text(data = AIC_conexp_bbinom_text, 
            aes(x = 100, 
                y = AIC_conexp_bbinom_text$tjust, 
                label = AIC_conexp_bbinom_text$label, hjust = 0),
            size = 6) +
  geom_text(aes(x = 100,
                y = 0.40,
                label = paste("Total AIC: ", round(total_AIC_conexp_bbinom)),
                hjust = 0),
            size = 6) +
  labs(subtitle = "Deterministic: CNE; Stochastic: Beta-Binomial distribution", 
       x = "Distance", 
       y = "Proportion positive") +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_x_continuous(limits = c(0, (max_pos + 60))) +   # Set the maximum x-axis value lower than what is measured to increase readability
  scale_color_manual(name = "Year", values = c("black", "red", "green3", "blue", "orange", "magenta3")) +
  scale_size_manual(values = 1.00) +
  guides(size = FALSE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))
