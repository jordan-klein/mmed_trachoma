### Trachoma Project ####
######################################################################
## Section 1: Creating a function to describe the behavior of an ODE
######################################################################

rm(list=ls())                   # Clear all variables and functions
library(deSolve)                # Load library to be used for numerical integration
library(tidyverse)
library(data.table)
library(stringr)

#S1 = susceptible 
#Ip1 = Infectious TF/TI negative
#It1 = Infectious TF/TI positve
#D1 = No longer infectious
#Si = Susceptible but previously infected
#Ipi = Infectious TF/TI negative (previously infected)
#Iti = Infectious TF/TI positive (previously infected)
#Di = No long infectious (previously infected)

##-- Model 1: SIID model--##
##
## siid() 
siid <- function(t, y, parms){
  with(c(as.list(y), parms),{
    
    N <- S1 + Ip1 + It1 + D1
    lambda <- beta*(Ip1 + It1 + Ipi + Iti)/N
    
    dS1dt <- b - (lambda + mu)*S1
    dIp1dt <- lambda*S1 - (alpha + mu)*Ip1
    dIt1dt <- alpha*Ip1 - (gamma1 + mu)*It1
    dD1dt <- gamma1*It1 - (omega1 + (lambda*sigma) + mu)*D1
    
    dSidt <- (omega1 + omegai)*D1 - (lambda + mu)*Si
    dIpidt <- lambda*Si - (alpha + mu)*Ipi
    dItidt <- (alpha*Ipi) + ((lambda*sigma)*D1) + ((lambda*sigma)*Di) - (gammai + mu)*Iti 
    dDidt <- gammai*Iti - (omegai + mu + (lambda*sigma)*Di)
    
    list(c(dS1dt, dIp1dt, dIt1dt, dD1dt, 
           dSidt, dIpidt, dItidt, dDidt))
  })
}



##Parameters
## function for infectivity with repeated infections

omega =  rep (0, 100)	## rate of leaving IA (weeks)- functional form

nun = 1/((7.7*7)/log(2)) ####2.71*30## pars that go into the function for recovery

nu0 = 0.0033 #####1.09*30 ## pars that go into the function for recovery

inf_red = 0.4509

for(j in 1: 100)
{ 
  omega[j] <-  (nu0-nun)*exp(-inf_red*(j-1))+ nun ### 1/(8.3*7) ##	 ## 	###
}

####################################################
## function for infectivity with repeated infections

gamma =  rep (0, 100)	## rate of leaving IA (weeks)- functional form

nun_A = 1/((1*7)/log(2)) ####2.71*30## pars that go into the function for recovery

nu0_A = 0.0050 ## 1/((20.7*7)/log(2)) #####1.09*30 ## pars that go into the function for recovery

inf_red_A = 0.3065

for(k in 1: 100)
{ 
  gamma[k] <-  (nu0_A - nun_A )*exp(-inf_red_A*(k-1))+ nun_A ###  1/((20.7*7)/log(2)) ## 1/(7.14*7) ##	1/(36.1*7) ###
}

#### Parameter values
values <- c(b = 37.4, 
            mu = .00007, 
            beta = .076, 
            sigma = .5, 
            alpha = 1/14, 
            omega1 = omega[1], 
            gamma1 = gamma[1], 
            omegai = mean(omega), 
            gammai = mean(gamma))

#### Time
time.out <- seq(0, 365*10, 1)

#### Population
N0 <- 100000 

pop.SIID <- c(S1 = N0/4,
              Ip1 = N0/8,
              It1 = N0/8,
              D1 = 0, 
              Si = N0/4, 
              Ipi = N0/8, 
              Iti = N0/8, 
              Di = 0)      


## Now let's see what happens if we plug our inputs into lsoda()...
ts <- data.table(lsoda(
  y = pop.SIID,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = siid,                   # Function to evaluate
  parms = values                # Vector of parameters
))

ts.long <- melt(ts, id.vars = 'time') %>% 
  mutate(compartment = str_sub(variable, start = 1, end = -2), 
         repeat_infection = case_when(grepl("1", variable) ~ "Naive", 
                                      grepl("i", variable) ~ "Repeat"))


(ggplot(ts.long)
  + aes(x = time, y = value, color = compartment, linetype = repeat_infection)
  + geom_line()
)

