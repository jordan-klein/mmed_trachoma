
#S = susceptible 
#Ip1 = Infectious TF/TI negative
#It1 = Infectious TF/TI positve
#D1 = No longer infectious
#Si = Susceptible but previously infected
#Ipi = Infectious TF/TI negative (previously infected)
#Iti = Infectious TF/TI positive (previously infected)
#Di = No long infectious (previously infected)

dS1dt = b - (lambda + mu)*S1
dIp1dt = lambda*S1 - (alpha + mu)*Ip1
dIt1dt = alpha*Ip1 - (gamma1 + mu)*It1
dD1dt = gamma1*It1 - (omega1 + (lambda*sigma) + mu)*D1
dSidt = (omega1 + omegai)*D1 - (lambda + mu)*Si
dIpidt = lambda*Si - (alpha + mu)*Ipi
dItidt = (alpha*Ipi) + ((lambda*sigma)*D1) + ((lambda*sigma)*Di) - (gammai + mu)*Iti 
dDidt = gammai*Iti - (omegai + mu + (lambda*sigma)*Di)

lamba = (beta*(Ip1+It1+Ipi+Iti)/N)

##Parameters
  #b = births constant [per day]
  #mu = deaths constant [per day]
  #lambda = rate of transmission [1/days]
  #alpha = rate at which people move develop clinical signs [1/days]
  #gamma1 = rate at which people become non-infectious/clear infection during initial infection [1/days]
  #sigma = rate of infection [person per person per day]
  #omega1 = recovery rate for initial infection [1/days]
  #gammai = rate at which people become non-infectious/clear infection (among reinfected) [1/days]
  #beta = force of infection 





## function for infectivity with repeated infections

r_IA =  rep (0, 100)	## rate of leaving IA (weeks)- functional form

nun = 1/((7.7*7)/log(2)) ####2.71*30## pars that go into the function for recovery

nu0 = 0.0033 #####1.09*30 ## pars that go into the function for recovery

inf_red = 0.4509

for(j in 1: 100)
{ 
  r_IA[j] <-  (nu0-nun)*exp(-inf_red*(j-1))+ nun ### 1/(8.3*7) ##	 ## 	###
}
r_IA
####################################################
## function for infectivity with repeated infections

r_A =  rep (0, 100)	## rate of leaving IA (weeks)- functional form

nun_A = 1/((1*7)/log(2)) ####2.71*30## pars that go into the function for recovery

nu0_A = 0.0050 ## 1/((20.7*7)/log(2)) #####1.09*30 ## pars that go into the function for recovery

inf_red_A = 0.3065

for(k in 1: max(N_infs))
{ 
  r_A[k] <-  (nu0_A - nun_A )*exp(-inf_red_A*(k-1))+ nun_A ###  1/((20.7*7)/log(2)) ## 1/(7.14*7) ##	1/(36.1*7) ###
}