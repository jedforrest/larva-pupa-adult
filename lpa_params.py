# This file contains the maximum likelihood parameter estimates, covariance estates, etc. from
#
#    ``Chaotic Dynamics in an Insect Population''
#    by R. F. Costantino, R. A. Desharnais,* J. M. Cushing, and Brian Dennis
#    Science 275.5298 (1997): 389-391.
#
# All rates are per two weeks.
birth_rate = 6.598
cannibalism_of_eggs_by_larva_rate = 1.209e-2
cannibalism_of_eggs_by_adult_rate = 1.155e-2
larva_mortality_rate = 0.2055
cannibalism_of_pupa_by_adult_rate = 4.700e-3
adult_mortality_rate = 7.629e-3
 
# Aliases which make my variable names roughly agree with notation in "Chaotic Dynamics..."
b = birth_rate
c_el = cannibalism_of_eggs_by_larva_rate
c_ea = cannibalism_of_eggs_by_adult_rate
mu_l = larva_mortality_rate
c_pa = cannibalism_of_pupa_by_adult_rate
mu_a = adult_mortality_rate

# See paragraph ``We experimentally set the adult mortality rate...''
initial_larva_population = 250
initial_pupa_population = 5
initial_adult_population = 100
initial_populations = (initial_larva_population, initial_pupa_population, initial_adult_population)

# TODO: Add constants from footnote on covariance estimates of control/maninpulated
