#!/usr/bin/env pythona
# lpa_simulation.py -- simulate system described in Constantino et al.  

from math import exp
import larva_pupa_adult_params as P

def step(l, p, a, round_population_count=True):
    # Eqn (1)
    l_new = P.b * a * exp(-( P.c_el*l + P.c_ea*a))
    # Eqn (2)
    p_new = l*(1 - P.mu_l)
    # Eqn (3)
    a_new = a*(1 - P.mu_a) + p* exp(-( P.c_pa*a ))

    # Do we assume population must be an integer?
    # If round_population_count=True On the 100th gen: 99 24.000000 19.000000 376.000000
    # If round_population_count=False On the 100th gen: 99 22.765323 18.092817 384.228320
    if round_population_count:
        l_new = round(l_new)
        p_new = round(p_new)
        a_new = round(a_new)

    return l_new, p_new, a_new

with open('lpa-data.txt', 'w') as lpa_data_file:
    # Simulation
    number_of_generations = 100 
    larva_population = P.initial_larva_population
    pupa_population = P.initial_pupa_population
    adult_population = P.initial_adult_population
    for gen in range(number_of_generations):
        print("%d %f %f %f" % (gen, larva_population, pupa_population, adult_population), file=lpa_data_file)
        larva_population, pupa_population, adult_population = step(larva_population, pupa_population, adult_population)
