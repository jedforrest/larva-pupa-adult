#!/usr/bin/env pythona
# lpa_simulation.py -- simulate system described in Constantino et al.  

from math import exp
import larva_pupa_adult_params as P

def step(l, p, a):
    # Eqn (1)
    l_new = P.b * a * exp(-( P.c_el*l + P.c_ea*a))
    # Eqn (2)
    p_new = l*(1 - P.mu_l)
    # Eqn (3)
    a_new = a*(1 - P.mu_a) + p* exp(-( P.c_pa*a ))
    return l_new, p_new, a_new

# Simulation
number_of_generations = 10
larva_pop = P.initial_larva_population
pupa_pop = P.initial_pupa_population
adult_pop = P.initial_adult_population
for gen in range(number_of_generations):
    print("%d,%s,%s,%s" % (gen, str(larva_pop), str(pupa_pop), str(adult_pop)))
    larva_pop, pupa_pop, adult_pop = step(larva_pop, pupa_pop, adult_pop)
