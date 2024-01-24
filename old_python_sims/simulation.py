#!/usr/bin/python3
"""
simulation.py -- Simulate larva pupa adult ssytem described in Constantino et al.  
"""
from math import exp
import lpa_params


def step(l, p, a, **params):
    b = params.get("b", lpa_params.b)
    c_el = params.get("c_el", lpa_params.c_el)
    c_ea = params.get("c_ea", lpa_params.c_ea)
    c_pa = params.get("c_pa", lpa_params.c_pa)
    mu_l = params.get("mu_l", lpa_params.mu_l)
    mu_a = params.get("mu_a", lpa_params.mu_a)

    l_new = b * a * exp(-(c_el * l + c_ea * a))
    p_new = l * (1 - mu_l)
    a_new = a * (1 - mu_a) + p * exp(-(c_pa * a))

    return (l_new, p_new, a_new)


def simulate(n, initial_populations=None, **kwargs):
    """
    Simulate system with n prolongations and an intial state. Use keyword arguments to set parameter values.

    Retrurns 
    """
    if initial_populations is None:
        initial_populations = lpa_params.initial_populations
    dataset = [list(initial_populations)]
    while len(dataset) < n:
        dataset.append(step(*dataset[-1], **kwargs))
    L, P, A = list(zip(*dataset))
    return list(L), list(P), list(A)


if __name__ == "__main__":
    import sys
    n = 4
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    dataset = simulate(n)
    for i, (l, p, a) in enumerate(zip(*dataset)):
        print(f'{l},{p},{a}')
