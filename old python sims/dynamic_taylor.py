"""
dynamic_taylor.py

We make LPA a polynomial system by approximating the exponential function using taylor at a fixed basepoint
"""
from math import exp
from time import time
from sys import argv

from sympy import expand, symbols

from homotopy_continuation import solve
import lpa_params
from simulation import simulate

from random import random

DEBUG=False

def factorial(n):
    if n == 0:
        return 1
    return n * factorial(n-1)

def taylor_approx(degree, basepoint):
    """
    >>> e = taylor_approx(10, 0)
    >>> round(e(1), 5)
    2.71828
    """
    taylor_input = symbols("taylor_input")
    taylor_poly = expand(exp(basepoint)*sum((taylor_input-basepoint)**k/factorial(k) for k in range(degree+1)))
    def exp_taylor(x):
        return taylor_poly.subs(taylor_input, x)
    return exp_taylor

def pretty_poly(p):
    """Make sympy poly better for poly input in other CAS"""
    return str(expand(p).as_expr()).replace("**", "^")

def debug(*args, color=33):
    if DEBUG:
        print(f"\033[{color};1m", *args, "\033[0m")

def larva_poly_system_withapprox(degree, adult_pop_series, larva_pop_series, c_el_approx, c_ea_approx):
    b, c_el, c_ea = symbols(f"b c_el c_ea")
    assert len(adult_pop_series) == len(larva_pop_series), "Adult and larva population time series should have same length."
    polynomials = []
    for i, (a, l) in enumerate(zip(adult_pop_series, larva_pop_series)):
        if i == len(adult_pop_series)-1:
            break
        exp_input = - (c_el_approx*l + c_ea_approx*a)
        exp_approx = taylor_approx(degree, exp_input)
        l_next = larva_pop_series[i+1]
        poly = expand(b * a * exp_approx(-(c_el * l + c_ea * a)) - l_next)
        polynomials.append(poly)
        debug("basepoint approx = ", exp_input, color=32)
        debug("poly approx error = ", poly.subs(b, lpa_params.b).subs(c_el, lpa_params.c_el).subs(c_ea, lpa_params.c_ea), color=33)
    return polynomials


def adult_poly_system_withapprox(degree, adult_pop_series, pupa_pop_series, c_pa_approx):
    c_pa, mu_a = symbols(f"c_pa mu_a")
    #mu_a_min, mu_a_max = mu_a_interval
    assert len(adult_pop_series) == len(pupa_pop_series), "Adult and pupa population time series should have same length."
    polynomials = []
    for i, (a, p) in enumerate(zip(adult_pop_series, pupa_pop_series)):
        if i == len(adult_pop_series)-1:
            break
        # a_new = a * (1 - mu_a) + p * exp(-(c_pa * a))
        exp_input = -(c_pa_approx*p)
        exp_approx = taylor_approx(degree, exp_input)
        a_next = adult_pop_series[i+1]
        poly = expand(a*(1-mu_a) + p * exp_approx(-c_pa*p) - a_next)
        polynomials.append(poly)
        debug("basepoint approx = ", exp_input, color=32)
        debug("poly approx error = ", poly.subs(mu_a, lpa_params.mu_a).subs(c_pa, lpa_params.c_pa), color=33)
    return polynomials


### The default values of off make sense for the default values in lpa_params. We should update this.
def gensys(prolog, degree, off=[0.0001,0.0001,0.0001]):
        L, P, A = simulate (prolog+1)
        lpoly = larva_poly_system_withapprox(degree, A, L, lpa_params.c_el+off[0], lpa_params.c_ea+off[1])


        print("howoff = {};".format(off))
        
        truesolLarva = "solL = [{} ; {} ; {}];".format(lpa_params.c_el,lpa_params.c_ea,lpa_params.b)

        print(truesolLarva)
        
        mysys = "fL = ["
        for poly in lpoly:
            strpoly = str(poly)
            strpoly = strpoly.replace("**","^")
            strpoly = strpoly.replace("c_el","c[1]")
            strpoly = strpoly.replace("c_ea","c[2]")
            mysys +=  strpoly.replace("b","c[3]")
            mysys += "; "

        mysys = mysys[:-3]
        mysys += "];"
        print(mysys)
            
        
        apoly = adult_poly_system_withapprox(degree, A, P, lpa_params.c_pa+off[2])

        truesolAdult = "solA = [{} ; {}];".format(lpa_params.c_pa,lpa_params.mu_a)

        print(truesolAdult)
        
        mysys = "fA = ["
        for poly in apoly:
            strpoly = str(poly)
            strpoly = strpoly.replace("**","^")
            strpoly = strpoly.replace("c_pa","d[1]")
            mysys +=  strpoly.replace("mu_a","d[2]")
            mysys += "; "

        mysys = mysys[:-3]
        mysys += "];"
        print(mysys)
