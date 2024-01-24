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
        poly = expand(a*(1-mu_a) + p * exp_approx(-c_pa*a) - a_next)
        polynomials.append(poly)
        debug("basepoint approx = ", exp_input, color=32)
        debug("poly approx error = ", poly.subs(mu_a, lpa_params.mu_a).subs(c_pa, lpa_params.c_pa), color=33)
    return polynomials


### The default values of off make sense for the default values in lpa_params. We should update this.
def gensys(prolog, degree, offsum=[0,0,0], offmul = [1,1,1], printout=False):
        L, P, A = simulate (prolog+1)

        truesolLarva = "\ntruesolL = [{} ; {} ; {}];".format(lpa_params.c_el,lpa_params.c_ea,lpa_params.b)
        truesolAdult = "\ntruesolA = [{} ; {}];".format(lpa_params.c_pa,lpa_params.mu_a)

        lpoly = larva_poly_system_withapprox(degree, A, L, offmul[0]*lpa_params.c_el+offsum[0], offmul[1]*lpa_params.c_ea+offsum[1])
        L_sys = "fL = ["
        for poly in lpoly:
            L_sys +=  str(poly).replace("**","^") + "; "

        L_sys = L_sys[:-3]
        L_sys += "];"
 
        apoly = adult_poly_system_withapprox(degree, A, P, offmul[2]*lpa_params.c_pa+offsum[2])
        A_sys = "fA = ["
        for poly in apoly:
            A_sys +=  str(poly).replace("**","^") + "; "

        A_sys = A_sys[:-3]
        A_sys += "];"

        if printout:
            print("howoffsum = {};".format(offsum))
            print("howoffmul = {};".format(offmul))
            print(truesolLarva)
            print(L_sys)
            print(truesolAdult)
            print(A_sys)

        return L_sys, A_sys

def all_poly_system_withapprox(degree, L, P, A, offmul = [1,1,1]):

    lpoly = larva_poly_system_withapprox(degree, A, L, offmul[0]*lpa_params.c_el, offmul[1]*lpa_params.c_ea)
    L_sys = "fL = ["
    for poly in lpoly:
        L_sys +=  str(poly).replace("**","^") + "; "
    L_sys = L_sys[:-3]
    L_sys += "];"

    apoly = adult_poly_system_withapprox(degree, A, P, offmul[2]*lpa_params.c_pa)
    A_sys = "fA = ["
    for poly in apoly:
        A_sys +=  str(poly).replace("**","^") + "; "
    A_sys = A_sys[:-3]
    A_sys += "];"

    return L_sys, A_sys
