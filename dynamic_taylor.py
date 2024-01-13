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

def larva_poly_system(degree, adult_pop_series, larva_pop_series, b_interval, c_el_interval, c_ea_interval):
    b, c_el, c_ea = symbols(f"b c_el c_ea")
    b_min, b_max = b_interval
    c_el_min, c_el_max = c_el_interval
    c_ea_min, c_ea_max = c_ea_interval
    assert len(adult_pop_series) == len(larva_pop_series), "Adult and larva population time series should have same length."
    polynomials = []
    for i, (a, l) in enumerate(zip(adult_pop_series, larva_pop_series)):
        if i == len(adult_pop_series)-1:
            break
        exp_input_min = -(c_el_max*l + c_ea_max*a)
        exp_input_max = -(c_el_min*l + c_ea_min*a)
        exp_input_middle = (exp_input_max + exp_input_min) / 2
        exp_approx = taylor_approx(degree, exp_input_middle)
        l_next = larva_pop_series[i+1]
        poly = expand(b * a * exp_approx(-(c_el * l + c_ea * a)) - l_next)
        polynomials.append(poly)
        debug("basepoint range = ", exp_input_min, exp_input_max, color=32)
        debug("poly approx error = ", poly.subs(b, lpa_params.b).subs(c_el, lpa_params.c_el).subs(c_ea, lpa_params.c_ea), color=33)
    return polynomials

### def pupa_poly_system(degree, larva_pop_series, pupa_pop_system, mu_l_interval):
###     mu_l = symbols(f"mu_l")
###     assert len(larva_pop_series) == len(pupa_pop_system)
###     for i, (l, p) in enumerate(zip(larva_pop_series, pupa_pop_system)):
###         if i == len(larva_pop_series)-1
###     # TODO
###     pass


def adult_poly_system(degree, adult_pop_series, pupa_pop_series, c_pa_interval, mu_a_interval):
    c_pa, mu_a = symbols(f"c_pa mu_a")
    c_pa_min, c_pa_max = c_pa_interval
    #mu_a_min, mu_a_max = mu_a_interval
    assert len(adult_pop_series) == len(pupa_pop_series), "Adult and pupa population time series should have same length."
    polynomials = []
    for i, (a, p) in enumerate(zip(adult_pop_series, pupa_pop_series)):
        if i == len(adult_pop_series)-1:
            break
        # a_new = a * (1 - mu_a) + p * exp(-(c_pa * a))
        exp_input_min = -(c_pa_max*p)
        exp_input_max = -(c_pa_min*p)
        exp_input_middle = (exp_input_max + exp_input_min) / 2
        exp_approx = taylor_approx(degree, exp_input_middle)
        a_next = adult_pop_series[i+1]
        poly = expand(a*(1-mu_a) + p * exp_approx(-c_pa*p) - a_next)
        polynomials.append(poly)
        debug("basepoint range = ", exp_input_min, exp_input_max, color=32)
        debug("poly approx error = ", poly.subs(mu_a, lpa_params.mu_a).subs(c_pa, lpa_params.c_pa), color=33)
    return polynomials




if __name__ == "__main__":
    from random import random
    for degree in range(1, 10):
        print(('degree', degree))
        default_interval_min = lambda: 0.757 + random()/100
        default_interval_max = lambda: 1.22 + random()/100
        default_intervals = {
                "b":    [default_interval_min()*lpa_params.b,    default_interval_max()*lpa_params.b],
                "c_el": [default_interval_min()*lpa_params.c_el, default_interval_max()*lpa_params.c_el],
                "c_ea": [default_interval_min()*lpa_params.c_ea, default_interval_max()*lpa_params.c_ea],
                "c_pa": [default_interval_min()*lpa_params.c_pa, default_interval_max()*lpa_params.c_pa],
                "mu_a": [default_interval_min()*lpa_params.mu_a, default_interval_max()*lpa_params.mu_a],
        }
        L, P, A = simulate(4)
        print(repr("larva poly system"))
        lpoly = larva_poly_system(degree, A, L,
                b_interval=default_intervals["b"],
                c_el_interval=default_intervals["c_el"],
                c_ea_interval=default_intervals["c_ea"])
        polynomials = []
        for poly in lpoly:
            polynomials.append(pretty_poly(poly))
        print(repr(("True-solution", {"b": lpa_params.b, "c_el": lpa_params.c_el, "c_ea": lpa_params.c_ea})))
        start_time = time()
        for i, solution in enumerate(solve(["b", "c_el", "c_ea"], polynomials)):
            print(repr((i, solution)))
        print(("total-time", time() - start_time))

        apoly = adult_poly_system(degree, A, P,
                c_pa_interval=default_intervals["c_pa"],
                mu_a_interval=default_intervals["mu_a"])
        print("adult poly system")
        polynomials = []
        for poly in apoly[:2]:
            polynomials.append(pretty_poly(poly))
        print("True solution", {"c_pa": lpa_params.c_pa, "mu_a": lpa_params.mu_a})
        start_time = time()
        for i, solution in enumerate(solve(["c_pa", "mu_a"], polynomials)):
            print(repr((i, solution)))
        print(("total-time", time() - start_time))
