#!/usr/bin/env python3

from collections import namedtuple
import numpy as np


LPAParameters = namedtuple("LPAParameters", "b c_el c_ea mu_l c_pa mu_a")
# These parameters come from footnote 15 in the 1997 article
Params = LPAParameters(
        b=6.598,
        c_el=1.209e-2,
        c_ea=1.155e-2,
        mu_l=0.2055,
        c_pa=4.700e-3,
        mu_a=7.629e-3)

def make_variance(n=1):
    mean = [0, 0, 0]
    s11=0.3412
    s22=0.2488
    s33=1.627e-4
    s12=7.312e-2
    s13=-1.719e-3
    s23=3.374e-4
    cov = [
            [s11, s12, s13],
            [s12, s22, s23],
            [s13, s23, s33]
    ]
    return np.random.multivariate_normal(mean, cov, n)


def step(state):
    L, P, A = state
    L_ = Params.b*A*np.exp(-Params.c_el*L - Params.c_ea*A)
    P_ = L*(1 - Params.mu_l)
    A_ = A*(1 - Params.mu_a) + P*np.exp(-Params.c_pa*A)
#    return np.round(L_), np.round(P_), np.round(A_)
    return L_, P_, A_

N = 10000
L, P, A = np.zeros(N), np.zeros(N), np.arange(0, N, 1)
for g in range(10000):
    print("gen =", g)
    print("L =", L)
    print("P =", P)
    print("A =", A)
    print()
    L, P, A = step((L, P, A))
