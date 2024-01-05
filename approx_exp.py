from sympy import exp, series, lcm, Rational, Interval, QQ, ZZ, cos, pi, N, expand, nsimplify, Poly
from sympy.polys.polyfuncs import interpolate

def QQ_to_ZZ_poly(p):
    common_denominator = lcm([Rational(x).denominator for x in p.coeffs()])
    return (p*common_denominator).as_poly(domain=ZZ), common_denominator

from sympy import *
def exp_taylor(variable, basepoint, degree):
    """Calculate taylor series, dropping error term and representing as a polynomial over QQ"""
    p = series(exp(variable), x=variable, x0=basepoint, n=degree).removeO().as_poly(x)
    coeffs =[Rational(N(c)) for c in p.coeffs()] # !!!
    return Poly(coeffs, x)


def affine_pm_one(interval, x):
    """Send [-1, 1] to [a, b]"""
    u =  (x + 1) * Rational(1, 2) # on [0, 1]
    return u*interval.measure + interval.start

def chebyshev_nodes(interval, n):
    """
    See: https://en.wikipedia.org/wiki/Chebyshev_nodes
    """
    nodes = []
    for k in range(1, n+1):
        b_k = cos((k - 1)*pi/(n-1))
        x_k = affine_pm_one(interval, b_k)
        nodes.append(x_k)
    return nodes

def exp_taylor_interval(interval, variable, degree):
    """Approximate exp on an interval using taylor around middle. """
    median = Rational(interval.start + interval.end, 2)
    return exp_taylor(variable, median, degree)

def exp_chebyshev_nodes(interval, variable, degree):
    """Approximate exp on an interval using chebyshev nodes. """
    points = []
    for x_n in chebyshev_nodes(interval, degree):
        x_n = Rational(N(x_n)) #!!!
        y_n = Rational(N(exp(x_n))) # !!!
        points.append((x_n, y_n))
    return expand(interpolate(points, variable)).as_poly(domain=QQ)

if __name__ == '__main__':
    from sympy.abc import x
    f1 = exp_taylor_interval(Interval(-1, 0), x, 5)
    f2 = exp_chebyshev_nodes(Interval(-1, 0), x, 5)

    print('taylor = ', f1)
    print('cheby = ', f2)
    print("t, exp(t), taylor(t), cheby(t)")
    for i in range(10):
        t = i/10 - 1
        print(t, exp(t), N(f1.subs(x, t)), N(f2.subs(x, t)))
