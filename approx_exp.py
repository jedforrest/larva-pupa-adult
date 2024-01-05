from sympy import exp, series, lcm, Rational, Interval, QQ, ZZ, cos, pi, N, expand, nsimplify
from sympy.polys.polyfuncs import interpolate

def QQ_to_ZZ_poly(p):
    common_denominator = lcm([Rational(x).denominator for x in p.coeffs()])
    return (p*common_denominator).as_poly(domain=ZZ), common_denominator

def exp_taylor(variable, basepoint, degree):
    """
    Calculate taylor series, dropping error term and representing as a polynomial over QQ
    """
    return series(exp(variable), x=variable, x0=basepoint, n=degree).removeO().as_poly()




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
        nodes.append(Rational(N(x_k))) # !!!
    return nodes

def exp_taylor_interval(interval, variable, degree):
    """Approximate exp on an interval using taylor around middle. """
    median = Rational(interval.start + interval.end, 2)
    return exp_taylor(variable, median, degree)

def exp_chebyshev_nodes(interval, variable, degree):
    """Approximate exp on an interval using chebyshev nodes. """
    points = []
    for x_n in chebyshev_nodes(interval, degree):
        y_n = Rational(N(exp(x_n))) # !!!
        points.append((x_n, y_n))
        print(x_n, y_n)
    return expand(interpolate(points, variable)).as_poly(domain=QQ)

if __name__ == '__main__':
    from sympy.abc import x
    print(exp_chebyshev_nodes(Interval(-3, 2), x, 5))
    print(exp_taylor_interval(Interval(-3, 2), x, 5))
