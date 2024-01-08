'''
approximate_exponential.py

We approximate the exponential function
'''
import sympy
from fractions import Fraction
from sympy import exp, symbols, Interval, Rational
from sympy import exp, series, lcm, Rational, Interval, QQ, ZZ, cos, pi, N, expand, nsimplify, Poly
#from sympy.polys.polyfuncs import interpolate

DEGREE = 8
def rational_approx(p, x):
    coeffs =[Rational(N(c)) for c in p.as_poly(x).coeffs()]
    return Poly(coeffs, x)

def exp_taylor(variable, basepoint, degree):
    """Calculate taylor series, dropping error term and representing as a polynomial over QQ"""
    p = series(exp(variable), x=variable, x0=basepoint, n=degree).removeO().as_poly(variable)
    print("exp taylor", p)
    return p


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
    interval_center = N(interval.start + interval.end) / 2
    print("> Taylor interval", interval, interval_center)
    return exp_taylor(variable, interval_center, degree)

def exp_chebyshev_nodes(interval, variable, degree):
    """Approximate exp on an interval using chebyshev nodes. """
    points = []
    for x_n in chebyshev_nodes(interval, degree):
        y_n = exp(x_n)
        points.append((x_n, y_n))
    p = expand(interpolate(points, variable)).as_poly(domain=QQ)
    return rational_approx(p, variable)


def poly_to_julia(p):
    return str(expand(p.as_expr())).replace('**', '^')

data = []
with open('lpa_simulation.output.csv') as csv:
    for line in csv:
        i, l, p, a = line.split(',')
        i = int(i)
        l = float(l)
        p = float(p)
        a = float(a)
        data.append((i, l, p, a))


interval = Interval(-10, 0)

import larva_pupa_adult_params as PARAMS

PARAM_symbols = symbols(f'b c_el c_ea c_pa mu_a mu_l')
b, c_el, c_ea, c_pa, mu_a, mu_l = PARAM_symbols
param_substitution = [
            (b, PARAMS.b),
            (c_el, PARAMS.c_el),
            (c_ea, PARAMS.c_ea),
            (c_pa, PARAMS.c_pa),
            (mu_a, PARAMS.mu_a),
#            (mu_l, PARAMS.mu_l),
]
cannibalism_of_eggs_by_larva_rate = 1.209e-2
cannibalism_of_eggs_by_adult_rate = 1.155e-2
cannibalism_of_pupa_by_adult_rate = 4.700e-3

percentage = 1.10

c_el_interval = Interval((1-percentage)*PARAMS.c_el, (1+percentage)*PARAMS.c_el)
c_ea_interval = Interval((1-percentage)*PARAMS.c_ea, (1+percentage)*PARAMS.c_ea)
c_pa_interval = Interval((1-percentage)*PARAMS.c_pa, (1+percentage)*PARAMS.c_pa)

system_symbols = set(s for s in PARAM_symbols if s != mu_l)
l_equations = []
a_equations = []
equations = []
for (i, l_float, p_float, a_float) in data:
    print(f'generation = {i}')
    print(i, l, p, a)

    # We define symbols for python
    l_new, l, p_new, p, a_new, a = symbols(f'l_{i+1} l_{i} p_{i+1} p_{i} a_{i+1} a_{i}')

    # larva formula
    delta_l = (b * a * exp(-( c_el*l + c_ea*a))) - l_new
    
    exp1_input_minimum = -(c_el_interval.end*l_float + c_ea_interval.end*a_float)
    exp1_input_maximum = -(c_el_interval.start*l_float + c_ea_interval.start*a_float)
    exp1_interval = Interval(exp1_input_minimum, exp1_input_maximum)
   
    from sympy.abc import t
    exp1_poly = exp_taylor_interval(exp1_interval, t, DEGREE)
    print("exp1 interval", exp1_interval)
    print("exp1 poly", exp1_poly)
    delta_l_taylor = (b * a * exp1_poly.subs(t, -( c_el*l + c_ea*a))) - l_new

    # adult formula
    delta_a = a*(1 - mu_a) + p * exp(-( c_pa*a )) - a_new

    exp2_input_minimum = -(c_pa_interval.end*a_float)
    exp2_input_maximum = -(c_pa_interval.start*a_float)
    exp2_interval = Interval(exp2_input_minimum, exp2_input_maximum)
    exp2_poly = exp_taylor_interval(exp2_interval, t, DEGREE)
    print("exp2 interval", exp2_interval)
    print("exp2 poly", exp2_poly)
    delta_a_taylor = a*(1 - mu_a) + p * exp2_poly.subs(t, -( c_pa*a )) - a_new

    if i < 3:
        equations.append(str(expand(delta_l_taylor.as_expr())).replace('**', '^'))
    if i < 2:
        equations.append(str(expand(delta_a_taylor.as_expr())).replace('**', '^'))
    if i >= 3:
        break

    # store them to input them in Julia
    system_symbols.update([l_new, l, p_new, p, a_new, a])

    l_new_float = data[i+1][1]
    p_new_float = data[i+1][2]
    a_new_float = data[i+1][3]

    # replace the 6 experimental data variables with float
    data_substitution = [
        (l_new, l_new_float),
        (l, l_float),
        (p_new, p_new_float),
        (p, p_float),
        (a_new, a_new_float),
        (a, a_float),
    ]

    for variable, variable_value in data_substitution:
        if variable in [l, p, a] or i == 2:
            equations.append(variable - variable_value)

    substitution = data_substitution + param_substitution
    delta_float = delta_l.subs(substitution), delta_a.subs(substitution)
    print("float error", delta_float)
    my_delta_float = delta_l_taylor.subs(substitution), delta_a_taylor.subs(substitution)
    print("taylor error", my_delta_float)
    print("input to exp function", (-( c_el*l + c_ea*a)).subs(substitution), (-( c_pa*a )).subs(substitution))
    print()
    if i == 3:
        break

filename = "hc_input.jl"
htop_jl_input = open(filename, 'w')

print('using HomotopyContinuation', file=htop_jl_input)
system_symbols = sorted(str(x) for x in system_symbols)
print('@var %s' % (' '.join(system_symbols)), file=htop_jl_input)
print('f = System([', file=htop_jl_input)
for equation in equations:
    print("    %s," % equation, file=htop_jl_input)
print('])', file=htop_jl_input)
print('result = solve(f)', file=htop_jl_input)
print('''
open("julia_output.txt", "w") do io
	println(io, real_solutions(result))
end
''', file=htop_jl_input)
