# LPA model with taylor approximated exp function
using AbstractAlgebra
using LinearAlgebra
using TaylorSeries

S, params = polynomial_ring(QQ, ["b", "cel", "cea", "cpa", "μl", "μa"])
R, u = polynomial_ring(S, ["L" , "P", "A"])

l, p, a, = u
b, cel, cea, cpa, μl, μa = params

function exp_taylor_polyn(R::Ring, series_order::Int)
    # TaylorSeries is not compatible with PolynomialRings
    # so convert series approximation to polynomial with same coefficients
    exp_taylor = convert(Taylor1{Rational{Int}}, taylor_expand(exp, 0; order=series_order))
    coeffs = getcoeff.(exp_taylor, 0:series_order)
    return polynomial(R, coeffs)
end

exp_taylor = exp_taylor_polyn(R, 2)
exp_taylor(4)
exp_taylor(2a)
exp_taylor(a+l)
exp_taylor(b*p*a)
exp_taylor(l*p*a)


du = zeros(R, length(u))
LPA_taylor(du, u, params; exp_taylor=exp_taylor_polyn(R, 2))
du

du[1](1,2,3)

# TODO
# - remove excess test code
# - ring elements l0, l1, l2 ... a3
# - connect with data from simulation
# - AA rref
# - solve!

function LPA_taylor(du, u, p; exp_taylor)
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    # exp_taylor = taylor_expand(exp, 0; order)

    du[1] = b * A * exp_taylor(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp_taylor(-cpa * A) + (1 - μa) * A
    return nothing
end

function prolongate_LPA(LPA_model, vars, params; nsteps=1, order=1)
    du = zeros(Expression, 3)
    u = collect(zip(L, P, A))
    prolongations = Expression[]

    for i in 1:nsteps
        du = LPA_model(du, u[i], params, i; order) .- u[i+1]
        append!(prolongations, du)
    end
    return prolongations
end

# Simulation script
nsteps = 5
vars = @var L[0:nsteps] P[0:nsteps] A[0:nsteps]
params = @var b cel cea cpa μl μa

# generated in other file
include("LPA_simulations.jl")
data = round.(Int, LPAdata[:,1:nsteps+1])


eqns = prolongate_LPA(LPA_taylor, vars, params; nsteps, order=5)

eqns_eval = HomotopyContinuation.evaluate(eqns, data_map)
eqns_expand = expand.(eqns_eval)

lpa = [L; P; A]

eqns_expand[1:6]

F = System(eqns_expand[1:6]; variables=params, parameters=lpa)
