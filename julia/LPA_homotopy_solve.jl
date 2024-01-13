# LPA model with taylor approximated exp function
# Solving using HomotopyContinuation
using HomotopyContinuation
using TaylorSeries
using LinearAlgebra
using RowEchelon: rref_with_pivots

# HC symbolic variables and parameters
HC_vars = @var L[0:nsteps] P[0:nsteps] A[0:nsteps]
HC_params = @var b cel cea cpa μl μa

function LPA_taylor!(out, unp1, un, p; exp_func=exp)
    b, cel, cea, cpa, μl, μa = p # parameters
    L1, P1, A1 = unp1 # u_n+1
    L, P, A = un # u_n

    out[1] = -L1 + b * A * exp_func(-cel * L - cea * A)
    out[2] = -P1 + (1 - μl) * L
    out[3] = -A1 + P * exp_func(-cpa * A) + (1 - μa) * A
    return nothing
end

function prolongate_LPA(vars, params; nsteps=3, order=2)
    out = zeros(Expression, 3)
    u = collect(zip(vars...))
    prolongations = Expression[]

    exp_taylor = convert(Taylor1{Rational{Int}}, taylor_expand(exp, 0; order))

    for i in 1:nsteps
        LPA_taylor!(out, u[i+1], u[i], params; exp_func=exp_taylor)
        append!(prolongations, out)
    end
    return prolongations
end

function create_homotopy_system(eqns, params, data_map)
    eqns_eval = HomotopyContinuation.evaluate(eqns, data_map)
    F = System(eqns_eval; variables=params)
    J = jacobian(F)

    # random large integer heuristic to find full rank sub matrix of Jacobian
    large_ints = rand(DiscreteUniform(10^3, 10^4), 6)
    J2 = HomotopyContinuation.evaluate(J, params => large_ints)
    # pivots are a set of independent columns
    _, pivots = rref_with_pivots(J2')

    System(eqns_eval[pivots]; variables=params)
end
