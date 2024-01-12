# LPA model with taylor approximated exp function
# Solving using HomotopyContinuation
using HomotopyContinuation
using TaylorSeries
using LinearAlgebra
using RowEchelon: rref_with_pivots

function LPA_taylor(du, u, p, t; order=order)
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    exp_taylor = taylor_expand(exp, 0; order)

    du[1] = b * A * exp_taylor(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp_taylor(-cpa * A) + (1 - μa) * A
    return du
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

# Simulation script
nsteps = 3
vars = @var L[0:nsteps] P[0:nsteps] A[0:nsteps]
params = @var b cel cea cpa μl μa
eqns = prolongate_LPA(LPA_taylor, vars, params; nsteps, order=1)

# generated in other file
include("LPA_simulations.jl")
data = round.(Int, LPAdata[:,1:nsteps+1])

# create and solve system
# can't always find full rank matrix using heuristic
F = create_homotopy_system(eqns,
    [b, cel, cea, cpa, μl, μa],
    [L; P; A] => [data...]
)

results = HomotopyContinuation.solve(F)
res = solutions(results)

sampled_params
