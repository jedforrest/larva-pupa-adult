# LPA model with taylor approximated exp function
# Solving using HomotopyContinuation
using HomotopyContinuation
using TaylorSeries
using LinearAlgebra
using RowEchelon: rref_with_pivots

# order of taylor series expansion
# order = 1
# exp taylor approximation with rational coefficients
# exp_taylor = convert(Taylor1{Rational{Int}}, taylor_expand(exp, 0; order))
exp_taylor = taylor_expand(exp, 0; order=2)

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

    # use jacobian to find full rank sub matrix
    J = jacobian(F)
    J_subm = find_submatrix(J, params)

    System([J_subm...])
end

function find_submatrix(J, params)
    # random large integer heuristic approach
    large_ints = rand(DiscreteUniform(10^2, 10^3), 6)
    J2 = HomotopyContinuation.evaluate(J, params => large_ints)

    a, pivots = rref_with_pivots(J2')
    J[pivots, :]
end

# Simulation script
nsteps = 5
vars = @var L[0:nsteps] P[0:nsteps] A[0:nsteps]
params = @var b cel cea cpa μl μa
eqns = prolongate_LPA(LPA_taylor, vars, params; nsteps, order=5)

# generated in other file
include("LPA_simulations.jl")
data = round.(Int, LPAdata[:,1:nsteps+1])

F = create_homotopy_system(eqns,
    [b, cel, cea, cpa, μl, μa],
    [L; P; A] => [data...]
)

# Unable to get solutions from this
HomotopyContinuation.solve(F)
