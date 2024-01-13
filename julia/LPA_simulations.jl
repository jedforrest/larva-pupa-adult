# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)
using OrdinaryDiffEq
using Random, Distributions

include("parameters.jl")
include("LPA_homotopy_solve.jl")

# Original LPA model (exp version)
function LPA(du, u, p, t)
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    println(-cel * L - cea * A)
    println(-cpa * A)
    println()
    du[1] = b * A * exp(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp(-cpa * A) + (1 - μa) * A
end

function run_simulation(model, u0, params; steps=10)
    prob = DiscreteProblem(model, u0, (0., steps), params)
    sol = OrdinaryDiffEq.solve(prob, FunctionMap())
end

# simulation settings
steps = 10
prolongs = 3
te_order = 2#k # from parameters.jl

# Generates a matrix of values for LPA at t = 0, 1, ..., N
original_u0, original_params
LPA_sol = run_simulation(LPA, original_u0, original_params; steps)
LPAdata = LPA_sol[:,1:prolongs+1]
LPAdata = round.(Rational, LPAdata)

# Homotopy Solver
eqns = prolongate_LPA(HC_vars, HC_params; nsteps=3, order=te_order)
data = LPAdata[:,1:nsteps+1]

# create and solve homotopy system
# can't always find full rank matrix using heuristic
data_map = vcat(HC_vars...) => [data...]

eqns_eval = HomotopyContinuation.evaluate(eqns, data_map)

F = System(eqns_eval; variables=[HC_params...])
J = jacobian(F)
to_number.(J)

R, params = polynomial_ring(QQ, ["b", "cel", "cea", "cpa", "μl", "μa"])

B = QQ[
    1 2
    3 4
]

# TODO
# Convert J expr to MatSpace by replacing variables with QQ
# Then caclulate rank with rref_rational

typeof(B)
eltype(B)
RB=R.(B)
rref_rational(RB)
rref_rational(B)

qq_cea = params[3]
HomotopyContinuation.evaluate(J, cea => qq_cea)
HomotopyContinuation.coefficients(J[1], cea)
QQ.(J)

M, c = exponents_coefficients(J[1], [HC_params...])
M
c

# random large integer heuristic to find full rank sub matrix of Jacobian
large_ints = rand(DiscreteUniform(10^3, 10^4), 6)
J2 = HomotopyContinuation.evaluate(J, params => large_ints)
# pivots are a set of independent columns
_, pivots = rref_with_pivots(J2')

F = create_homotopy_system(eqns,
    [HC_params...],
    data_map
)

results = HomotopyContinuation.solve(F)
pred_params = real(solutions(results))[1]
