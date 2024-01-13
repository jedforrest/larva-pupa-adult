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
te_order = 2

# Generates a matrix of values for LPA at t = 0, 1, ..., N
original_u0, original_params
LPA_sol = run_simulation(LPA, original_u0, original_params; steps)
LPAdata = LPA_sol[:,1:prolongs+1]

# Homotopy Solver
eqns = prolongate_LPA(LPA_taylor; nsteps=3, order=te_order)
data = LPAdata[:,1:nsteps+1]

# create and solve homotopy system
# can't always find full rank matrix using heuristic
data_map = vcat(HC_vars...) => [data...]
F = create_homotopy_system(eqns,
    [b, cel, cea, cpa, μl, μa],
    data_map
)

results = HomotopyContinuation.solve(F)
pred_params = real(solutions(results))[1]
