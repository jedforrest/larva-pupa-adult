# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)
using OrdinaryDiffEq
using Random, Distributions
using TaylorSeries
using LinearAlgebra

include("parameters.jl")
include("LPA_homotopy_solve.jl")

#----------------------------------------------
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
#----------------------------------------------
# simulation settings
steps = 10
prolongs = 3
te_order = 3#k # from parameters.jl

# Generates a matrix of values for LPA at t = 0, 1, ..., N
original_u0, original_params
LPA_sol = run_simulation(LPA, original_u0, original_params; steps)
LPAdata = LPA_sol[:,1:prolongs+1]
LPAdata = Rational.(round.(LPAdata))

#----------------------------------------------
# Homotopy Solver
eqns = prolongate_LPA(HC_vars, HC_params; nsteps=3, order=te_order)
data = LPAdata[:,1:nsteps+1]

# create and solve homotopy system
data_map = vcat(HC_vars...) => [data...]
eqns_eval = HomotopyContinuation.evaluate(eqns, data_map)

F = System(eqns_eval; variables=[HC_params...]);

# Too expensive to solve; one known full rank subsystem is:
# [1, 2, 3, 4, 6, 7]
# F2 = full_rank_subsystem(F)
pivots = [1,2,3,4,6,7]
F2 = System(eqns_eval[pivots]);

res = HomotopyContinuation.solve(F2)
HomotopyContinuation.results(res, only_real=true)
sol = real(solutions(res, only_real=true))[1]

pred_params = real(solutions(res, only_real=true))[1]
