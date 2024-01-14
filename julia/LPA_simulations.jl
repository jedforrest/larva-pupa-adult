# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)
using OrdinaryDiffEq
using HomotopyContinuation
using Random, Distributions
using TaylorSeries
using LinearAlgebra
using DelimitedFiles

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
#----------------------------------------------
# simulation settings
steps = 3 # simulation steps aka prolongs
Ntaylor = 7 # max taylor approx.
Nsims = 100 # sims per parameter set

# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
original_params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

# HC symbolic variables and parameters
HC_vars = @var L[0:steps] P[0:steps] A[0:steps]
HC_params = @var b cel cea cpa μl μa

# results file headings
open("simulation_results.csv", "a") do io
    write(io, "i,n,")
    write(io, "b,cel,cea,cpa,μl,μa,")
    write(io, "b_pred,cel_pred,cea_pred,cpa_pred,μl_pred,μa_pred")
    write(io, "\n")
end

for i in 1:Nsims
    # parameter p is sampled from a range [p ± 25%]
    sampled_params = [rand(Uniform(0.75p, 1.25p)) for p in original_params]

    # generates a matrix of values for LPA at t = 0, 1, ..., N
    LPA_sol = run_simulation(LPA, original_u0, sampled_params; steps)
    LPAdata = Rational.(round.(LPA_sol[:,:]))

    for n in 1:Ntaylor
        # prolongate with taylor approximation order n
        eqns = prolongate_LPA(HC_vars, HC_params; nsteps=steps, order=n)

        # create and solve homotopy system
        data_map = vcat(HC_vars...) => [data...]
        eqns_eval = HomotopyContinuation.evaluate(eqns, data_map)
        F = System(eqns_eval; variables=[HC_params...]);

        # too expensive to run 'full_rank_subsystem(F)'; instead using known full rank column set
        pivots = [1,2,3,4,6,7]
        F2 = System(eqns_eval[pivots]);

        res = HomotopyContinuation.solve(F2)
        real_res = real_solutions(res)
        real_res = length(real_res) == 0 ? fill(0,6) : real_res[1]

        # save results
        open("simulation_results.csv", "a") do io
            write(io, "$i,$n,")
            writedlm(io, [real_res; sampled_params]', ',')
        end
    end
end
