# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)
using OrdinaryDiffEq
using HomotopyContinuation
using Random, Distributions
using TaylorSeries
using LinearAlgebra
using DataFrames, CSV

include("LPA_models.jl")

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
#----------------------------------------------
# results file headings
df = DataFrame([
    "param_set" => Int[],
    "n_taylor" => Int[],
    "true_parameters" => Vector{Float64}[],
    "pred_parameters" => Vector{Float64}[],
])

# presample all parameter sets
# parameter p is sampled from a range [p ± 25%]
sampled_params = [[rand(Uniform(0.75p, 1.25p)) for p in original_params] for i in 1:Nsims];

for i in 1:Nsims
    print("Simulation $i/$Nsims\nTaylor(n):")

    # generates a matrix of values for LPA at t = 0, 1, ..., N
    true_params = sampled_params[i]
    LPA_sol = run_simulation(LPA!, original_u0, true_params; steps)
    data_map = vcat((HC_vars[i] .=> sol[i,:] for i in eachindex(HC_vars))...)

    for n in 1:Ntaylor
        print(" $n")
        # prolongate with taylor approximation order n
        eqns = prolongate_LPA(HC_vars, HC_params; nsteps=steps, order=n)

        # create and solve homotopy system
        eqns_eval = HomotopyContinuation.evaluate(eqns, data_map...)
        F = System(eqns_eval; variables=[HC_params...]);

        # too expensive to run 'full_rank_subsystem(F)'
        # instead using known full rank column set
        pivots = [1,2,3,4,6,7]
        F2 = System(eqns_eval[pivots]);

        res = HomotopyContinuation.solve(F2)
        real_res = real_solutions(res)
        real_res = length(real_res) == 0 ? fill(0,6) : real_res[1]

        # save results
        push!(df, (i, n, true_params, real_res))
    end
    println("\n")
end

# write results to CSV file
CSV.write("tables/simulation_results.csv", df)
