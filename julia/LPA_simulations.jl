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
include("LPA_homotopy_solve.jl")

#----------------------------------------------
# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
original_params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

# simulation settings
steps = 3 # simulation steps aka prolongs
Ntaylor = 15 # max taylor approx.
Nsims = 100 # sims per parameter set

# HC symbolic variables and parameters
HC_vars = @var L[0:steps] P[0:steps] A[0:steps]
HC_params = @var b cel cea cpa μl μa

# presample all parameter sets
# parameter p is sampled from a range [p ± 25%]
sampled_params = [[rand(Uniform(0.75p, 1.25p)) for p in original_params] for i in 1:Nsims];
#----------------------------------------------
# results file headings
df = DataFrame([
    "param_set" => Int[],
    "taylor_n" => Int[],
    "num_solutions" => Int[],
    "num_real_solutions" => Int[],
    "true_parameters" => Vector{Float64}[],
    "pred_parameters" => Vector{Float64}[],
])
allowmissing!(df, :pred_parameters)

@time for n in 1:Ntaylor
    print("[Sim $n/$Ntaylor]")
    # prolongate with taylor approximation order n
    eqns = prolongate_LPA(HC_vars, HC_params; nsteps=steps, order=n)

    # Create homotopy system from equations using known full rank column set
    pivots = [1,2,3,4,6,7]
    F = System(eqns[pivots]; variables=[HC_params...], parameters=[L..., P..., A...])

    # Parameter homotopy method
    # Solve the system in the Complex domain to find maximum number of solutions
    # Use these solutions to find the Real solutions we want according to our data
    # Our 'parameters' are the data from the simulation
    p0 = 100 .* rand(ComplexF64, length(data))
    result_p0 = HomotopyContinuation.solve(F, target_parameters = p0)

    for i in 1:Nsims
        i%10 == 0 && print(" $i")

        # generate a matrix of values for LPA at t = 0, 1, ..., N
        true_params = sampled_params[i]
        LPA_sol = run_simulation(LPA!, original_u0, true_params; steps)
        data = vcat(LPA_sol[:,:]'...)

        # Solve HC system with data from simulation
        res = HomotopyContinuation.solve(
            F,
            solutions(result_p0);
            start_parameters = p0,
            target_parameters = data,
            transform_result = (r,p) -> real_solutions(r),
            flatten = true
        )
        pred_params = nreal(res) > 0 ? real_solutions(res)[1] : missing

        # save results
        push!(df, (i, n, nsolutions(res), nreal(res), true_params, pred_params))
    end
    println()
end

# write results to CSV file
CSV.write("tables/simulation_results.csv", df)
