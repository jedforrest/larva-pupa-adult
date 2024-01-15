# Single run simulation for testing purposes
# Use in interactive mode

using OrdinaryDiffEq
using HomotopyContinuation
using Random, Distributions
using TaylorSeries
using LinearAlgebra
using DataFrames, CSV

include("LPA_models.jl")
include("LPA_homotopy_solve.jl")

#----------------------------------------------
# simulation settings
steps = 3 # simulation steps aka prolongs


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
    "taylor_n" => Int[],
    "num_solutions" => Int[],
    "num_real_solutions" => Int[],
    "true_parameters" => Vector{Float64}[],
    "pred_parameters" => Vector{Float64}[],
])

# presample all parameter sets
# parameter p is sampled from a range [p ± 25%]
sampled_params = [rand(Uniform(0.75p, 1.25p)) for p in original_params]

# generates a matrix of values for LPA at t = 0, 1, ..., N
true_params = original_params
LPA_sol = run_simulation(LPA!, original_u0, true_params; steps)
data_map = vcat((HC_vars[i] .=> LPA_sol[i,:] for i in eachindex(HC_vars))...)

# prolongate with taylor approximation order n
n = 9
eqns = prolongate_LPA(HC_vars, HC_params; nsteps=steps, order=n)

# create and solve homotopy system
# using known full rank column set
pivots = [1,2,3,4,6,7]
F = System(eqns[pivots]; variables=[HC_params...], parameters=[L..., P..., A...])

data = vcat(LPA_sol[:,:]'...)

# res = parameter_homotopy_solve(F, data)
    # Parameter homotopy method
    # Solve the system in the Complex domain to find maximum number of solutions
    # Use these solutions to find the Real solutions we want according to our data
    # Our 'parameters' are the data from the simulation
    p0 = 100 .* rand(ComplexF64, length(data))
    @time result_p0 = HomotopyContinuation.solve(F, target_parameters = p0)

    @time HomotopyContinuation.solve(
        F,
        solutions(result_p0);
        start_parameters = p0,
        target_parameters = data,
        transform_result = (r,p) -> real_solutions(r),
        flatten = true
    )

res
real_solutions(res)[1]
nreal(res)
nresults(res)

res

pred_params = nreal(res) > 0 ? real_solutions(res)[1] : missing


# save results
push!(df, (1, n, nsolutions(res), nreal(res), true_params, pred_params))



# write results to CSV file
CSV.write("tables/simulation_results.csv", df)


@var x y a b
F = System([x^2 - a, x * y - a + b]; variables=[x,y], parameters =[a,b])
start_solutions = [[1, 1]]
HomotopyContinuation.solve(F, start_solutions; start_parameters=[1, 0], target_parameters=[2, 5])
