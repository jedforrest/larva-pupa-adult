# Single run simulation for testing purposes
# Use in interactive mode

using OrdinaryDiffEq
using HomotopyContinuation
using Random, Distributions
using TaylorSeries
using LinearAlgebra
using DataFrames, CSV
using Symbolics

include("LPA_models.jl")
include("taylorseries_patch.jl")

#----------------------------------------------
# simulation settings
steps = 3 # simulation steps aka prolongs

# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
original_params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

# HC symbolic variables and parameters
# HC_vars = @var L[0:steps] P[0:steps] A[0:steps]
# HC_params = @var b cel cea cpa μl μa
Sym_vars = @variables L[0:steps] P[0:steps] A[0:steps]
Sym_params = @variables b cel cea cpa μl μa

# presample all parameter sets
# parameter p is sampled from a range [p ± 25%]
sampled_params = [rand(Uniform(0.75p, 1.25p)) for p in original_params]
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

i = 1
n = 1

# prolongate with taylor approximation order n
p_interval(x) = RealInterval(0, round(10. ^ ceil(log10(x)), sigdigits=1))
param_intervals = p_interval.(original_params)
eqns = prolongate_LPA(Sym_vars, Sym_params, param_intervals; nsteps=steps, order=n)

data

subs = substitute(eqns, Dict(Sym_params .=> true_params))

eqns_subs = substitute(eqns, Dict([L..., P..., A...] .=> data))

hc_eqns = convert_to_HC_expression.(eqns_subs)

# I shouldn't be getting negative numbers in my solution
# this wasn't happening before
# something wrong in the shift algorithm

F = System(hc_eqns)

pivots = [1,2,3,4,6,7]
F = System(hc_eqns[pivots])

p0 = 100 .* rand(Float64, 12)
# result_p0 = HomotopyContinuation.solve(F, target_parameters = p0)
result_p0 = HomotopyContinuation.solve(F)
res = real_solutions(result_p0)[1]

verify_solution(eqns, sym_vars, data, Sym_params, res)

# generate a matrix of values for LPA at t = 0, 1, ..., N
true_params = sampled_params
LPA_sol = run_simulation(LPA!, original_u0, true_params; steps)
data = vcat(LPA_sol[:,:]'...)

# # determine rank of jacobian
# data_map = Dict([L..., P..., A...] .=> data)
# eqn_subs = substitute(eqns, data_map)
# # J = Symbolics.jacobian(eqn_subs, Sym_params)
# # rank(J)

F = System(hc_eqns[pivots])
HomotopyContinuation.solve(F)

res = HomotopyContinuation.solve(
    F,
    solutions(result_p0);
    start_parameters = p0,
    target_parameters = data,
    transform_result = (r,p) -> real_solutions(r),
    flatten = true
)

pivots = [1,2,3,4,6,7]
F = System(eqns[pivots]; variables=[Sym_params...], parameters=[L..., P..., A...])

p0 = 100 .* rand(ComplexF64, 12)


result_p0 = HomotopyContinuation.solve(F, target_parameters = p0)

verify_solution(eqns, sym_vars, data, Sym_params, true_params)


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
