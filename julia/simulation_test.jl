# Single run simulation for testing purposes
# Use in interactive mode

using OrdinaryDiffEq
using HomotopyContinuation
using Random, Distributions
using TaylorSeries
using LinearAlgebra
using DataFrames, CSV
using Symbolics
using RowEchelon

include("LPA_models.jl")
include("taylorseries_patch.jl")

LPA_sol = run_simulation(LPA!, original_u0, original_params; steps=3)

#----------------------------------------------
# simulation settings
nsteps = 3 # simulation steps aka prolongs

# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
original_params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

# HC symbolic variables and parameters
# HC_vars = @var L[0:steps] P[0:steps] A[0:steps]
# HC_params = @var b cel cea cpa μl μa
sym_vars = @variables L[0:steps] P[0:steps] A[0:steps]
sym_params = @variables b cel cea cpa μl μa
sym_vars_flat = [L..., P..., A...]

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
n = 5
true_params = sampled_params

# prolongate with taylor approximation order n
p_interval(p) = RealInterval(rationalize.((0.75p, 1.25p))...)
param_intervals = p_interval.(original_params)
midpoints = map(p_i -> (p_i.ub - p_i.lb) / 2, param_intervals)

eqns = prolongate_LPA(sym_vars, sym_params, midpoints; nsteps, order=n);

# generate a matrix of values for LPA at t = 0, 1, ..., N
LPA_sol = run_simulation(LPA!, original_u0, true_params; steps)
data = vcat(LPA_sol[:, :]'...)

eqns_subs = substitute(eqns, Dict(sym_vars_flat .=> data));

J = Symbolics.jacobian(eqns_subs, sym_params);
# J = Symbolics.sparsejacobian(eqns_subs, sym_params)

largeints = 2 .^ (1:6)
J2 = substitute(J, Dict(sym_params .=> largeints))
J3 = map(j -> Float64(j.val), J2')

res, piv = rref_with_pivots(J3)
res
piv

hc_eqns = convert_to_HC_expression.(eqns_subs)
pivots = [1, 2, 3, 4, 6, 7]
F = System(hc_eqns[pivots])
J = jacobian(F)
J2 = to_number.(J)


res = HomotopyContinuation.solve(F)
real_solutions(res)[1]
