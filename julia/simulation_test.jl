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
# include("LPA_homotopy_solve.jl")

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
param_intervals = fill(RealInterval(0, 1), 6)
eqns = prolongate_LPA(Sym_vars, Sym_params, param_intervals; nsteps=steps, order=n)
eqns[1]


exp_a = taylor_expand(exp, 0; order=n)
exp_a(adult_eqn)

exp_a = taylor_expand(exp, substitute(adult_eqn, Sym_params .=> centres); order=n)
exp_a(adult_eqn)

# generate a matrix of values for LPA at t = 0, 1, ..., N
true_params = sampled_params
LPA_sol = run_simulation(LPA!, original_u0, true_params; steps)
data = vcat(LPA_sol[:,:]'...)

# determine rank of jacobian
data_map = Dict([L..., P..., A...] .=> data)
eqn_subs = substitute(eqns, data_map)
# J = Symbolics.jacobian(eqn_subs, Sym_params)
# rank(J)

eqns


p0 = rand(DiscreteUniform(0,10), 12)

eqns_subs = substitute(eqns, Dict([L..., P..., A...] .=> p0))
eqns_expand = Symbolics.expand.(eqn_subs)

eqn1 = eqns_subs[1]
eqn_vars = Symbolics.get_variables(eqn1)

function convert_to_HC_expression(eqn)
    eqn_vars = Symbolics.get_variables(eqn)
    eqn_syms = Symbolics.tosymbol.(eqn_vars)

    hc_vars = HomotopyContinuation.Variable.(eqn_syms)
    var_map = Dict(eqn_vars .=> hc_vars)
    sub = substitute(eqn, var_map)
end
eqn_1 = convert_to_HC_expression(eqn1)
typeof(eqn_1)

hc_eqns = convert_to_HC_expression.(eqn_expand)

F = System(hc_eqns)



pivots = [1,2,3,4,6,7]
F = System(eqns[pivots]; variables=[Sym_params...], parameters=[L..., P..., A...])

p0 = 100 .* rand(ComplexF64, 12)
result_p0 = HomotopyContinuation.solve(F, target_parameters = p0)



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
