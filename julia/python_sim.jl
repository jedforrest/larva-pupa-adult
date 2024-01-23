using OrdinaryDiffEq
using HomotopyContinuation
using Random, Distributions
using TaylorSeries
# using LinearAlgebra
using DataFrames, CSV

include("LPA_models.jl")
include("taylorseries_patch.jl")
#----------------------------------------------
using PythonCall
using CondaPkg

# if packages are missing install them in Julia with:
# using CondaPkg
# ]conda add some_package
@py import sys
@py import sympy
@py import importlib

sys.path.append(pwd()*"/../old_python_sims")

@py import dynamic_taylor
importlib.reload(dynamic_taylor) # for interactive code reloading

gensys = pyimport("dynamic_taylor" => "gensys")
all_poly_sys = pyimport("dynamic_taylor" => "all_poly_system_withapprox")

function gensys_to_hc_jl(py_sys, vars=nothing)
    # convert python gensys output to homotopy continuation expression in julia
    if isnothing(vars)
        @var b c_el c_ea c_pa mu_l mu_a
    end
    str_tuple = pyconvert(Tuple, py_sys)
    return eval.(Meta.parse.(str_tuple))
end

nsteps = 3
n = 5
py_sys = gensys(nsteps, 1, offmul=[1.0,1.0,1.0]);
hc_vars = @var b c_el c_ea c_pa mu_l mu_a
fL, fA = gensys_to_hc_jl(py_sys, hc_vars);

fL
fA

L_data, P_data, A_data = run_simulation(LPA!, original_u0, original_params; nsteps)
c_pa_approx = original_params[4]

py_sys = all_poly_sys(n, L_data, P_data, A_data, offmul=[1, 1, 1])
fL, fA = gensys_to_hc_jl(py_sys)

fL_sys = System(fL)
fA_sys = System(fA[1:2])
result = HomotopyContinuation.solve(fL_sys)
result = HomotopyContinuation.solve(fA_sys)
# solutions(result)
real_solutions(result)
original_params

function filter_solutions(result, paramnames, p_intervals)
    numreal = nreal(result)

    if numreal == 1
        return real_solutions(result)[1]
    elseif numreal > 1
        # filter
        p_int = collect(p_intervals[paramnames])
        for res in real_solutions(result)
            if all(in.(res, p_int))
                return res
            end
        end
    end
    return missing # default
end

#----------------------------------------------
# simulation settings
nsteps = 3 # simulation steps aka prolongs
Ntaylor = 7 # max taylor approx.
Nsims = 10 # sims per parameter set
interval_ranges = [0.05, 0.1, 0.2, 0.25, 0.5]

# presample all parameter sets
# parameter x is sampled from range [x ± p%]
create_intervals(x, p) = RealInterval((1 - p) * x, (1 + p) * x)
sample_interval(x_int) = rand(Uniform(x_int.lb, x_int.ub))

hc_vars = @var b c_el c_ea c_pa mu_l mu_a

param_names = [:b :c_el :c_ea :c_pa :mu_l :mu_a]
paramtuple(x) = NamedTuple{Tuple(param_names)}(Tuple(x))
#----------------------------------------------
param_intervals = create_intervals.(original_params, 0.05)
paramtuple(original_params)

# results file headings
df = DataFrame([
    "sim_num" => Int[],
    "interval_range" => Float64[],
    "taylor_n" => Int[],
    "L_num_real_solutions" => Int[],
    "L_pred_parameters" => Vector{Float64}[],
    "L_all_real_solutions" => Vector{Vector{Float64}}[],
    "A_num_real_solutions" => Int[],
    "A_pred_parameters" => Vector{Float64}[],
    "A_all_real_solutions" => Vector{Vector{Float64}}[],
    "sampled_parameters" => Vector{Float64}[],
    "offmul" => Vector{Float64}[],
    # "sampled_u0" => Vector{Float64}[],
])
allowmissing!(df, r"pred_parameters")

@time for I_range in interval_ranges
    print("[I=$I_range]")

    # parameter interval is sampled from range [x ± I_range%]
    param_intervals = paramtuple(create_intervals.(original_params, I_range))
    u0_intervals = original_u0 #create_intervals.(original_u0, I_range)

    # used to centre taylor exp in all_poly_sys
    offmul = rand(Uniform(1 - I_range, 1 + I_range), 6)

    for i in 1:Nsims
        print(" $i")

        # sample each new parameters and ICs from its own interval
        sampled_params = sample_interval.(collect(param_intervals))
        sampled_u0 = original_u0#sample_interval.(u0_intervals)

        # generate simulated data for LPA at t = 0, 1, ..., N
        L_data, P_data, A_data = run_simulation(LPA!, sampled_u0, sampled_params; nsteps)

        for n in 1:Ntaylor
            # println("N=$n]")

            # poly subsystems for L, P, A
            py_sys = all_poly_sys(n, L_data, P_data, A_data, offmul)
            L_sys, A_sys = gensys_to_hc_jl(py_sys)

            # Solve HC system and keep real solutions
            result = HomotopyContinuation.solve(L_sys)
            L_num_real = nreal(result)
            L_pred_params = filter_solutions(result, Symbol.(variables(L_sys)), param_intervals)
            L_all_real_solutions = real_solutions(result)

            # 3 eqns 2 vars - overdetermined
            result = HomotopyContinuation.solve(A_sys[1:2])
            A_num_real = nreal(result)
            A_pred_params = filter_solutions(result, Symbol.(variables(A_sys)), param_intervals)
            A_all_real_solutions = real_solutions(result)

            # save results
            push!(df, (i, I_range, n,
                L_num_real, L_pred_params, L_all_real_solutions,
                A_num_real, A_pred_params, A_all_real_solutions,
                sampled_params, offmul))
        end
    end
    println()
end

df

CSV.write("tables/simulation_results_python.csv", df)
