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

sys.path.append(pwd()*"\\old python sims")

importlib.reload(dynamic_taylor) # for interactive code reloading
gensys = pyimport("dynamic_taylor" => "gensys")
# pysolve = pyimport("homotopy_continuation" => "solve")

function gensys_to_hc_jl(py_sys, vars=nothing)
    # convert python gensys output to homotopy continuation expression in julia
    if isnothing(vars)
        @var b c_el c_ea c_pa mu_l mu_a
    end
    str_tuple = pyconvert(Tuple, py_sys)
    eqn_tuple = eval.(Meta.parse.(str_tuple))
    return System.(eqn_tuple)
end

nsteps = 3
n = 1
py_sys = gensys(nsteps, 1, offmul=[1.0,1.0,1.0]);
hc_vars = @var b c_el c_ea c_pa mu_l mu_a
fL, fA = gensys_to_hc_jl(py_sys, hc_vars);

fL
fA

fL_sys = System(fL)
fA_sys = System(fA)
result = HomotopyContinuation.solve(fL_sys)
# result = HomotopyContinuation.solve(fA_sys)
# solutions(result)
real_solutions(result)
original_params

#----------------------------------------------
# using Symbolics
# nsteps=3
# sym_vars = @variables L[0:nsteps] P[0:nsteps] A[0:nsteps]
# sym_params = @variables b cel cea cpa μl μa
# sym_vars_flat = [L..., P..., A...]
# n=1
# midpoints = [0.,0.,0.,0.,0.,0.]
# eqns = prolongate_LPA(sym_vars, sym_params, midpoints; nsteps, order=n)
# eqns = Symbolics.expand.(eqns)

#-----
# u = collect(zip(sym_vars...))
# out = zeros(Num, 3)
# prolongations = Num[]

# for i in 0:nsteps-1
#     LPA!(out, u[i], sym_params, nothing)
#     append!(prolongations, out .- u[i+1])
# end
# eqns = prolongations
# eqns = Symbolics.expand.(prolongations)

# eqns = prolongate_LPA(sym_vars, sym_params, midpoints; nsteps, order=n)
# eqns = Symbolics.expand.(eqns)

#-----

# sample new parameters and ICs for given interval
sampled_params = original_params#sample_interval.(param_intervals)
sampled_u0 = original_u0#sample_interval.(u0_intervals)

# generate simulated data for LPA at t = 0, 1, ..., N
L_data, P_data, A_data = run_simulation(LPA!, sampled_u0, sampled_params; nsteps)

eqns_3x3 = reshape(eqns, 3, :)
L_eqns = eqns_3x3[1,:]
P_eqns = eqns_3x3[2,:]
A_eqns = eqns_3x3[3,:]

data_map = [collect(L) .=> L_data; collect(P) .=> P_data; collect(A) .=> A_data]

L_subs = substitute(L_eqns, Dict(data_map))
P_subs = substitute(P_eqns, Dict(data_map))
A_subs = substitute(A_eqns, Dict(data_map))

fL_sys
fA_sys

#----------------------

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
#----------------------------------------------
param_intervals = create_intervals.(original_params, 0.05)


# results file headings
df = DataFrame([
    "data_set_num" => Int[],
    "interval_range" => Float64[],
    "taylor_n" => Int[],
    "num_solutions" => Int[],
    "num_real_solutions" => Int[],
    "pred_parameters" => Vector{Float64}[],
    "sampled_parameters" => Vector{Float64}[],
    "sampled_u0" => Vector{Float64}[],
])
allowmissing!(df, :pred_parameters)

@time for I_range in interval_ranges
    param_intervals = create_intervals.(original_params, I_range)
    u0_intervals = original_u0 #create_intervals.(original_u0, I_range)

    # halfway between upper and lower bounds of parameter intervals
    # midpoints = map(p_i -> (p_i.ub - p_i.lb) / 2, param_intervals)

    # sample new parameters and ICs for given interval
    sampled_params = sample_interval.(param_intervals)
    sampled_u0 = original_u0#sample_interval.(u0_intervals)

    interval = RealInterval(1 - I_range, 1 + I_range)
    offmul = rand(Uniform(interval.lb, interval.ub), 3)

    for n in 1:Ntaylor
        println("[I=$I_range, N=$n]")
        # prolongate with taylor approximation order n
        # eqns = prolongate_LPA(sym_vars, sym_params, midpoints; nsteps=steps, order=n)


        py_sys = gensys(nsteps, n, offmul=offmul);
        hc_vars = @var b c_el c_ea c_pa mu_l mu_a
        fL, fA = gensys_to_hc_jl(py_sys, hc_vars);
        F = fL # for now

        for i in 1:Nsims
            # print(" $i")



            # generate simulated data for LPA at t = 0, 1, ..., N
            # LPA_sol = run_simulation(LPA!, sampled_u0, sampled_params; steps)
            # data = [LPA_sol[1, :]; LPA_sol[2, :]; LPA_sol[3, :]]

            # Create homotopy system from equations using known full rank column set
            # eqns_with_data = substitute(eqns, Dict(sym_vars_flat .=> data))
            # hc_eqns = convert_to_HC_expression.(eqns_with_data)
            # pivots = [1,2,3,4,6,7]
            # F = System(hc_eqns[pivots])

            pred_params = missing

            # Solve HC system and keep real solutions
            result = HomotopyContinuation.solve(F)
            println("Num real: ", nreal(result))
            if nreal(result) == 1
                pred_params = real_solutions(result)[1]
            elseif nreal(result) > 1
                # filter
                p_int = param_intervals[[1,3,2]]
                for res in real_solutions(result)
                    if all(in.(res, p_int))
                        pred_params = res
                        break
                    end
                end

                # println("none found")
                # pred_params = missing # if nothing fits
            # else
                # pred_params = missing
            end

            # save results
            push!(df, (i, I_range, n, nsolutions(result), nreal(result),
                pred_params, sampled_params, sampled_u0))
        end
        println()
    end
end

df

CSV.write("tables/simulation_results_python.csv", df)
