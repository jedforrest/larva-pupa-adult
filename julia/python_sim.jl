using OrdinaryDiffEq
using HomotopyContinuation
using Random, Distributions
using DataFrames, CSV, Dates

include("LPA_models.jl")
#----------------------------------------------
using PythonCall
using CondaPkg

# if packages are missing install them in Julia with:
# using CondaPkg
# ]conda add some_package
@py import sys
@py import sympy
@py import importlib

sys.path.append(pwd()*"/../python_sims")
sys.path.append(pwd()*"./python_sims")

@py import dynamic_taylor
importlib.reload(dynamic_taylor) # for interactive code reloading
all_poly_sys = pyimport("dynamic_taylor" => "all_poly_system_withapprox")

function python_sys_to_hc(py_sys, vars=nothing)
    # convert python gensys output to homotopy continuation expression in julia
    str_tuple = pyconvert(Tuple, py_sys)
    return eval.(Meta.parse.(str_tuple))
end

#----------------------------------------------
# b c_ea c_el mu_l c_pa mu_a
original_params = NamedTuple([
    :b => 6.598,
    :c_ea => 1.155e-2,
    :c_el => 1.209e-2,
    :mu_l => 0.2055,
    :c_pa => 4.7e-3,        
    :mu_a => 7.629e-3,
])
# L0, P0, A0
original_u0 = [250., 5., 100.]

# HomotopyContinuation Variables
hc_vars = @var b c_ea c_el mu_l c_pa mu_a

# parameter tuples
PTuple = NamedTuple{(:b, :c_ea, :c_el, :mu_l, :c_pa, :mu_a)}

# parameter p is sampled from range [x ± r%]
function create_intervals(r, interval_centre=values(original_params))
    PTuple([RealInterval((1 - r) * p, (1 + r) * p) for p in interval_centre])
end

function sample_intervals(intervals)
    PTuple([rand(Uniform(int.lb, int.ub)) for int in values(intervals)])
end

function solve_and_filter_solutions(sys, intervals)
    result = HomotopyContinuation.solve(sys, show_progress = false)

    # only interested in real solutions
    all_real_sols = real_solutions(result)

    # filter solution(s) within intervals
    real_sol = filter_solutions(all_real_sols, intervals)

    return real_sol, all_real_sols
end

function filter_solutions(real_sols, intervals)
    for rsol in real_sols
        if all(in.(rsol, values(intervals)))
            return rsol
        end
    end
    return fill(0., length(intervals)) # default
end

#----------------------------------------------
Random.seed!(1234);

# simulation settings
nsteps = 3 # max simulation steps aka prolongs
Ntaylor = 10 # max taylor approx.
Nsims = 100 # sims per parameter set
interval_ranges = [0.05, 0.1, 0.2, 0.25, 0.5]

#----------------------------------------------
# results file headings
df = DataFrame([
    "sim_num" => Int[],
    "interval_range" => Float64[],
    "taylor_n" => Int[],
    "L_num_real_solutions" => Int[],
    "P_num_real_solutions" => Int[],
    "A_num_real_solutions" => Int[],
    "pred_parameters" => Any[],
    "sampled_parameters" => Any[],
    "param_intevals" => Any[],
    "L_all_real_solutions" => Vector{Vector{Float64}}[],
    "P_all_real_solutions" => Vector{Vector{Float64}}[],
    "A_all_real_solutions" => Vector{Vector{Float64}}[],
    "offmul" => Vector{Float64}[],
]);

# function solve_overdetermined_L(fL)
    
# function solve_L(f; relativeLogPresKer = 2)
    
    
#     degC3 = maximum(degree.(monomials(f[1])))-1
#     splitDeg(nn) = [div(nn,degC3) mod(nn,degC3)]
#     ltpoints(nn) = getLatticePoints(newtonPolytope(f[1]^splitDeg(nn)[1]*(c[1] + 200*c[2] - 10)^splitDeg(nn)[2],c))

#     # we estimate the degree of regularity using the approximated HS
#     t = Taylor1(Float64, size(f,1)*degC3)
#     HS = (1-t^degC3)^size(f,1)/((1-t^degC3)*(1-t)^2)
#     candReg = [i for i in range(1,10) if HS[i] <= 0][1] 

#     # In case our guess is not good, we keep increasing the deg of regularity
#     for k in range(1,degC3+1)
#         # We construct the monomials for the matrix
#         A₀ = ltpoints(1)
#         D = ltpoints(candReg)
#         E0 = ltpoints(candReg-1)
#         Ei = ltpoints(candReg-degC3)
#         E = pushfirst!([Ei for i = 1:size(f, 1)],E0);
#         # this creates the matrix from which we want to compute the kernel
#         Sylv = EigenvalueSolver.getRes(f, EigenvalueSolver.exptomon(D, c),  [EigenvalueSolver.exptomon(e, c) for e ∈ E[2:end]], c; complex = false)
#         # we use svd to compute the kernel
#         svdobj = GenericSVD.svd(transpose(Sylv), full = true);

#         sortedDiff = sort([[log10(svdobj.S[i])-log10(svdobj.S[i+1]),i] for i in range(1,size(svdobj.S,1)-1)])
#         if sortedDiff[end][2] == size(svdobj.S,1)-1 && sortedDiff[end][1] - sortedDiff[end-1][1] > relativeLogPresKer
#             N = transpose(svdobj.V[:, end:end])
#             ee = exponents.(EigenvalueSolver.exptomon(D, c));
#             # we compute the approximation by inverting the monomial map

#             c1pos = findall(x->x==[1,0,1], ee)[1]
#             c2pos = findall(x->x==[0,1,1], ee)[1]
#             candApprox = [N[c1pos]/N[2],N[c2pos]/N[2],N[2]/N[1]]

#             return candApprox[1], candApprox[2], candApprox[3]
#         else
#             candReg += 1
#         end
#     end
#     error("More than one solution, increase precision")
# end

    
namevars = "classic"

@time for I_range in interval_ranges
    print("[I=$I_range]")

    # parameter interval is sampled from range [p ± I_range%]
    param_intervals = create_intervals(I_range)

    # used to centre taylor exp in all_poly_sys
    offmul = rand(Uniform(1 - I_range, 1 + I_range), 3)

    for i in 1:Nsims
        print(" $i")

        # sample each new parameters and ICs from its own interval
        sampled_params = sample_intervals(param_intervals)
        sampled_u0 = original_u0 # may add u0 sampling

        # generate simulated data for LPA at t = 0, 1, ..., N
        L_data, P_data, A_data = run_simulation(LPA!, sampled_u0, sampled_params; nsteps)

        for n in 1:Ntaylor
            # poly subsystems for L, P, A
            py_sys = all_poly_sys(n, L_data, P_data, A_data, offmul, namevars)
            L_sys, P_sys, A_sys = python_sys_to_hc(py_sys)

            # Solve HC system and keep real solutions
            L_sym_vars = Symbol.(variables(L_sys))
            L_pred, L_all_real = solve_and_filter_solutions(L_sys, param_intervals[L_sym_vars])

            # 3 eqns 1 var - overdetermined
            P_sym_vars = Symbol.(variables(P_sys))
            P_pred, P_all_real = solve_and_filter_solutions([P_sys[1]], param_intervals[P_sym_vars])

            # 3 eqns 2 vars - overdetermined
            A_sym_vars = Symbol.(variables(A_sys))
	    A_pred, A_all_real = solve_and_filter_solutions(A_sys[1:2], param_intervals[A_sym_vars])
            # save results (as named tuples)
            pred_params = NamedTuple([
                L_sym_vars .=> L_pred;
                P_sym_vars .=> P_pred;
                A_sym_vars .=> A_pred;
            ])

            push!(df, (i, I_range, n,
                length(L_all_real),
                length(P_all_real),
                length(A_all_real),
                pred_params,
                sampled_params,
                param_intervals,
                L_all_real,
                P_all_real,
                A_all_real,
                offmul))
        end
    end
    println()
end

df

# UTC year-month-day-hour
timestamp() = Dates.format(now(UTC), "yy-mm-ddTHH")
CSV.write("tables/simulation_results_$(timestamp()).csv", df)
