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
include("taylorseries_patch.jl")

#----------------------------------------------
# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cea, cel, μl, cpa, μa
original_params = [6.598, 1.155e-2, 1.209e-2, 0.2055, 4.7e-3, 7.629e-3]

# simulation settings
steps = 3 # simulation steps aka prolongs
Ntaylor = 1 # max taylor approx.
Nsims = 1 # sims per parameter set
interval_ranges = [0.05]#, 0.1, 0.2, 0.25, 0.5]
centre_exps = true

# HC symbolic variables and parameters
sym_vars = @variables L[0:steps] P[0:steps] A[0:steps]
sym_params = @variables b cea cel μl cpa μa
sym_vars_flat = [L..., P..., A...]

# presample all parameter sets
# parameter x is sampled from range [x ± p%]
create_intervals(x, p) = RealInterval((1 - p) * x, (1 + p) * x)
sample_interval(x_int) = rand(Uniform(x_int.lb, x_int.ub))
#----------------------------------------------
# results file headings
df = DataFrame([
    "interval_range" => Float64[],
    "taylor_n" => Int[],
    "num_solutions" => Int[],
    "num_real_solutions" => Int[],
    "pred_parameters" => Vector{Number}[],
    "sampled_parameters" => Vector{Float64}[],
    "sampled_u0" => Vector{Float64}[],
])

@time for I_range in interval_ranges
    param_intervals = create_intervals.(original_params, I_range)
    u0_intervals = create_intervals.(original_u0, I_range)

    centres = if centre_exps
        # halfway between upper and lower bounds of parameter intervals
        map(p_i -> (p_i.ub - p_i.lb) / 2, param_intervals)
    else
        # no centering
        zeros(length(sym_params))
    end
    
    println("Centres for Taylor:")
    println(centres)
    println("===================")

    for n in 1:Ntaylor
        print("[I=$I_range, N=$n]")
        # prolongate with taylor approximation order n
        eqns = prolongate_LPA(sym_vars, sym_params, centres; nsteps=steps, order=n)

        for i in 1:Nsims
            print(" $i")

            # sample new parameters and ICs for given interval
            sampled_params = sample_interval.(param_intervals)
            sampled_u0 = sample_interval.(u0_intervals)

            println("Sampled parameters:")
            println(sampled_params)
            println("==================")
  
       	    println("Sampled IC:")
       	    println(sampled_u0)
       	    println("==================")             

            # generate simulated data for LPA at t = 0, 1, ..., N
            LPA_sol = run_simulation(LPA!, sampled_u0, sampled_params; steps)
            data = [LPA_sol[1, :]; LPA_sol[2, :]; LPA_sol[3, :]]
            println("Data:")
            println(data)
            println("================")

            println("Equations before data:")
            println(eqns)
            println("=================")

            # Create homotopy system from equations using known full rank column set

            eqns_with_data = substitute(eqns, Dict(sym_vars_flat .=> data))
            hc_eqns = convert_to_HC_expression.(eqns_with_data)
            pivots = [1,2,3,4,6,7]
            F = System(hc_eqns[pivots])

            # Solve HC system and keep real solutions (or first complex)
            println("System of polynomials to solve:")
            println(F)
            println("=========================")
            res = HomotopyContinuation.solve(F#;
#                show_progress = false#,
      #          stop_early_cb = r -> is_real(r) # stop once first real solution is found
            )
            pred_params = nreal(res) > 0 ? real_solutions(res)[1] : solutions(res)[1]
            println(real_solutions(res))
            # save results
            push!(df, (I_range, n, nsolutions(res), nreal(res), 
                pred_params, sampled_params, sampled_u0))
        end
        println()
    end
end

df

# write results to CSV file
CSV.write("tables/simulation_results_no_centre.csv", df)
