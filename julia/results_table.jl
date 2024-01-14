# eventually read from a file
# for now I'm running LPA_simulations.jl interactively
using Plots

param_names = string.(HC_params)

# from other files
original_params
pred_params

original_params .- pred_params

# simulate with predicted parameters
# compare with original simulation
steps = 20
u0 = original_u0
lpa_sol_original = run_simulation(LPA, u0, original_params; steps)
lpa_sol_pred = run_simulation(LPA, u0, pred_params; steps)

plot(plot(lpa_sol_original), plot(lpa_sol_pred); layout=(:,1))
