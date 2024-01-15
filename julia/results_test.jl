# Test predicted parameters in the original LPA sim algorithm
using DataFrames, CSV, Chain
using Plots

include("LPA_models.jl")

# results_table.jl
df = load_data("tables/simulation_results.csv")

idx = 90
true_params = df.true_p[idx]
pred_params = df.pred_p[idx]

# simulate with predicted parameters
# compare with original simulation
steps = 20
u0 = original_u0
lpa_sol_true = run_simulation(LPA!, u0, true_params; steps)
lpa_sol_pred = run_simulation(LPA!, u0, pred_params; steps)

labels = ["L" "P" "A"]
plot(
    plot(lpa_sol_true; labels, legend=:right),
    plot(lpa_sol_pred; labels, legend=:right);
    layout=(:,1)
)
