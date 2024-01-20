# eventually read from a file
# for now I'm running LPA_simulations.jl interactively
using DataFrames, CSV, Chain

param_names = [:b :cel :cea :cpa :μl :μa]
paramtuple(x) = NamedTuple{Tuple(param_names)}(Tuple(x))

function load_data(filename)
    evalparse(str::String) = eval(Meta.parse(str))
    @chain filename begin
        DataFrame(CSV.File(_))
        subset(:taylor_n => n -> .!iseven.(n)) # no real solutions found for even n
        transform([:sampled_parameters, :pred_parameters] .=> ByRow(evalparse) .=> [:true_p, :pred_p])
    end
end

#------------------------------------------------------------
df = load_data("tables/simulation_results_centering.csv")

df_RMSE = @chain df begin
    select(:taylor_n, :interval_range, [:true_p, :pred_p] => ByRow((x, y) -> paramtuple((x - y) .^ 2)) => AsTable)
    groupby([:taylor_n, :interval_range])
    combine(nrow => :count, Not(:taylor_n, :interval_range) .=> sqrt ∘ mean => x -> x * "_RMSE")
end

CSV.write("tables/error_table_centering_odd_only.csv", df_RMSE)
