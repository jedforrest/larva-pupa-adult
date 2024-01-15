# eventually read from a file
# for now I'm running LPA_simulations.jl interactively
using DataFrames, CSV, Chain

param_names = [:b :cel :cea :cpa :μl :μa]
paramtuple(x) = NamedTuple{Tuple(param_names)}(Tuple(x))

function load_data(filename)
    df = DataFrame(CSV.File(filename))
    evalparse(str::String) = eval(Meta.parse(str))
    select(df, :n_taylor,
        [:true_parameters, :pred_parameters] .=> ByRow(evalparse) .=> [:true_p, :pred_p])
end

#------------------------------------------------------------
df = load_data("tables/simulation_results.csv")

df_RMSE = @chain df begin
    select(:n_taylor, [:true_p, :pred_p] => ByRow((x,y) -> paramtuple((x-y).^2)) => AsTable)
    groupby(:n_taylor)
    combine(Not(:n_taylor) .=> sqrt∘mean => x -> x*"_RMSE")
end

CSV.write("tables/error_table.csv", df_RMSE)
