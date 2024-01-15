# eventually read from a file
# for now I'm running LPA_simulations.jl interactively
using DataFrames, CSV, Chain

param_names = [:b :cel :cea :cpa :μl :μa]
paramtuple(x) = NamedTuple{Tuple(param_names)}(Tuple(x))

function load_data(filename)
    evalparse(str::String) = eval(Meta.parse(str))
    @chain filename begin
        DataFrame(CSV.File(_))
        dropmissing # no real solutions found for even n
        transform([:true_parameters, :pred_parameters] .=> ByRow(evalparse) .=> [:true_p, :pred_p])
    end
end

#------------------------------------------------------------
df = load_data("tables/simulation_results.csv")

df_RMSE = @chain df begin
    select(:taylor_n, [:true_p, :pred_p] => ByRow((x,y) -> paramtuple((x-y).^2)) => AsTable)
    groupby(:taylor_n)
    combine(nrow => :count, Not(:taylor_n) .=> sqrt∘mean => x -> x*"_RMSE")
end

CSV.write("tables/error_table.csv", df_RMSE)
