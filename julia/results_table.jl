# eventually read from a file
# for now I'm running LPA_simulations.jl interactively
using DataFrames, CSV, Chain

param_names = [:b :cel :cea :cpa :μl :μa]
paramtuple(x) = NamedTuple{Tuple(param_names)}(Tuple(x))
PTuple = NamedTuple{(:b, :c_el, :c_ea, :c_pa, :mu_l, :mu_a)}

function load_data(filename)
    evalparse(str::String) = eval(Meta.parse(str))
    @chain filename begin
        DataFrame(CSV.File(_))
        transform([:sampled_parameters, :pred_parameters] .=>
            ByRow(evalparse) .=> [:true_p, :pred_p])
    end
end

#------------------------------------------------------------
df = load_data("tables/simulation_results_24-01-24T12.csv")

names(df)

x=df.true_p[1]
y=df.pred_p[1]
x-y

df_RMSE = @chain df begin
    select(:taylor_n, :interval_range, [:true_p, :pred_p] => ByRow((x, y) -> ((x - y) .^ 2)) => AsTable)
    groupby([:taylor_n, :interval_range])
    combine(nrow => :count, Not(:taylor_n, :interval_range) .=> sqrt ∘ mean => x -> x * "_RMSE")
end

CSV.write("tables/error_table_results_python_100_random.csv", df_RMSE)
