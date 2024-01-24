# eventually read from a file
# for now I'm running LPA_simulations.jl interactively
using DataFrames, CSV, Chain, Statistics

param_names = [:b :cea :cel :μl :cpa :μa]
paramtuple(x) = NamedTuple{Tuple(param_names)}(Tuple(x))
PTuple = NamedTuple{(:b, :c_ea, :c_el, :mu_l, :c_pa, :mu_a)}

function load_data(filename)
    evalparse(str::String) = collect(eval(Meta.parse(str)))
    @chain filename begin
        DataFrame(CSV.File(_))
        transform([:sampled_parameters, :pred_parameters] .=>
            ByRow(evalparse) .=> [:true_p, :pred_p])
    end
end

#------------------------------------------------------------
df = load_data("tables/simulation_results_24-01-24T20.csv")

names(df)

x=df.true_p[1]
y=df.pred_p[1]
x-y

df_mean = @chain df begin
    select(:taylor_n, :interval_range, [:true_p, :pred_p] => ByRow((x, y) -> (abs.(x - y)./x)) => AsTable)
    groupby([:taylor_n, :interval_range])
    combine(nrow => :count, Not(:taylor_n, :interval_range) .=> mean => x -> x * "_mean")
end

CSV.write("tables/error_table_results_means.csv", df_mean)
