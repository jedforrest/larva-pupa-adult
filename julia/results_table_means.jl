# eventually read from a file
# for now I'm running LPA_simulations.jl interactively
using DataFrames, CSV, Chain, Statistics, Printf, Latexify

function load_data(filename)
    evalparse(str::String) = eval(Meta.parse(str))
    @chain filename begin
        DataFrame(CSV.File(_))
        @aside println.(names(_)) # column names
        transform([:sampled_parameters, :pred_parameters] .=>
            ByRow(evalparse) .=> [:true_p, :pred_p])
        transform([:sim_num, :taylor_n, :interval_range] .=> ByRow(Symbol); renamecols=false)
    end
end

#------------------------------------------------------------
df = load_data("tables/simulation_results_fixed_center_24-01-24T21.csv");

# check that the parameter orders of true_p and pred_p are the same
@assert keys(df.true_p[1]) == keys(df.pred_p[1])
PTuple = NamedTuple{keys(df.true_p[1])}

# reduce df to only the columns used for easier processing and debugging
# also 'collect' numbers only from param tuples
select!(df, :sim_num, :taylor_n, :interval_range,
    [:true_p, :pred_p] .=> ByRow(collect) .=> [:true_p, :pred_p])

# relative error
rel_error(x, y) = abs.(x - y)./x

df_mean_median = @chain df begin
    select(:taylor_n, :interval_range, [:true_p, :pred_p] => ByRow((x, y) -> PTuple(rel_error(x,y))) => AsTable)
    groupby([:taylor_n, :interval_range])
    combine(
        Not(:taylor_n, :interval_range) .=> mean => x -> "Mean($x)",
        Not(:taylor_n, :interval_range) .=> median => x -> "Med($x)")
    rename(:taylor_n => "Taylor(n)", :interval_range => "Interval")
end

CSV.write("tables/error_table_means_medians.csv", df_mean_median)

# write tables in latex output
open("tables/latex_tables.tex", "w") do f
    write(f, "Relative Mean and Median (Scientific Notation)\n")
    write(f, latexify(df_mean_median; env=:table, fmt="%.1e"))
    write(f, "Relative Mean and Median (Percentage Notation)\n")
    write(f, latexify(df_mean_median; env=:table, fmt=x->@sprintf("%.1f", 100x)))
end
