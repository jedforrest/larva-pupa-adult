# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)
using OrdinaryDiffEq
using Random, Distributions

include("parameters.jl")
include("LPA_homotopy_solve.jl")

# Original LPA model (exp version)
function LPA(du, u, p, t)
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    println(-cel * L - cea * A)
    println(-cpa * A)
    println()
    du[1] = b * A * exp(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp(-cpa * A) + (1 - μa) * A
end

function run_simulation(model, u0, params; steps=10)
    prob = DiscreteProblem(model, u0, (0., steps), params)
    sol = OrdinaryDiffEq.solve(prob, FunctionMap())
end

# simulation settings
steps = 10
prolongs = 3
te_order = 1#k # from parameters.jl

# Generates a matrix of values for LPA at t = 0, 1, ..., N
original_u0, original_params
LPA_sol = run_simulation(LPA, original_u0, original_params; steps)
LPAdata = LPA_sol[:,1:prolongs+1]
LPAdata = round.(Rational, LPAdata)

# Homotopy Solver
eqns = prolongate_LPA(HC_vars, HC_params; nsteps=3, order=te_order)
data = LPAdata[:,1:nsteps+1]

# create and solve homotopy system
# can't always find full rank matrix using heuristic
data_map = vcat(HC_vars...) => [data...]

eqns_eval = HomotopyContinuation.evaluate(eqns, data_map)

#----------------------------------------------
F = System(eqns_eval; variables=[HC_params...])
J = jacobian(F)
j1 = J[1,1]

function expression_to_polyn_ring()

end

j1
variables(j1)
HomotopyContinuation.coefficients(j1, variables(j1))

p = ["b", "cel", "cea", "cpa", "μl", "μa"]
R, params_qq = polynomial_ring(QQ, p.*"_qq")
b_qq, cel_qq, cea_qq, cpa_qq, μl_qq, μa_qq = params_qq

string.(J)
replacements = p .=> params_qq
Jstr = replace.(string.(J), replacements...)
# Jr = replace(string(J), replacements...)
Jpar = Meta.parse.(Jstr)
typeof.(Jpar)
Jeval = eval.(Jpar)
typeof.(Jeval)
Jeval[1]
rank(Jeval)

Jstr = string(J)
Jstr = replace(Jstr, "Expression" => "R")
ex = Meta.parse(Jstr)
Jeval = eval(ex)
typeof(Jeval)
Jeval[1,1]
Jeval[1,1] |> typeof

rank(Jeval)
J2 = transpose(Jeval)
r, A, d = rref_rational(J2)
r
A
Ad = A/d
is_rref(Ad)

Ad[:,1]

Ad
J2

# TODO
# determine columns in Ad which are pivots
# then select columns from F
# try to solve F in AA with solve_rational
# if not, return to HC to solve


#----------------------------------------------

F = create_homotopy_system(eqns,
    [HC_params...],
    data_map
)

results = HomotopyContinuation.solve(F)
pred_params = real(solutions(results))[1]
