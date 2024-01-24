# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)

#----------------------------------------------
function lagrange_error_bound(f, c; a=0, tol=1e-3, max_k=19)
    # for function f(x) x∈[a,c] determine order k such that erro |ε| < tol
    # See: https://en.wikipedia.org/wiki/Taylor%27s_theorem#Estimates_for_the_remainder
    M = ceil(f(c))
    for k in 1:max_k
        ε = M*(abs(c-a)^(k+1))/factorial(k+1)
        if ε < tol
            return k, ε
        end
    end
    error("max k=$max_k reached")
end

c = 4
k, ε = lagrange_error_bound(exp, c; a=0, tol=1e-3)

#----------------------------------------------
# Original LPA model (exp version)
function LPA!(du, u, p, t)
    b, c_el, c_ea, mu_l, c_pa, mu_a = p
    L, P, A = u

    du[1] = b * A * exp(-c_el * L - c_ea * A)
    du[2] = (1 - mu_l) * L
    du[3] = P * exp(-c_pa * A) + (1 - mu_a) * A
end

# Taylor series approximated LPA model
function LPA_taylor!(out, unp1, un, p; order)
    b, c_el, c_ea, mu_l, c_pa, mu_a = p
    L1, P1, A1 = unp1 # u_n+1
    L, P, A = un # u_n

    exp_taylor = taylor_expand(exp, 0; order)

    # implicit form (output should be ≈ 0)
    out[1] = -L1 + b * A * exp_taylor(-c_el * L - c_ea * A)
    out[2] = -P1 + (1 - mu_l) * L
    out[3] = -A1 + P * exp_taylor(-c_pa * A) + (1 - mu_a) * A
    return nothing
end

# Taylor series approximated LPA model with centering at each step
# Centering is determined by user inputed intervals (need not be exact)
function LPA_taylor_centred!(out, unp1, un, p, midpoints; order)
    b, c_el, c_ea, mu_l, c_pa, mu_a = p
    L1, P1, A1 = unp1 # u_n+1
    L, P, A = un # u_n

    # centres are 'initial guesses' for parameters within some interval
    l_midpoints = substitute(-c_el * L - c_ea * A, p .=> midpoints)
    a_midpoints = substitute(-c_pa * A, p .=> midpoints)

    # implicit form (output should be ≈ 0)
    out[1] = -L1 + b * A * exp_shift(-c_el * L - c_ea * A, l_midpoints; order)
    out[2] = -P1 + (1 - mu_l) * L
    out[3] = -A1 + P * exp_shift(-c_pa * A, a_midpoints; order) + (1 - mu_a) * A
    return nothing
end

exp_shift(x, c; order) = taylor_expand(exp, c; order)(x - c)

#----------------------------------------------
function prolongate_LPA(vars, params, midpoints; nsteps=3, order)
    u = collect(zip(vars...))
    out = zeros(Num, 3)
    prolongations = Num[]

    for i in 0:nsteps-1
        LPA_taylor_centred!(out, u[i+1], u[i], params, midpoints; order)
        append!(prolongations, out)
    end
    return prolongations
end

function convert_to_HC_expression(eqn)
    eqn_vars = Symbolics.get_variables(eqn)
    hc_vars = HomotopyContinuation.Variable.(Symbolics.tosymbol.(eqn_vars))
    sub = substitute(eqn, Dict(eqn_vars .=> hc_vars))
end

function run_simulation(model, u0, params; nsteps=10)
    prob = DiscreteProblem(model, u0, (0., nsteps), params)
    sol = OrdinaryDiffEq.solve(prob, FunctionMap())
    return [sol[i, :] for i in 1:length(u0)]
end

function verify_solution(eqns, sym_vars, data, sym_params, sol_params)
    # verify that the parameters found are a solution to the equations + data
    # exact solution would return all zeros
    res = substitute(eqns, Dict(sym_vars .=> data))
    res = substitute(res, Dict(sym_params .=> sol_params))
    return res
end

#----------------------------------------------
