# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)

# Original parameters
# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
original_params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

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
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    du[1] = b * A * exp(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp(-cpa * A) + (1 - μa) * A
end

# Taylor series approximated LPA model
function LPA_taylor!(out, unp1, un, p; exp_func=exp)
    b, cel, cea, cpa, μl, μa = p # parameters
    L1, P1, A1 = unp1 # u_n+1
    L, P, A = un # u_n

    out[1] = -L1 + b * A * exp_func(-cel * L - cea * A)
    out[2] = -P1 + (1 - μl) * L
    out[3] = -A1 + P * exp_func(-cpa * A) + (1 - μa) * A
    return nothing
end

# Taylor series approximated LPA model with centering at each step
# Centering is determined by user inputed intervals (need not be exact)
function LPA_taylor_centred!(out, unp1, un, p; centres, order)
    b, cel, cea, cpa, μl, μa = p # parameters
    L1, P1, A1 = unp1 # u_n+1
    L, P, A = un # u_n

    # centres are 'initial guesses' for parameters within some interval
    l_centre = substitute(-cel * L - cea * A, p .=> centres)
    a_centre = substitute(-cpa * A, p .=> centres)

    # exp_shift(-cel * L - cea * A, l_centre; order)
    # exp_shift(-cpa * A, a_centre; order)

    out[1] = -L1 + b * A * exp_shift(-cel * L - cea * A, l_centre; order)
    out[2] = -P1 + (1 - μl) * L
    out[3] = -A1 + P * exp_shift(-cpa * A, a_centre; order) + (1 - μa) * A
    return nothing
end

exp_shift(x, c; order) = taylor_expand(exp, c; order)(x-c)

#----------------------------------------------
function prolongate_LPA(vars, params, param_intervals; nsteps=3, order=2)
    u = collect(zip(vars...))
    out = zeros(Num, 3)
    prolongations = Num[]

    # halfway between upper and lower bounds of intervals
    centres = map(p_i -> (p_i.ub - p_i.lb)/2, param_intervals)

    for i in 0:nsteps-1
        LPA_taylor_centred!(out, u[i+1], u[i], params; centres, order)
        append!(prolongations, out)
    end
    return prolongations
end



function run_simulation(model, u0, params; steps=10)
    prob = DiscreteProblem(model, u0, (0., steps), params)
    sol = OrdinaryDiffEq.solve(prob, FunctionMap())
end

#----------------------------------------------
