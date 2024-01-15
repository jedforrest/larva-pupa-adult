# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)

# Original parameters
# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
original_params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

#----------------------------------------------
# Original LPA model (exp version)
function LPA!(du, u, p, t)
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    du[1] = b * A * exp(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp(-cpa * A) + (1 - μa) * A
end

function LPA_taylor!(out, unp1, un, p; exp_func=exp)
    b, cel, cea, cpa, μl, μa = p # parameters
    L1, P1, A1 = unp1 # u_n+1
    L, P, A = un # u_n

    out[1] = -L1 + b * A * exp_func(-cel * L - cea * A)
    out[2] = -P1 + (1 - μl) * L
    out[3] = -A1 + P * exp_func(-cpa * A) + (1 - μa) * A
    return nothing
end

#----------------------------------------------
function run_simulation(model, u0, params; steps=10)
    prob = DiscreteProblem(model, u0, (0., steps), params)
    sol = OrdinaryDiffEq.solve(prob, FunctionMap())
end

function prolongate_LPA(vars, params; nsteps=3, order=2)
    out = zeros(Expression, 3)
    u = collect(zip(vars...))
    prolongations = Expression[]

    exp_taylor = convert(Taylor1{Rational{Int}}, taylor_expand(exp, 0; order))

    for i in 1:nsteps
        LPA_taylor!(out, u[i+1], u[i], params; exp_func=exp_taylor)
        append!(prolongations, out)
    end
    return prolongations
end

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
