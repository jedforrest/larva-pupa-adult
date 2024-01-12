# LPA model with taylor approximated exp function
using AbstractAlgebra
using LinearAlgebra
using TaylorSeries

function exp_taylor_polyn(R::Ring, series_order::Int)
    # TaylorSeries is not compatible with PolynomialRings
    # so convert series approximation to polynomial with same coefficients
    exp_taylor = convert(Taylor1{Rational{Int}}, taylor_expand(exp, 0; order=series_order))
    coeffs = getcoeff.(exp_taylor, 0:series_order)
    return polynomial(R, coeffs)
end

function LPA_taylor(out, du, u, p; exp_func=exp)
    b, cel, cea, cpa, μl, μa = p # parameters
    L1, P1, A1 = du # u_n+1
    L, P, A = u # u_n

    out[1] = -L1 + b * A * exp_func(-cel * L - cea * A)
    out[2] = -P1 + (1 - μl) * L
    out[3] = -A1 + P * exp_func(-cpa * A) + (1 - μa) * A
    return nothing
end

function prolongate_LPA(; nsteps=3, order=1)
    S, params = polynomial_ring(QQ, ["b", "cel", "cea", "cpa", "μl", "μa"])
    R, l, p, a = polynomial_ring(S, :L => 0:3, :P => 0:3, :A => 0:3)
    out = zeros(R, 3)
    prolongations = []

    exp_taylor = exp_taylor_polyn(R, order)

    for i in 1:nsteps
        u = [l[i], p[i], a[i]]
        du = [l[i+1], p[i+1], a[i+1]]
        LPA_taylor(out, du, u, params; exp_func=exp_taylor)
        append!(prolongations, out)
    end
    return prolongations
end

function substitute_lpa_data(equations, data)
    data_flat = reshape(data,1,:)
    [eqn(data_flat...) for eqn in equations]
end

nsteps=3
equations = prolongate_LPA(nsteps=nsteps, order=1)

# data generated in other file
include("LPA_simulations.jl")
data = R.(round.(LPAdata[:,1:nsteps+1]))

substitute_lpa_data(equations, data)
