# LPA model with taylor approximated exp function
# Solving using HomotopyContinuation
using HomotopyContinuation
using TaylorSeries

# order of taylor series expansion
te_order = 1
# exp taylor approximation with rational coefficients
# exp_taylor = convert(Taylor1{Rational{Int}}, taylor_expand(exp, 0; order))
exp_taylor = taylor_expand(exp, 0; order=te_order)

function LPA_taylor(du, u, p, t; order=te_order)
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    exp_taylor = taylor_expand(exp, 0; order)

    du[1] = b * A * exp_taylor(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp_taylor(-cpa * A) + (1 - μa) * A
    return du
end


@var L P A p[1:6]
u0 = [L, P, A]
M = zeros(Expression,3,3)
for i in 1:3
    du_expr = zeros(Expression,3)
    exprs = LPA_taylor(du_expr, u0, p, 0)
    M[:,i] = exprs
    u0 = exprs
end

M

# Jacobian
F = System(M; variables=[L,P,A], parameters=p)
jacobian(F, [L,P,A], p)
