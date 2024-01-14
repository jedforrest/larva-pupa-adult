# LPA model with taylor approximated exp function
# Solving using HomotopyContinuation
using HomotopyContinuation
using TaylorSeries
using LinearAlgebra
# using AbstractAlgebra

# HC symbolic variables and parameters
nsteps = 3
HC_vars = @var L[0:nsteps] P[0:nsteps] A[0:nsteps]
HC_params = @var b cel cea cpa μl μa

#----------------------------------------------
function LPA_taylor!(out, unp1, un, p; exp_func=exp)
    b, cel, cea, cpa, μl, μa = p # parameters
    L1, P1, A1 = unp1 # u_n+1
    L, P, A = un # u_n

    out[1] = -L1 + b * A * exp_func(-cel * L - cea * A)
    out[2] = -P1 + (1 - μl) * L
    out[3] = -A1 + P * exp_func(-cpa * A) + (1 - μa) * A
    return nothing
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
function full_rank_subsystem(F::System)
    sys_vars = variables(F)
    R, ring_vars = polynomial_ring(QQ, string.(sys_vars).*"_qq")
    J = jacobian(F)
    m, n = size(J)

    # Convert Expression to PolynomialRing so that we can use rref and determine pivots
    J_r = matrix(R, [expr_to_ring(J[i,j], sys_vars, ring_vars) for i in 1:m, j in 1:n])

    # reduced row echelon form
    J_r = m > n ? transpose(J_r) : J_r
    _, A, d = rref_rational(J_r)
    pivot_cols = pivot_columns(A, d)

    # new subsystem
    System(expressions(F)[pivot_cols])
end

function expr_to_ring(expr, sys_vars, ring_vars)
    if iszero(expr)
        return 0
    else
        expons, coeffs = exponents_coefficients(expr, sys_vars)
        return create_polyn(ring_vars, expons, coeffs)
    end
end

function create_polyn(vars, expons, coeffs)
    # should work regardless of type of 'vars'
    sum(coeff * prod(vars .^ col) for (col, coeff) in zip(eachcol(expons), coeffs))
end

function pivot_columns(A, d)
    # pivot column if A[i,i]=1 and A[j,i] = 0, i!=j
    is_rref(A) || error("Matrix must be in reduced row echelon form")
    _A = Matrix(A) # must be a Julia Matrix type
    m, n = size(_A)
    dIm = d.*collect(I(m)) # Identity matrix * denominator
    pivot_cols = []
    for i in 1:n
        if _A[:,i] in eachcol(dIm)
            push!(pivot_cols, i)
        end
    end
    return pivot_cols
end
#----------------------------------------------
