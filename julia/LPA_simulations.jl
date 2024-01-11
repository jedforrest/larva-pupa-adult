# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)
using OrdinaryDiffEq
using Random, Distributions

# Original LPA model (exp version)
function LPA(du, u, p, t)
    b, cel, cea, cpa, μl, μa = p
    L, P, A = u

    du[1] = b * A * exp(-cel * L - cea * A)
    du[2] = (1 - μl) * L
    du[3] = P * exp(-cpa * A) + (1 - μa) * A
end

function generate_synthetic_data(model, u0, params; steps=10)
    prob = DiscreteProblem(model, u0, (0., steps), params)
    sol = OrdinaryDiffEq.solve(prob, FunctionMap())
    sol[:,:]
end

# solving the original problem with provided parameters
# L0, P0, A0
u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

# Generates a matrix of values for LPA at t = 0, 1, ..., N
N = 40
LPAdata = generate_synthetic_data(LPA, u0, params; steps=N)
LPAdata[:,1] # L0, P0, A0
LPAdata[:,2] # L1, P1, A1
LPAdata[:,3] # L2, P2, A2
LPAdata[:,4] # L3, P3, A3


# sample new parameters and ICs
# parameter p is sampled from a range [p ± tol%]
sample_parameter(p; tol) = rand(Uniform(p*(1-tol), p*(1+tol)))
sampled_params = sample_parameter.(params; tol=0.3)

# sample ICs (u0) from a fixed integer distribution
sample_ICs(lb, ub) = Float64.(rand(DiscreteUniform(lb, ub)))
sampled_u0 = [sample_ICs(225, 275), sample_ICs(0, 10), sample_ICs(80, 120)]

# Generates a matrix of values for LPA at t = 0, 1, ..., N
N = 40
LPAdata = generate_synthetic_data(LPA, sampled_u0, sampled_params; steps=N)
LPAdata[:,1:4]
