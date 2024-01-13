# "Larvae-Pupae-Adult" (LPA) model of the the flour beetle
# Model from 'Modeling Populations of Flour Beetles' (Robertson, 2009)
# Parameter values from 'Chaotic Dynamics in an Insect Population' (Constantino et al., 1997)

# Original parameters
# L0, P0, A0
original_u0 = [250., 5., 100.]
# b, cel, cea, cpa, μl, μa
original_params = [6.598, 1.209e-2, 1.155e-2, 4.7e-3, 0.2055, 7.629e-3]

# sample new parameters and ICs
# parameter p is sampled from a range [p ± tol%]
sample_parameter(p; tol=0.3) = rand(Uniform(p*(1-tol), p*(1+tol)))
sampled_params = sample_parameter.(original_params; tol=0.3)

# sample ICs (u0) from a fixed integer distribution
sample_ICs(ic; lb=ic-5, ub=ic+5) = Float64.(rand(DiscreteUniform(lb, ub)))
sampled_u0 = sample_ICs.(original_u0)
