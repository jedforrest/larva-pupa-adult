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

k = 2 # order
t = Taylor1(k)
exp(t)

function lagrange_error_bound(f, c; a=0, tol=1e-3, max_k=19)
    # for function f(x) x∈[a,c] determine order k such that erro |ε| < tol
    # See: https://en.wikipedia.org/wiki/Taylor%27s_theorem#Estimates_for_the_remainder
    M = ceil(exp(c))
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
