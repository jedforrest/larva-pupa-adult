# larva-pupa-adult

We are using the Larva-Pupa-Adult (LPA) equation as a test case for designing
novel algorithms for solving parameter estimation problems for transcendental
systems which are amenable to polynomial approximation.

This repository documents this work.

## Index

* [larva_pupa_adult_params.py](/larva_pupa_adult_params.py) -- parameters as they are in the paper constantino_chaotic_dynamics_science_1997.pdf
* [lpa_simulation.py](/lpa_simulation.py) -- run the simulation using floating point numbers
* [lpa_simulation.output.csv](/lpa_simulation.output.csv) -- output of lpa_simulation.py

## Python Setup
From the command line:

```bash
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```
## Julia Setup

```julia
julia> using Pkg;
julia> Pkg.add("DifferentialEquations")
```

TODO(yberman) What other requirents are there to run Josh's 

## Approximating exponential function


Given an interval $[a, b]$ we would like to approximate the exponential function. 

$$\min_{p} \max_{x \in [a, b]} |\exp(x) - p(x)|$$

* An issue we have encountered with Taylor expansion is over-constrained unsolvable systems. We would
like to filter the equations used according to their Jacobian rank.

### Jacobian testing

### Questions

* Why not just use [sgd](https://scikit-learn.org/stable/modules/sgd.html) on the parameters?
