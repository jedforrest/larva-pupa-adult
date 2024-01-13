# larva-pupa-adult

We are using the Larva-Pupa-Adult (LPA) equation as a test case for designing
novel algorithms for solving parameter estimation problems for transcendental
systems which are amenable to polynomial approximation.

This repository documents this work.

## Index

* [lpa_params.py](/lpa_params.py) -- parameters as they are in the paper constantino_chaotic_dynamics_science_1997.pdf
* [simulation.py](/simulation.py) -- run the simulation using floating point numbers
* [homotopy_continuation.py](/homotopy_continuation.py) -- call HC.jl see Julia setup below
* [dynamic_taylor.py](/dynamic_taylor.py) -- run dynamic taylor algorithm

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

## Julia Setup

TODO(yberman)
