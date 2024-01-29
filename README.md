# larva-pupa-adult

We are using the Larva-Pupa-Adult (LPA) equation as a test case for designing
novel algorithms for solving parameter estimation problems for transcendental
systems which are amenable to polynomial approximation.

This repository documents this work.

## Index

* [python_sim_fixed_center.jl](/julia/python_sim_fixed_center.jl) -- main code to run in Julia that carries out testing of the quality of parameter estimation as described in Section 4 of the paper "Symbolic-numeric algorithm for parameter estimation in discrete-time
 models with exp" by Yosef Berman, Joshua Forrest, Matthew Grote, Alexey Ovchinnikov, and Sonia Rueda
* [simulation_results_fixed_center_24-01-24T21.csv](/tables/simulation_results_fixed_center_24-01-24T21.csv) -- table of randomized parameter estimation experiments as described in the above paper. Such a file is produced by the above Julia code
* [results_table_means.jl](/julia/results_table_means.jl) -- code to run in Julia to output mean and median errors based on the "simulation_results_fixed_center_xxx.csv" file
* [error_table_means_medians.csv](/tables/error_table_means_medians.csv) -- table with errors as described in the above paper produced by the above Julia code
* [lpa_params.py](/old_python_sims/lpa_params.py) -- parameters as they are in the paper "Chaotic Dynamics in an Insect Population" by Costantino et al
* [simulation.py](/old_python_sims/simulation.py) -- run the simulation using floating point numbers
* [homotopy_continuation.py](/old_python_sims/homotopy_continuation.py) -- call HC.jl see Julia setup below
* [dynamic_taylor.py](/old_python_sims/dynamic_taylor.py) -- run dynamic taylor algorithm

## Python Setup
From the command line:

```bash
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```
## Julia Setup

To install all dependencies from a Julia terminal run:
```julia-repl
julia> using Pkg
julia> Pkg.instantiate()
```
This installs all missing packages from `Project.toml`. Alternatively, from the Pkg REPL run:
```julia-repl
(larva-pupa-adult) pkg> instantiate
```

## Over-determined systems solving

This approach was based on:
[1] Bender, M.R. and Telen, S., 2021. Yet another eigenvalue algorithm for solving polynomial systems. arXiv preprint arXiv:2105.08472.
