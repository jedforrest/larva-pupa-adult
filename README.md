# larva-pupa-adult

We are using the Larva-Pupa-Adult (LPA) equation as a test case for designing
novel algorithms for solving parameter estimation problems for transcendental
systems which are amenable to polynomial approximation.

This repository documents this work.

## Index

* [lpa_simulation.py](/lpa_simulation.py) builds the data in
  [lpa-data.txt](/lpa-data.txt) using the parameters defined in
  [larva_pupa_adult_params.py](/larva_pupa_adult_params.py).
* [larva_pupa_adult_params.py](/larva_pupa_adult_params.py) serves as "true" values of parameters based on quantities discovered at end of Constantino paper.

## LPA Equations

Suppose you have three time series

* $`L = (l_1, l_2, \ldots)`$
* $`P = (p_1, p_2, \ldots)`$
* $`A = (a_1, a_2, \ldots)`$

representing the demographics of a sample

TODO(yberman) Either finish or remove this description of the problem.

## Setup:

This work is being done in Python. Install required packages with:

```bash
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

## MSolve

[MSolve](https://msolve.lip6.fr/) is our tool of choice for solving systems of equations. To build MSolve:

```
git clone git@github.com:algebraic-solving/msolve.git
cd msolve/
./autogen.sh 
./configure 
make 
```

Msolve can be used to solve polynomial systems of equations. Consider $`x^2 + y^2 = 4`$ and $`x = y`$. See [msolve-input-example.txt](/msolve-input-example.txt). The command `msolve -f msolve-input-example.txt` outputs:

```
[0, [1,
[[[-481231938336009023090067544955250113855 / 2^128, -240615969168004511545033772477625056927 / 2^127], [-2219290601644883707169587795849264957011310523479721104423 / 2^190, -1109645300822441853584793897924632478505655261739860552211 / 2^189]], [[240615969168004511545033772477625056927 / 2^127, 481231938336009023090067544955250113855 / 2^128], [1109645300822441853584793897924632478505655261739860552211 / 2^189, 2219290601644883707169587795849264957011310523479721104423 / 2^190]]]
]]:
```

There are two solutions real solutions

* `[[-481231938336009023090067544955250113855 / 2^128, -240615969168004511545033772477625056927 / 2^127], [-2219290601644883707169587795849264957011310523479721104423 / 2^190, -1109645300822441853584793897924632478505655261739860552211 / 2^189]]`
* `[[240615969168004511545033772477625056927 / 2^127, 481231938336009023090067544955250113855 / 2^128], [1109645300822441853584793897924632478505655261739860552211 / 2^189, 2219290601644883707169587795849264957011310523479721104423 / 2^190]]`


## Approximating exponential function


Given an interval $[a, b]$ we would like to approximate the exponential function. 

$$\min_{p} \max_{x \in [a, b]} |\exp(x) - p(x)|$$

* An issue we have encountered with Taylor expansion is over-constrained unsolvable systems. We would
like to filter the equations used according to their Jacobian rank.

### Jacobian testing

### Questions

* Why not just use [sgd](https://scikit-learn.org/stable/modules/sgd.html) on the parameters?
