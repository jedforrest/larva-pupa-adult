# larva-pupa-adult

We are using the Larva-Pupa-Adult (LPA) equation as a test case for designing
novel algorithms for solving parameter estimation problems for transcendental
systems which are amenable to polynomial approximation.

This repository is documenting this work.

## LPA Equations

Suppose you have three time series

* $`L = (l_1, l_2, \ldots)`$
* $`P = (p_1, p_2, \ldots)`$
* $`A = (a_1, a_2, \ldots)`$

representing the demographics of a sample

$L$ is the number of larva, $P$ is the

## Setup

This work is being done in Python. Install packages with:

```bash
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

## MSolve

[MSolve](https://msolve.lip6.fr/) is our tool of choice for solving `


### Buildling MSolve
```
git clone git@github.com:algebraic-solving/msolve.git
cd msolve/
./autogen.sh 
./configure 
make 
```

### Using MSolve

Msolve can be used to solve polynomial systems of equations.

## Approximating exponential function


Given an interval $[a, b]$ we would like to approximate the exponential function. 

$$\min_{p} \max_{x \in [a, b]} |\exp(x) - p(x)|$$

* An issue we have encountered with Taylor expansion is over-constrained unsolvable systems. We would
like to filter the equations used according to their Jacobian rank.

### Jacobian testing

### Questions

* Why not just use [sgd](https://scikit-learn.org/stable/modules/sgd.html) on the parameters?
