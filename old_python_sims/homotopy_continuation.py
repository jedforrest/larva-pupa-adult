#!/usr/bin/env python3
# homotopy_continuation.py -- Python wrapper around www.juliahomotopycontinuation.org

from subprocess import run
from textwrap import dedent
import os

def to_julia_src(variables, system, julia_output_filename):
    variables_str = ' '.join(sorted(str(x) for x in variables))
    system_eqn = ', '.join(str(eqn) for eqn in system)
    return dedent(f'''
    using HomotopyContinuation
    @var {variables_str}
    f = System([{system_eqn}])
    result = solve(f)
    open("{julia_output_filename}", "w") do io
        println(io, real_solutions(result))
    end
    ''')


def process_output(variables, filename):
    solution_dicts = []
    with open(filename) as julia_output:
        file_content = julia_output.read()
        if file_content == 'Vector{Float64}[]\n':
            return []
        sol_lst = eval(file_content)
        for sol in sol_lst:
            d = {}
            for var, x in zip(variables, sol):
                d[var] = x
            solution_dicts.append(d)
    return solution_dicts

def solve(variables, system):
    julia_source_filename = 'hom_cont_input.jl'
    julia_output_filename = 'julia_output.txt'
    if os.path.exists(julia_source_filename):
        yn = input(f'The file `{julia_source_filename}` exists so the program cannot proceed. Delete? [y/n] ')
        assert yn == "y", f'please delete {julia_source_filename}'
        os.remove(julia_source_filename)
    if os.path.exists(julia_output_filename):
        yn = input(f'The file `{julia_output_filename}` exists so the program cannot proceed. Delete? [y/n] ')
        assert yn == "y", f'please delete {julia_output_filename}'
        os.remove(julia_output_filename)
    with open(julia_source_filename, 'w') as hc_in:
        hc_in.write(to_julia_src(variables, system, julia_output_filename))
    completed_process = run(['julia', julia_source_filename])
    assert completed_process.returncode == 0, f'julia failed. Try running input {julia_source_filename} directly'
    solutions = process_output(variables, julia_output_filename)
    os.remove(julia_source_filename)
    os.remove(julia_output_filename)
    return solutions


if __name__ == '__main__':
    solutions = solve(['x', 'y'], ['x^2 + y^2 - 4', 'x - y'])
    print(solutions)
