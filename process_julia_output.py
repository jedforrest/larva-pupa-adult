import larva_pupa_adult_params as PARAMS
solutions = eval(open('julia_output.txt').read())
varnames = 'a_0 a_1 a_2 a_3 b c_ea c_el c_pa l_0 l_1 l_2 l_3 mu_a p_0 p_1 p_2 p_3'.split()
for i, solution in enumerate(solutions):
    d = {}
    for variable, value in zip(varnames, solution):
        d[variable] = value
    for c in 'alp':
        for j in range(4):
            del d[f'{c}_{j}']
    print(i, d)

d = {
        'b': PARAMS.b,
        'c_ea': PARAMS.c_ea,
        'c_el': PARAMS.c_el,
        'c_pa': PARAMS.c_pa,
        'mu_a': PARAMS.mu_a,
}
print(d)
