import gurobipy as gp
from gurobipy import GRB


def solve_vaccine_problem(k, peptides, pep_count, alleles, freq_vector, B_matrix, min_hits=1, populations=['World'], path=''):
    # create new model
    m = gp.Model('naive_vaccine_design')

    # create peptide variables
    x = m.addVars(peptides, vtype=GRB.BINARY, name='peptide')

    # create hit variables
    hit = m.addVars(alleles, vtype=GRB.BINARY, name='allele_hit')

    # maximize population coverage
    m.setObjective(gp.quicksum(hit[allele] * freq_vector[allele][pop] for pop in populations for allele in alleles),
                   GRB.MAXIMIZE)

    # add constraints
    m.addConstrs(x[i] <= 1 for i in peptides)

    m.addConstr(gp.quicksum(x[i] * len(i) for i in peptides) <= k, 'Sum of peptide lengths smaller k')

    m.addConstrs(hit[a] * min_hits <= gp.quicksum(x[i] * B_matrix[a][i] for i in peptides) for a in alleles)

    m.optimize()
    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_concat.lp')
    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_concat.sol')
    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_concat.json')

    chosen_pep = []
    if m.Status == GRB.OPTIMAL:
        sol = m.getAttr('X', x)
        for peptide in x:
            print('%s: %g' % (peptide, sol[peptide]))
            if sol[peptide]:
                chosen_pep.append(peptide)
    
    chosen_pep = set(chosen_pep)
    with open(path + 'pep_out/' + pep_count + '_chosen_peptides_concat.txt', 'w') as file:
        file.write('\n'.join(chosen_pep))
