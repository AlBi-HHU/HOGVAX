import os
import numpy
import random
import draw_hog
import gurobipy as gp
import networkx as nx
from gurobipy import GRB


def log(string):
    if logging_enabled:
        print(string)


# subtour elimination version 2
def subtour_elim_boekler_paper(model, where):
    global save_for_later
    if where == GRB.Callback.MIPSOL:
        vals = model.cbGetSolution(model._x_edges)
        graph = model._graph
        # find cycle
        subtours = find_subtour(vals, graph)
        # sum over all incoming edges of node v in W <= sum of all incoming edges to node set W
        for subtour in sorted(subtours):
            subtour = sorted(list(subtour), key=int)
            log('## Subtour ' + ', '.join(subtour))
            save_for_later.append(subtour)
            model.cbLazy(1 <= gp.quicksum(model._x_edges[node_outer, node_inner] for node_outer, node_inner in sorted(model._x_edges)
                                          if node_inner in subtour and node_outer not in subtour),
                         'At least one incoming edge to subtour ' + '_'.join(subtour) + ' from the outside')


def find_subtour(vals, graph):
    # chosen edges
    edges = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] >= 0.5)
    log('Subtours edges')
    log(edges)
    subgraph = graph.edge_subgraph(edges).copy()
    # finding strongly connected components with networkx contains a random factor: algorithm starts from random node
    # generate sorted list of strongly connected components to remove random factor of scc algorithm
    scc = sorted(nx.strongly_connected_components(subgraph), key=len, reverse=True)
    subtours = []
    for c in scc:
        if '0' not in c:
            subtours.append(c)
    return subtours


def hogvax(k, alleles, freq_vector, B_matrix, leaves, pep_count, graph, path, min_hits=1, populations=['World'], optional_peptides=[], maximize_peptides=False, logging=False, coloring=False):
    global logging_enabled
    logging_enabled = logging

    # create new model
    m = gp.Model('ivp_on_hog')
    # set time limit to 48 hours
    m.setParam('TimeLimit', 172800)
    m.setParam('Seed', 42)
    random.seed(42)
    numpy.random.seed(42)
    # only for debugging, lp files becomes very large
    # if not os.path.exists(path + 'lp_out/'):
    #     os.mkdir(path + 'lp_out/')

    # create hit variables
    log('Create hit variables')
    hit = gp.tupledict()
    for allele in alleles:
        hit[allele] = m.addVar(vtype=GRB.BINARY, name='_'.join(allele).replace(':', '-'))

    # create node variables
    log('Create node variables for peptide nodes')
    x_nodes = gp.tupledict()
    for node in leaves.values():
        x_nodes[node] = m.addVar(vtype=GRB.BINARY, name=graph.nodes[node]['string'])

    # create edge variables
    log('Create edge variables')
    x_edges = gp.tupledict()
    for node_a, node_b in sorted(graph.edges):
        x_edges[node_a, node_b] = m.addVar(vtype=GRB.INTEGER,
                                           name=graph.nodes[node_a]['string'] + '_' + graph.nodes[node_b]['string'])

    # the objective function is to maximize population coverage
    log('Set objective function')
    m.setObjective(gp.quicksum(hit[allele] * freq_vector[allele][pop] for pop in populations for allele in alleles),
                   GRB.MAXIMIZE)

    # subject to the following constraints
    log('Add >>At least one outgoing edge from root<< constraint')
    m.addConstr(x_edges.sum('0', '*') >= 1, 'At least one outgoing edge from root')

    # optionally force include peptides in vaccine
    if optional_peptides:
        for p in optional_peptides:
            n = leaves[p]
            m.addConstr(x_nodes[n] >= 1, 'Must include peptide ' + p)

    log('Add >>Sum of edge weights below threshold k<< constraint')
    m.addConstr(gp.quicksum(x_edges[node_a, node_b] * graph.edges[node_a, node_b]['weight']
                            for node_a, node_b in sorted(graph.edges)) <= k, 'Sum of edge weights below threshold k')

    log('Add flow conservation for nodes')
    for node in sorted(graph.nodes):
        m.addConstr(x_edges.sum('*', node) == x_edges.sum(node, '*'),
                    'Flow conservation at node ' + str(node))

    log('Add peptide node visiting constraint')
    for pep, node in leaves.items():
        m.addConstr(x_nodes[node] <= x_edges.sum('*', node), 'Node visiting constraint for ' + pep)

    '''
    without this update function, the data will be cached in some intermediate storage that might require more memory
    than the final storage to which the ilp will be transferred by update() and optimize().
    '''
    m.update()
    log('Add allele hit constraint')
    for allele in alleles:
        # check for peptides (leaf nodes) with non-zero entry in B matrix for current allele
        hit_leaves = list(filter(lambda x: 1.0 == B_matrix[allele][x], leaves.keys()))
        # only sum over such peptides instead of summing over all peptides -> only beneficial for ba threshold > 0.xx
        m.addConstr(hit[allele] * min_hits <= gp.quicksum(x_edges.sum('*', leaves[leaf]) for leaf in hit_leaves),
                    'Peptide-hit for allele ' + '_'.join(allele))
    m.update()

    # sub tour elimination with callbacks
    m._x_edges = x_edges
    m._graph = graph
    m.Params.LazyConstraints = 1
    # keep subtour to add them after solving to lp file: lazy constraints are not added to lp file otherwise
    global save_for_later
    save_for_later = []
    chosen_pep = []
    # Solve
    log('Now optimize:')
    m.optimize(subtour_elim_boekler_paper)

    # 2nd optimization: maximize number of peptides
    if maximize_peptides:
        # require optimal solution for first objective function
        assert m.Status == GRB.Status.OPTIMAL
        opt_coverage = m.getObjective()
        # set new objective function
        # maximize number of selected peptides
        m.setObjective(gp.quicksum(x_nodes), GRB.MAXIMIZE)

        # add previous objective function with objective value as constraint
        m.addConstr(gp.quicksum(hit[allele] * freq_vector[allele][pop] for pop in populations for allele in alleles) >= opt_coverage.getValue(),
                    'Add maximal population coverage as constraint')
        m.optimize(subtour_elim_boekler_paper)

    solution_edges = []
    if m.Status == GRB.OPTIMAL:
        sol = m.getAttr('X', x_edges)
        for edge in x_edges:
            node_a, node_b = edge
            # print('%s -> %s: %g' % (node_a, node_b, sol[node_a, node_b]))
            # variable tolerance of gurobi is 1e-5
            traversals = round(sol[node_a, node_b])
            for i in range(traversals):
                solution_edges.append(edge)
                if node_a in leaves.values():
                    chosen_pep.append(graph.nodes[node_a]['string'])
                if node_b in leaves.values():
                    chosen_pep.append(graph.nodes[node_b]['string'])

    if not os.path.exists(path + 'lp_out/'):
        os.mkdir(path + 'lp_out/')
    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_hog.sol')
    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_hog.json')

    if not os.path.exists(path + 'pep_out/'):
        os.mkdir(path + 'pep_out/')

    if coloring:
        draw_hog.draw_hog(graph, leaves, sol, path + 'pep_out/' + pep_count + '_colored_hog.png', True)

    # write chosen peptides to output file, remove empty string from root node
    if '' in chosen_pep:
        chosen_pep.remove('')
    chosen_pep = list(set(chosen_pep))
    with open(path + 'pep_out/' + pep_count + '_chosen_peptides_hog.txt', 'w') as file:
        file.write('\n'.join(chosen_pep))

    # create subgraph from edges chosen by ILP in order to find eulerian path to construct vaccine sequence
    # use multi graph here to add edges that are traversed more than once multiple times to the graph
    log('Create vaccine sequence')
    sub_hog = nx.MultiDiGraph()
    sub_hog.add_edges_from(solution_edges)
    log('sub graph is eulerian' + str(nx.is_eulerian(sub_hog)))
    if coloring:
        dot = nx.drawing.nx_pydot.to_pydot(sub_hog)
        dot.write_png(path + 'pep_out/' + pep_count + '_sub_hog.png')
        draw_hog.draw_hog(graph, leaves, sol, path + 'pep_out/' + pep_count + '_colored_sub_hog.png')

    hogvaxine = ''
    for e in nx.eulerian_circuit(sub_hog, source='0'):
        if graph.edges[e]['is_slink'] == 'yes':
            continue
        else:
            node_a, node_b = e
            length = graph.edges[e]['weight']
            node_seq = graph.nodes[node_b]['string']
            hogvaxine += node_seq[len(node_seq) - length:]

    with open(path + 'pep_out/' + pep_count + '_hogvaxine.txt', 'w') as out:
        out.write('> MHC optimized combined peptide vaccine sequence with overlaps\n')
        out.write(hogvaxine + '\n')
        out.write('> MHC optimized combined peptide vaccine sequence concatenated\n')
        out.write(''.join(chosen_pep))

    # after solving add lazy constraints to lp file
    for i, subtour in enumerate(save_for_later):
        m.addConstr(1 <= gp.quicksum(x_edges[node_outer, node_inner] for node_outer, node_inner in sorted(x_edges)
                                     if node_inner in subtour and node_outer not in subtour),
                    'At least one incoming edge to subtour ' + str(i) + ' from the outside')

    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_hog.lp')

    return chosen_pep
