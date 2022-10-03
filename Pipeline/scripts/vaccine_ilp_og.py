import draw_output
import gurobipy as gp
from gurobipy import GRB


def log(string):
    if logging_enabled:
        print(string)


# subtour elimination version 1
def subtour_elim_cycle(model, where):
    global save_for_later
    if where == GRB.Callback.MIPSOL:
        vals = model.cbGetSolution(model._x_edges)
        # find cycle
        subtour_edges = find_subtour(vals)
        # add sub tour elimination constr. for specific sub tour
        if subtour_edges:
            model.cbLazy(gp.quicksum(model._x_edges[i, j] for i, j in subtour_edges) <= len(subtour_edges) - 1)
            save_for_later.append(subtour_edges)


# subtour elimination version 2
def subtour_elim_boekler_paper(model, where):
    global save_for_later
    if where == GRB.Callback.MIPSOL:
        vals = model.cbGetSolution(model._x_edges)
        # find cycle
        subtour_edges = find_subtour(vals)
        # sum over all incoming edges of node v in W <= sum of all incoming edges to node set W
        if subtour_edges:
            subtour_nodes = [node for node, _ in subtour_edges]
            save_for_later.append(subtour_nodes)
            for node in subtour_nodes:
                model.cbLazy(gp.quicksum(model._x_edges.select('*', node)) <=
                             gp.quicksum(model._x_edges[node_outer, node_inner] for node_outer, node_inner in model._x_edges
                                         if node_inner in subtour_nodes and node_outer not in subtour_nodes))
                # first add constraint only for single node, if this does not work, add it for all nodes in cycle
                break


def find_subtour(vals):
    # chosen edges
    edges = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j])
    nodes = ['Source'] + [node for node, _ in edges if node != 'Source']
    while nodes:
        current, nextn = edges.select(nodes[0], '*')[0]
        nodes.remove(current)
        curr_tour_edges = [(current, nextn)]
        while nextn != 'Sink':
            if nextn not in nodes:
                return curr_tour_edges
            else:
                nodes.remove(nextn)
            current, nextn = edges.select(nextn, '*')[0]
            curr_tour_edges.append((current, nextn))
    return None


def solve_vaccine_problem(k, alleles, freq_vector, B_matrix, peptides, pep_count, graph, path, min_hits=1, populations=['World'], subtour_el='cycle', logging=False, coloring=False):
    global logging_enabled
    logging_enabled = logging

    # create new model
    m = gp.Model('vaccine_overlap_design')
    # m.setParam('TimeLimit', 36000)

    # create allele hit variables
    log('Create allele hit variables')
    hit = gp.tupledict()
    for allele in alleles:
        hit[allele] = m.addVar(vtype=GRB.BINARY, name='_'.join(allele).replace(':', '-'))

    # create edge variables
    log('Create edge variables')
    x_edges = gp.tupledict()
    for node_a, node_b in graph.edges():
        x_edges[node_a, node_b] = m.addVar(vtype=GRB.BINARY, name=node_a + '_' + node_b)

    # the objective function is to maximize population coverage
    log('Set objective function')
    m.setObjective(gp.quicksum(hit[allele] * freq_vector[allele][pop] for pop in populations for allele in alleles),
                   GRB.MAXIMIZE)

    # subject to the following constraints
    log('Add constraints')
    m.addConstr(x_edges.sum('Source', '*') == 1, 'Exactly one outgoing edge from the source')
    m.addConstr(x_edges.sum('*', 'Sink') == 1, 'Exactly one incoming edge to the sink')

    log('Add >>At most one outgoing edge constraint<<')
    for node_a in graph.nodes():
        if node_a == 'Source' or node_a == 'Sink':
            continue
        m.addConstr(x_edges.sum(node_a, '*') <= 1, 'At most one outgoing edge from ' + node_a)
        '''
        Due to flow conservation, we might not need constraint of
        'at most one incoming edge per node' -> skip it here;
        '''
        m.addConstr(x_edges.sum('*', node_a) == x_edges.sum(node_a, '*'), 'Flow conservation at node ' + node_a)

    log('Sum of edge weights below threshold k constraint')
    m.addConstr(gp.quicksum(x_edges[node_a, node_b] * graph.edges[node_a, node_b]['weight']
                            for node_a, node_b in graph.edges()) <= k, 'Sum of edge weights below threshold k')

    '''
    without this update function, the data will be cached in some intermediate storage that might require more memory
    than the final storage to which the ilp will be transferred by update() and optimize().
    '''
    m.update()
    log('Allele hit constraint')
    for allele in alleles:
        # check for those peptides with non-zero entry in B matrix for current allele
        hit_nodes = list(filter(lambda x: 1.0 == B_matrix[allele][x], peptides))
        # only sum over such peptides instead of summing over all peptides -> only beneficial for ba threshold > 0.xx
        m.addConstr(hit[allele] * min_hits <= gp.quicksum(x_edges.sum(node_a, '*') for node_a in hit_nodes),
                    'Peptide-hit for allele ' + '_'.join(allele))
        m.update()

    # sub tour elimination with callbacks
    m._x_edges = x_edges
    m.Params.LazyConstraints = 1
    # keep subtour to add them after solving to lp file: lazy constraints are not added to lp file otherwise
    global save_for_later
    save_for_later = []
    chosen_pep = []
    # Solve
    log('Now optimize:')
    if subtour_el == 'cycle':
        m.optimize(subtour_elim_cycle)
    else:
        m.optimize(subtour_elim_boekler_paper)
    if m.Status == GRB.OPTIMAL:
        sol = m.getAttr('X', x_edges)
        for node_a, node_b in x_edges:
            # print('%s -> %s: %g' % (node_a, node_b, sol[node_a, node_b]))
            if sol[node_a, node_b]:
                chosen_pep.append(node_a)
    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_og.sol')
    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_og.json')

    if coloring:
        draw_output.draw_colored_graph(graph, sol, path + 'pep_out/' + pep_count + '_colored_og.png')

    # write chosen peptides to output file
    with open(path + 'pep_out/' + pep_count + '_chosen_peptides_og.txt', 'w') as file:
        chosen_pep.remove('Source')
        file.write('\n'.join(chosen_pep))

    # after solving add lazy constraints to lp file
    for i, subtour in enumerate(save_for_later):
        if subtour_el == 'cycle':
            m.addConstr(gp.quicksum(m._x_edges[i, j] for i, j in subtour) <= len(subtour) - 1,
                        'Subtour elimination ' + str(i))
        else:
            for node in subtour:
                m.addConstr(gp.quicksum(m._x_edges.select('*', node)) <=
                             gp.quicksum(m._x_edges[node_outer, node_inner] for node_outer, node_inner in m._x_edges
                                         if node_inner in subtour and node_outer not in subtour),
                            'Subtour elimination ' + node)
                # first add constraint only for single node, if this does not work, add it for all nodes in cycle
                break

    m.write(path + 'lp_out/' + pep_count + '_vaccine_ilp_og.lp')
