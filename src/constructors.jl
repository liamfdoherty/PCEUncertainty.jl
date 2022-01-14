"""
`StochasticODEProblem` - initialize the data structure 
"""
function StochasticODEProblem(t_max, basis_type, basis_degree, collocation_size)
    basis = basis_type(basis_degree)
    collocation_nodes = basis_type(collocation_size).quad.nodes
    collocation_weights = basis_type(collocation_size).quad.weights
    return StochasticODEProblem(t_max, basis_type, basis_degree, basis, collocation_size, collocation_nodes, collocation_weights)
end