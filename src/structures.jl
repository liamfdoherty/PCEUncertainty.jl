"""
`StochasticODEProblem` - structure to contain information about an ODE problem with random coefficient

Fields:
`t_max` - end time for the ODE simulation
`basis` - data structure with the polynomial basis
`collocation_nodes` - set of nodes to perform the collocation with
`collocation_weights` - set of weights for the collocation nodes
"""
struct StochasticODEProblem
    t_max::Real
    basis
    collocation_nodes
    collocation_weights

    function StochasticODEProblem(t_max, basis_type, basis_degree)
        basis = basis_type(basis_degree)
        collocation_nodes = basis_type(basis_degree + 1).quad.nodes
        collocation_weights = basis_type(basis_degree + 1).quad.weights
        return new(t_max, basis, collocation_nodes, collocation_weights)
    end
end
