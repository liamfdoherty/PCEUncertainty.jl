"""
`evaluate_basis` - evaluate the polynomial basis (univariate or multivariate) at a vector of nodes and return a matrix with this information

### Fields:
`nodes` - set of nodes to evaluate the basis at
`basis` - vector of (univariate) basis functions. If length(basis) == 1, then the basis is univariate. Otherwise, the basis is a multivariate basis
constructed as a tensor product of univariate bases.
"""
function evaluate_basis(nodes, basis)
    @assert all([length(basis) == length(nodes[i]) for i in 1:length(nodes)]) # check that there are the right number of node components for the basis
    
    # evaluate all of the basis components at the nodes
    basis_components_evaluations = []
    for (i, basis_component) in enumerate(basis)
        nodes_i = [nodes[j][i] for j in 1:length(nodes)] # get the right components of the nodes for the basis and evaluate
        evaluation = evaluate(nodes_i, basis_component)
        push!(basis_components_evaluations, evaluation)
    end

    # get all the valid combinations of basis degree for the tensor product of univariate polynomials
    degrees = [basis_component.deg for basis_component in basis]
    degree_iterators = [0:degrees[i] for i in 1:length(degrees)]
    inds = vec(collect(Iterators.product(degree_iterators...)))
    
    # construct the columns of the Φ matrix, each corresponding to a basis function
    Φ = zeros(length(nodes), length(inds))
    for (i, basis_element_degrees) in enumerate(inds)
        basis_element_at_nodes = ones(length(nodes))
        for (j, variable_degree) in enumerate(basis_element_degrees)
            basis_element_at_nodes .*= basis_components_evaluations[j][:, variable_degree + 1]
        end
        Φ[:, i] = basis_element_at_nodes
    end
    
    return inds, Φ
end