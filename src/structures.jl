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

    # Constructor for a single uncertainty
    function StochasticODEProblem(t_max::Real, basis_type, basis_degree::Int)
        basis = basis_type(basis_degree)
        collocation_nodes = basis_type(basis_degree + 1).quad.nodes
        collocation_weights = basis_type(basis_degree + 1).quad.weights
        return new(t_max, basis, collocation_nodes, collocation_weights)
    end

    # Constructor for multiple uncertainties
    function StochasticODEProblem(t_max::Real, basis_type::Vector, basis_degree::Vector{Int})
        @assert length(basis_type) == length(basis_degree) "Vector of types and vector of degrees do not match!"
        degree = minimum(basis_degree)

        basis_components = [basis_type[i](degree) for i ∈ 1:length(basis_type)]
        basis = MultiOrthoPoly(basis_components, minimum(basis_degree))

        higher_basis = [basis_type[i](degree + 1) for i ∈ 1:length(basis_type)]
        
        collocation_nodes_components = [higher_basis[i].quad.nodes for i ∈ 1:length(basis_components)]
        collocation_nodes = [collect(node) for node in vec(collect(Iterators.product(collocation_nodes_components...)))]
        collocation_nodes = permutedims(hcat(collocation_nodes...))

        collocation_weights_components = [higher_basis[i].quad.weights for i ∈ 1:length(basis_components)]
        collocation_weights = [collect(weight) for weight in vec(collect(Iterators.product(collocation_weights_components...)))]
        collocation_weights = permutedims(hcat(collocation_weights...))
        return new(t_max, basis, collocation_nodes, collocation_weights)
    end
end
