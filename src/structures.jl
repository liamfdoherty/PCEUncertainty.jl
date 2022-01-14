"""
`StochasticODEProblem` - structure to contain information about an ODE problem with random coefficient

Fields:
`t_max` - end time for the ODE simulation
`basis_type` - type of polynomials used in the PCE basis
`basis_degree` - maximum degree of basis polynomials
`basis` - data structure with the polynomial basis
`collocation_size` - number of collocation points to use in the probability space
`collocation_nodes` - set of nodes to perform the collocation with
`collocation weights` - weights associated with the collocation nodes
"""
struct StochasticODEProblem{TB}
    t_max::Real
    basis_type::TB
    basis_degree::Int
    basis
    collocation_size::Int
    collocation_nodes::Vector{Float64}
    collocation_weights::Vector{Float64}
end
