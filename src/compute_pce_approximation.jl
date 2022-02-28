"""
`compute_pce_approximation` - compute the PCE approximation of the performance measure F at the quadrature nodes and the associated quadrature weights

### Fields:
`resolution` - degree of the basis to be used in the interpolation of the performance measure
`vandermonde` - whether or not to use inversion of the Vandermonde matrix for computation of PCE coefficients
`num_quadrature_nodes` - number of nodes to use for computing the risk-sensitive forms
"""
function compute_pce_approximation(resolution; vandermonde = false, num_quadrature_nodes = 255)
    # Set up the StochasticODEProblem
    t_max = 1.
    basis_type = [Uniform01OrthoPoly, Uniform01OrthoPoly];
    F(x) = x^2; basis_degree = [resolution, resolution]
    prob = StochasticODEProblem(t_max, basis_type, basis_degree)

    # Generate and plot the collocation points
    sols = generate_collocation_2d(prob);
    end_states = [sols[i].u[end] for i in 1:size(prob.collocation_nodes)[1]];

    # Compute the collocation solution for F
    end_states = F.(end_states)
    F̂ = generate_pce_coefficients(prob, end_states, vandermonde = vandermonde)

    # Compute the nodes and weights to be used for quadrature
    quadpoly = Uniform01OrthoPoly(num_quadrature_nodes)
    z₁_vals = quadpoly.quad.nodes; z₂_vals = quadpoly.quad.nodes
    z₁_weights = quadpoly.quad.weights; z₂_weights = quadpoly.quad.weights

    # Package the nodes and weights together
    z_vals_components = [z₁_vals, z₂_vals]
    z_vals = [collect(z_val) for z_val in vec(collect(Iterators.product(z_vals_components...)))]
    z_weights_components = [z₁_weights, z₂_weights]
    z_weights = [collect(z_weight) for z_weight in vec(collect(Iterators.product(z_weights_components...)))]
    z_weights = reshape(z_weights, (num_quadrature_nodes, num_quadrature_nodes))

    # Evaluate the basis and compute the interpolated performance measure at the quadrature nodes
    inds, Φⱼz = evaluate_basis(z_vals, prob.basis)
    for (j, col) in enumerate(1:size(Φⱼz)[2])
        degs = inds[j]
        Φⱼz[:, col] .*= binomial(2*degs[1], degs[1])*sqrt(2*degs[1] + 1)*binomial(2*degs[2], degs[2])*sqrt(2*degs[2] + 1) # This normalization is specific to U(0, 1)
    end
    F̃ = [dot(F̂, Φⱼz[row, :]) for row in 1:size(Φⱼz)[1]]
    F̃ = reshape(F̃, (num_quadrature_nodes, num_quadrature_nodes))
    return z_weights, F̃
end