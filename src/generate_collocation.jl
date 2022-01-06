using PCEUncertainty

function generate_collocation_points(N::Int)
    # Get the nodes and weights for the quadrature rule
    ortho_polynomials = GaussOrthoPoly(N)
    nodes = ortho_polynomials.quad.nodes
    weights = ortho_polynomials.quad.weights
    sols = []

    # Solve the ODE numerically at the collocation points over the probability space
    for i = 1:N
        # Set up the ODE
        f(u, p, t) = -p[1]*u
        u0 = 1.
        tspan = (0., 1.)
        p = nodes[i]

        # Solve the ODE
        prob = ODEProblem(f, u0, tspan, [p])
        sol = solve(prob)
        push!(sols, sol)
    end
    return nodes, weights, sols
end