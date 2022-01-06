using PCEUncertainty

function generate_collocation_points(N::Int, )
    # Get the nodes and weights for the quadrature rule
    p_vals = sqrt(2)*gausshermite(N)[1]
    sols = []

    # Solve the ODE numerically at the collocation points over the probability space
    for i = 1:N
        # Set up the ODE
        f(u, p, t) = -p[1]*u
        u0 = 1.
        tspan = (0., 1.)
        p = p_vals[i]

        # Solve the ODE
        prob = ODEProblem(f, u0, tspan, [p])
        sol = solve(prob)
        push!(sols, sol)
    end
    return p_vals, sols
end