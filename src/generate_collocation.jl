"""
`generate_collocation` - compute the solution to the ODE for a set of collocation points in the probability space
"""
function generate_collocation(prob::StochasticODEProblem)
    # Solve the ODE numerically at each of the collocation points over the probability space
    sols = []
    for i = 1:prob.collocation_size 
        # Set up the ODE
        f(u, p, t) = -p[1]*u
        u0 = 1.
        tspan = (0., prob.t_max)
        p = prob.collocation_nodes[i]

        # Solve the ODE
        ODEprob = ODEProblem(f, u0, tspan, [p])
        sol = solve(ODEprob)
        push!(sols, sol)
    end
    return sols
end