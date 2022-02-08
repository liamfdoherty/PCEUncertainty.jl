"""
`generate_collocation` - compute the solution to the ODE for a set of collocation points in the probability space
"""
function generate_collocation(prob::StochasticODEProblem)
    sols = []
    for i = 1:length(prob.collocation_nodes)
        # Set up the ODE
        f(u, p, t) = -p[1]*u
        u0 = 1.; tspan = (0., prob.t_max)
        p = prob.collocation_nodes[i]

        # Solve the ODE
        ODEprob = ODEProblem(f, u0, tspan, [p])
        sol = solve(ODEprob)
        push!(sols, sol)
    end
    return sols
end

"""
`generate_collocation2` - compute the solution to the ODE for a set of collocation points in the probability space; different from generate_collocation because this encodes
data about the underlying differential equation model.
"""
function generate_collocation2(prob::StochasticODEProblem)
    sols = []
    for i = 1:size(prob.collocation_nodes)[1]
        # Set up the ODE
        f(u, p, t) = -p[1]*u
        u0 = prob.collocation_nodes[i][2]; tspan = (0., prob.t_max)
        p = prob.collocation_nodes[i][1]

        # Solve the ODE
        ODEprob = ODEProblem(f, u0, tspan, [p])
        sol = solve(ODEprob)
        push!(sols, sol)
    end
    return sols
end