"""
`generate_pce_coefficients` - compute the PCE coefficients for the orthogonal representation of the interpolant
"""
function generate_pce_coefficients(prob::StochasticODEProblem, collocation_values)
    # Compute the v̂ⱼ (PCE) coefficients
    Φⱼz = evaluate(prob.collocation_nodes, prob.basis) # unnormalized (i.e., not accounting for n!)
    for (n, column) in enumerate(1:size(Φⱼz)[2])
        Φⱼz[:, column] = Φⱼz[:, column] ./ sqrt(factorial(n - 1))
    end
    v̂ = Φⱼz\collocation_values
    return v̂
end