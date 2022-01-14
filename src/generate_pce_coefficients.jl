"""
`generate_pce_coefficients` - compute the PCE coefficients for the orthogonal representation of the interpolant
"""
function generate_pce_coefficients(prob::StochasticODEProblem, collocation_values)
    # Compute the v̂ⱼ (PCE) coefficients
    Φⱼz = evaluate(prob.collocation_nodes, prob.basis)
    v̂ = Φⱼz\collocation_values
    return v̂
end