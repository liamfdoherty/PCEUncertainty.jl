"""
`generate_pce_coefficients` - compute the PCE coefficients for the orthogonal representation of the interpolant
"""
function generate_pce_coefficients(prob::StochasticODEProblem, collocation_values)
    # Evaluate the orthogonal system of polynomials at the collocation nodes
    Φⱼz = evaluate(prob.collocation_nodes, prob.basis)
    
    # Normalize the orthogonal system
    for column in 1:size(Φⱼz)[2]
        Φⱼz[:, column] ./= sqrt(sum(Φⱼz[:, column] .* Φⱼz[:, column] .* prob.collocation_weights))
    end

    # Compute and return the PCE coefficients
    v̂ = Φⱼz\collocation_values
    return v̂
end