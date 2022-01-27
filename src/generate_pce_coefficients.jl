"""
`generate_pce_coefficients` - compute the PCE coefficients for the orthogonal representation of the interpolant
"""
function generate_pce_coefficients(prob::StochasticODEProblem, collocation_values; vandermonde = false)
    # Evaluate the orthogonal system of polynomials at the collocation nodes
    Φⱼz = evaluate(prob.collocation_nodes, prob.basis)

    # Normalize the orthogonal system
    if typeof(prob.basis) <: MultiOrthoPoly
        for row in 1:size(Φⱼz)[1]
            Φⱼz[row, :] ./= sqrt(sum(Φⱼz[row, :] .* Φⱼz[row, :] .* prob.collocation_weights[:, 1] .* prob.collocation_weights[:, 2]))
        end
    else
        for column in 1:size(Φⱼz)[2]
            Φⱼz[:, column] ./= sqrt(sum(Φⱼz[:, column] .* Φⱼz[:, column] .* prob.collocation_weights))
        end
    end

    # Compute and return the PCE coefficients
    if vandermonde
        if typeof(prob.basis) <: MultiOrthoPoly
            v̂ = Φⱼz'\collocation_values # Artifact of the difference in how Φⱼz is stored for multivariate vs. univariate bases
        else
            v̂ = Φⱼz\collocation_values
        end
    else
        if typeof(prob.basis) <: MultiOrthoPoly
            v̂ = zeros(size(Φⱼz)[1])
            for (j, v̂ⱼ) in enumerate(v̂)
                v̂[j] = sum(collocation_values .* Φⱼz[j, :] .* prob.collocation_weights[:, 1] .* prob.collocation_weights[:, 2])
            end
        else
            v̂ = zeros(size(Φⱼz)[2])
            for (j, v̂ⱼ) in enumerate(v̂)
                v̂[j] = sum(collocation_values .* Φⱼz[:, j] .* prob.collocation_weights)
            end
        end
    end
    return v̂
end