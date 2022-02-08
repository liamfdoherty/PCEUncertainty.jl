"""
`generate_pce_coefficients` - compute the PCE coefficients for the orthogonal representation of the interpolant
"""
function generate_pce_coefficients(prob::StochasticODEProblem, collocation_values; vandermonde = false)
    # Evaluate the orthogonal system of polynomials at the collocation nodes
    inds, Φⱼz = evaluate_basis(prob.collocation_nodes, prob.basis)

    # Normalize the orthogonal system
    if length(prob.basis) > 1
        for column in 1:size(Φⱼz)[2]
            Φⱼz[:, column] ./= sqrt(sum(Φⱼz[:, column] .* Φⱼz[:, column] .* [w[1]*w[2] for w in prob.collocation_weights]))
        end
    else
        for column in 1:size(Φⱼz)[2]
            Φⱼz[:, column] ./= sqrt(sum(Φⱼz[:, column] .* Φⱼz[:, column] .* prob.collocation_weights))
        end
    end

    # Compute and return the PCE coefficients
    if vandermonde
        v̂ = Φⱼz\collocation_values
    else
        if length(prob.basis) > 1
            v̂ = zeros(size(Φⱼz)[1])
            for (j, v̂ⱼ) in enumerate(v̂)
                v̂[j] = sum(collocation_values .* Φⱼz[j, :] .* [w[1]*w[2] for w in prob.collocation_weights])
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