function compute_integrals(quadrature_values, quadrature_weights, c::Real)
    @assert c > 0 "c must be greater than zero!"
    Λc = 0.
    idx = 1
    for i = 1:length(quadrature_weights[1, :])
        for j = 1:length(quadrature_weights[2, :])
            Λc += exp(c*quadrature_values[idx]) * quadrature_weights[1, i] * quadrature_weights[2, j]
            idx += 1
        end
    end
    Λc = log(Λc)/c

    Λc¹ = 0.
    idx = 1
    for j = 1:length(quadrature_weights[2, :])
        exponent = 0.
        for i = 1:length(quadrature_weights[1, :])
            exponent += c*quadrature_values[idx]*quadrature_weights[1, i]
            idx += 1
        end
        Λc¹ += exp(exponent)*quadrature_weights[2, j]
    end
    Λc¹ = log(Λc¹)/c

    Λc² = 0.
    idx = 1
    for i = 1:length(quadrature_weights[1, :])
        argument = 0.
        for j = 1:length(quadrature_weights[2, :])
            argument += exp(c*quadrature_values[idx])*quadrature_weights[2, j]
            idx += 1
        end
        Λc² += log(argument)*quadrature_weights[1, i]
    end
    Λc² = Λc²/c

    return Λc, Λc¹, Λc²
end

function compute_integrals2(quadrature_values, quadrature_weights, c::Real)
    @assert c > 0 "c must be greater than zero!"
    Λc = 0.
    for idx = 1:size(quadrature_weights)[1]
        Λc += exp(c*quadrature_values[idx]) * quadrature_weights[idx, 1] * quadrature_weights[idx, 2]
    end
    Λc = log(Λc)/c

    return Λc
end