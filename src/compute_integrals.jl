"""
`compute_integrals` - numerically evaluates the risk-sensitive form Λ and the two hybrid forms Λ¹ and Λ² for fixed values of c.

### Fields:
`quadrature_values` - values of the performance measure at the quadrature nodes
`quadrature_weights` - weights associated with the quadrature nodes
`c` - numerical parameter associated with the integrals
"""
function compute_integrals(quadrature_values, quadrature_weights, c::Real)
    @assert c > 0 "c must be greater than zero!"
    Λc = 0.
    for j = 1:size(quadrature_weights)[2]
        for i in 1:size(quadrature_weights)[1]
            Λc += exp(c*quadrature_values[i, j]) * quadrature_weights[i, j][1] * quadrature_weights[i, j][2]
        end
    end
    Λc = log(Λc)/c

    Λc¹ = 0.
    for j = 1:size(quadrature_weights)[2]
        exponent = 0.
        for i = 1:size(quadrature_weights)[1]
            exponent += c*quadrature_values[i, j]*quadrature_weights[i, j][1]
        end
        Λc¹ += exp(exponent)*quadrature_weights[1, j][2]
    end
    Λc¹ = log(Λc¹)/c

    Λc² = 0.
    for i = 1:size(quadrature_weights)[1]
        inner_sum = 0.
        for j = 1:size(quadrature_weights)[2]
            inner_sum += exp(c*quadrature_values[i, j])*quadrature_weights[i, j][2]
        end
        Λc² += log(inner_sum)*quadrature_weights[i, 1][1]
    end
    Λc² = Λc²/c

    return Λc, Λc¹, Λc²
end
