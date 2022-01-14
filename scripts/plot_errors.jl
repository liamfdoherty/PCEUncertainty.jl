using Plots
using PCEUncertainty

t_max = 1.
basis_type = GaussOrthoPoly
degrees = collect(2:11)

# Compute the means and variances
means = []; variances = []
for degree in degrees
    local basis_degree = degree; local collocation_size = degree
    local prob = StochasticODEProblem(t_max, basis_type, basis_degree, collocation_size)
    local sols = generate_collocation(prob)
    local end_states = [sols[i].u[end] for i in 1:prob.collocation_size]
    local v̂ = generate_pce_coefficients(prob, end_states)
    
    # Add means
    push!(means, v̂[1])

    # Add variances
    push!(variances, sum(v̂.^2) - v̂[1]^2)
end

# True mean and variance computed with Wolfram Alpha
true_mean = exp(1/2)
true_variance = exp(2) - exp(1)

mean_errors = [abs(mean - true_mean) for mean in means]
variance_errors = [abs(variance - true_variance) for variance in variances]
plot(degrees, mean_errors, label = "Mean error", title = "Relative errors", yaxis = :log, ylims = (1e-7, 1e0))
xlabel!("P"); ylabel!("error")
plot!(degrees, variance_errors, label = "Variance error")