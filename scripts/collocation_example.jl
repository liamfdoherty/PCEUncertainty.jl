using Distributions
using Plots
using PolyChaos
using PCEUncertainty

# Set up the problem data structure
t_max = 1.
basis_type = GaussOrthoPoly; basis_degree = 5 # Check to see if basis_degree and collocation_size have to be the same
collocation_size = 5
prob = StochasticODEProblem(t_max, basis_type, basis_degree, collocation_size)

# Plot the true solution
z_vals = LinRange(-3, 3, 100)
truth(x) = exp(-x)
plt = plot(z_vals, truth.(z_vals), title = "Collocation points vs. Truth", label = "Truth")

# Compute and plot the collocation points
sols = generate_collocation(prob)
end_states = [sols[i].u[end] for i in 1:prob.collocation_size]
scatter!(prob.collocation_nodes, end_states, label = "Collocation Points")
xlabel!("z"); ylabel!("u(t = 1; z)")

# Compute and plot the interpolant
v̂ = generate_pce_coefficients(prob, end_states)
Φⱼz = evaluate(z_vals, prob.basis)
v = zeros(length(Φⱼz[:, 1]))
for i in 1:prob.basis_degree + 1
    v .+= v̂[i] .* Φⱼz[:, i]
end
plot!(z_vals, v, label = "Interpolated solution")

# Display all plots together
display(plt)