using Distributions
using Plots
using PCEUncertainty

N = 4

# Plot the true solution
x_vals = LinRange(-3, 3, 100)
truth(x) = exp(-x)
plt = plot(x_vals, truth.(x_vals), title = "Collocation points vs. Truth", label = "Truth")

# Compute and plot the collocation points
nodes, weights, sols = generate_collocation_points(N) # x values are Gauss-Hermite quadrature points
end_states = [sols[i].u[end] for i in 1:N]
scatter!(nodes, end_states, label = "Collocation Points")
xlabel!("z"); ylabel!("u(t = 1; z)")



display(plt)