using Distributions
using Plots
using PCEUncertainty

N = 5
dist = Uniform(-2, 2)

# Plot the true solution
x_vals = LinRange(-3, 3, 100)
truth(x) = exp(-x)
plt = plot(x_vals, truth.(x_vals), title = "Collocation points vs. Truth", label = "Truth")

# Compute and plot the collocation points

# p_vals, sols = generate_collocation_points(N, dist) # for if we want to sample collocation points from dist
p_vals, sols = generate_collocation_points(N) # x values are Gauss-Legendre quadrature points
end_states = [sols[i].u[end] for i in 1:N]
scatter!(p_vals, end_states, label = "Collocation Points")
xlabel!("z"); ylabel!("u(t = 1; z)")

# Compute and plot the interpolated solution



display(plt)