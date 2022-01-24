using Plots
using PolyChaos
using PCEUncertainty

# Set up the problem data structure
t_max = 1.
basis_type = GaussOrthoPoly

# Plot the interpolant for degrees 2,3,4,5
plots = []
for num_points in 2:5
    # Plot the true solution
    z_vals = LinRange(-3, 3, 100)
    truth(x) = exp(-x)
    subplot = plot(z_vals, truth.(z_vals), title = "N = $(num_points)", label = "Truth", color = :green)

    # Set up the problem with the given degree
    basis_degree = num_points - 1
    prob = StochasticODEProblem(t_max, basis_type, basis_degree)

    # Compute and plot the collocation points
    sols = generate_collocation(prob);
    end_states = [sols[i].u[end] for i in 1:length(prob.collocation_nodes)];
    scatter!(prob.collocation_nodes, end_states, label = "Collocation Points", color = :red)
    xlabel!("z"); ylabel!("u(t = 1; z)")

    # Compute and plot the interpolant
    v̂ = generate_pce_coefficients(prob, end_states)
    Φⱼz = evaluate(z_vals, prob.basis)
    for (n, column) in enumerate(1:size(Φⱼz)[2])
        Φⱼz[:, column] = Φⱼz[:, column] ./ sqrt(factorial(n - 1)) # This normalization is dependent on the underlying measure
    end
    v = zeros(length(Φⱼz[:, 1]))
    for i in 1:prob.basis.deg + 1
        v .+= v̂[i] .* Φⱼz[:, i]
    end
    plot!(z_vals, v, label = "Interpolated solution", color = :blue)

    push!(plots, subplot)
end

# Display all plots together
plt = plot(plots[1], plots[2], plots[3], plots[4], layout = 4)
display(plt)