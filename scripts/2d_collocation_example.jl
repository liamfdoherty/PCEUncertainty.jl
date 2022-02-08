using Plots
using LinearAlgebra
using PolyChaos
using PCEUncertainty

# Set up the StochasticODEProblem
t_max = 1.
basis_type = [Uniform01OrthoPoly, Uniform01OrthoPoly];
# F(x) = x^2; basis_degree = [7, 7]
F(x) = 0.8 ≤ x ≤ 1 ? 1 : 0; basis_degree = [11, 11]
prob = StochasticODEProblem(t_max, basis_type, basis_degree)

# Plot the true solution of the ODE
x_vals = LinRange(0, 1, 100); y_vals = LinRange(0, 1, 100)
truth(z₁, z₂) = z₂*exp(-z₁)
plt = surface(x_vals, y_vals, truth)

# Generate and plot the collocation points
sols = generate_collocation2(prob);
end_states = [sols[i].u[end] for i in 1:size(prob.collocation_nodes)[1]];
x_nodes = [node[1] for node in prob.collocation_nodes]; y_nodes = [node[2] for node in prob.collocation_nodes]
scatter!(collect(zip(x_nodes, y_nodes, end_states)), label = "Collocation Points")
display(plt)

# Compute the collocation solution for F
end_states = F.(end_states)
F̂ = generate_pce_coefficients(prob, end_states, vandermonde = false)

quadpoly = Uniform01OrthoPoly(2^8 - 1)
z₁_vals = quadpoly.quad.nodes; z₂_vals = quadpoly.quad.nodes
z₁_weights = quadpoly.quad.weights; z₂_weights = quadpoly.quad.weights

z_vals_components = [z₁_vals, z₂_vals]
z_vals = [collect(z_val) for z_val in vec(collect(Iterators.product(z_vals_components...)))]
# z_vals = permutedims(hcat(z_vals...))
z_weights_components = [z₁_weights, z₂_weights]
z_weights = [collect(z_weight) for z_weight in vec(collect(Iterators.product(z_weights_components...)))]
# z_weights = permutedims(hcat(z_weights_components...))

inds, Φⱼz = evaluate_basis(z_vals, prob.basis) # Need to normalize this
for (j, col) in enumerate(1:size(Φⱼz)[2])
    degs = inds[j] # degrees of polynomials in the tensor product for the jth basis element
    Φⱼz[:, col] .*= binomial(2*degs[1], degs[1])*sqrt(2*degs[1] + 1)*binomial(2*degs[2], degs[2])*sqrt(2*degs[2] + 1) # This normalization is specific to U(0, 1)
end
F̃ = [dot(F̂, Φⱼz[row, :]) for row in 1:size(Φⱼz)[1]]

# Compute the risk-sensitive integrals
c_vals = 0.01:1.:60.01
Λc_vals = []; Λc¹_vals = []; Λc²_vals = []
for c in c_vals
    Λc, Λc¹, Λc² = compute_integrals(F̃, z_weights, c)
    push!(Λc_vals, Λc)
    # push!(Λc¹_vals, Λc¹)
    # push!(Λc²_vals, Λc²)
end

plt2 = plot(c_vals, Λc_vals, title = "Λc values", xlabel = "c", label = "Λc", color = :red, legend = :outertopright)
# plot!(c_vals, Λc¹_vals, label = "Λc¹", color = :green)
# plot!(c_vals, Λc²_vals, label = "Λc²", color = :blue)
ylims!(-0.5, 1.); yticks!(0.:0.1:1.)
display(plt2)