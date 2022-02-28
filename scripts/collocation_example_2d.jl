using Plots
using QuadGK
using PCEUncertainty

# Compute the interpolated performance measure and associated weights for quadrature
z_weights, F̃ = compute_pce_approximation(7)

# Compute the risk-sensitive integrals numerically
c_vals = 0.01:1.:60.01
Λc_vals = []; Λc¹_vals = []; Λc²_vals = []
for c in c_vals
    Λc, Λc¹, Λc² = compute_integrals(F̃, z_weights, c)
    push!(Λc_vals, Λc)
    push!(Λc¹_vals, Λc¹)
    push!(Λc²_vals, Λc²)
end

# Compute the true risk-sensitive integrals with quadrature
true_Λc = []; true_Λc¹ = []; true_Λc² = []
for c in c_vals
    Λc = log(quadgk(z₁ -> quadgk(z₂ -> exp(c*z₂^2*exp(-2*z₁)), 0., 1.)[1], 0., 1.)[1])/c
    push!(true_Λc, Λc)

    Λc¹ = log(quadgk(z₂ -> exp(quadgk(z₁ -> c*z₂^2*exp(-2*z₁), 0., 1.)[1]), 0., 1.)[1])/c
    push!(true_Λc¹, Λc¹)

    Λc² = quadgk(z₁ -> log(quadgk(z₂ -> exp(c*z₂^2*exp(-2*z₁)), 0., 1.)[1]), 0., 1.)[1]/c
    push!(true_Λc², Λc²)
end

# Plot the numerically computed risk-sensitive integrals
plt = plot(c_vals, Λc_vals, title = "Λc values", xlabel = "c", label = "Λc", color = :red, legend = :outertopright)
plot!(c_vals, Λc¹_vals, label = "Λc¹", color = :purple)
plot!(c_vals, Λc²_vals, label = "Λc²", color = :orange)

# Plot the true risk-sensitive integrals computed with a quadrature package
plot!(c_vals, true_Λc, label = "True Λc", linestyle = :dash, color = :green)
plot!(c_vals, true_Λc¹, label = "True Λc¹", linestyle = :dash, color = :yellow)
plot!(c_vals, true_Λc², label = "True Λc²", linestyle = :dash, color = :blue)

ylims!(0., 1.); yticks!(0.:0.1:1.)
display(plt)
