using Plots
using PCEUncertainty

# Compute the true risk-sensitive integrals with quadrature
c_vals = [30., 60.]
true_Λc = []; true_Λc¹ = []; true_Λc² = []
for c in c_vals
    Λc = log(quadgk(z₁ -> quadgk(z₂ -> exp(c*z₂^2*exp(-2*z₁)), 0., 1.)[1], 0., 1.)[1])/c
    push!(true_Λc, Λc)

    Λc¹ = log(quadgk(z₂ -> exp(quadgk(z₁ -> c*z₂^2*exp(-2*z₁), 0., 1.)[1]), 0., 1.)[1])/c
    push!(true_Λc¹, Λc¹)

    Λc² = quadgk(z₁ -> log(quadgk(z₂ -> exp(c*z₂^2*exp(-2*z₁)), 0., 1.)[1]), 0., 1.)[1]/c
    push!(true_Λc², Λc²)
end

relative_errors_30 = []; relative_errors_60 = []; resolutions = collect(3:15)
for N in resolutions
    # Compute the interpolated performance measure and associated weights for quadrature
    local z_weights, F̃ = compute_pce_approximation(N)

    # Compute the risk-sensitive integrals at c = 30 and c = 60
    Λc_30, Λc¹_30, Λc²_30 = compute_integrals(F̃, z_weights, 30.)
    Λc_60, Λc¹_60, Λc²_60 = compute_integrals(F̃, z_weights, 60.)

    # Compute relative errors
    Λc_error_30 = abs(Λc_30 - true_Λc[1])/true_Λc[1]; Λc_error_60 = abs(Λc_60 - true_Λc[2])/true_Λc[2]
    Λc¹_error_30 = abs(Λc¹_30 - true_Λc¹[1])/true_Λc¹[1]; Λc¹_error_60 = abs(Λc¹_60 - true_Λc¹[2])/true_Λc¹[2]
    Λc²_error_30 = abs(Λc²_30 - true_Λc²[1])/true_Λc²[1]; Λc²_error_60 = abs(Λc²_60 - true_Λc²[2])/true_Λc²[2]
    push!(relative_errors_30, [Λc_error_30, Λc¹_error_30, Λc²_error_30])
    push!(relative_errors_60, [Λc_error_60, Λc¹_error_60, Λc²_error_60])
end

# Plot errors
plt1 = plot(resolutions, [error[1] for error in relative_errors_30], label = "Λc", title = "c = 30")
plot!(resolutions, [error[2] for error in relative_errors_30], label = "Λc¹")
plot!(resolutions, [error[3] for error in relative_errors_30], label = "Λc²")

plt2 = plot(resolutions, [error[1] for error in relative_errors_60], label = "Λc", title = "c = 60")
plot!(resolutions, [error[2] for error in relative_errors_60], label = "Λc¹")
plot!(resolutions, [error[3] for error in relative_errors_60], label = "Λc²")

plt = plot(plt1, plt2, layout = (2, 1))
display(plt)