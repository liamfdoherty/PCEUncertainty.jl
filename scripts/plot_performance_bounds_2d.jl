using Plots
using PCEUncertainty

# Acceptable tolerance for relative entropy
B = 0.0484

# Compute the risk-sensitive integrals numerically
c_vals = 1.01:1.:20.01
Λc_vals = []; Λc¹_vals = []; Λc²_vals = []
for c in c_vals
    Λc, Λc¹, Λc² = compute_integrals(F̃, z_weights, c)
    push!(Λc_vals, Λc)
    push!(Λc¹_vals, Λc¹)
    push!(Λc²_vals, Λc²)
end

# Compute performance bounds
performance_bounds_Λc = [B/c + Λc_vals[i] for (i, c) in enumerate(c_vals)]
performance_bounds_Λc¹ = [B/c + Λc¹_vals[i] for (i, c) in enumerate(c_vals)]
performance_bounds_Λc² = [B/c + Λc²_vals[i] for (i, c) in enumerate(c_vals)]

# Plot performance bounds
plt = plot(c_vals, performance_bounds_Λc, label = "B/c + Λc", title = "Performance Bounds", legend = :outertopright)
plot!(c_vals, performance_bounds_Λc¹, label = "B/c + Λc¹")
plot!(c_vals, performance_bounds_Λc², label = "B/c + Λc²")
xlabel!("c")
display(plt)