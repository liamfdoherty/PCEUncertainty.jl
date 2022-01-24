using Plots

z₁_vals = LinRange(0, 1, 100); z₂_vals = LinRange(0, 1, 100)

truth(z₁, z₂) = z₂*exp(-z₁) # True value of solution of ODE with 2 uncertainties (no application of F (i.e., no squaring or indicator function applied))
surface(z₁_vals, z₂_vals, truth) # Plot the true solution over all values of the uncertain parameters
scatter!((1/2, 1/2, 3)) # Testing the ability to scatter points of the form (parameter 1, parameter 2, final solution to ODE)