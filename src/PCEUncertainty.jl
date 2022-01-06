module PCEUncertainty

using Random, Distributions
using PolyChaos
using DifferentialEquations

include("generate_collocation.jl")
export generate_collocation_points

end #module