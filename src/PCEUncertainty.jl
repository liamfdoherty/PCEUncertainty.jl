module PCEUncertainty

using PolyChaos
using DifferentialEquations
using LinearAlgebra

include("structures.jl")
export StochasticODEProblem

include("generate_collocation.jl")
export generate_collocation, generate_collocation_2d

include("evaluate_basis.jl")
export evaluate_basis

include("generate_pce_coefficients.jl")
export generate_pce_coefficients

include("compute_pce_approximation.jl")
export compute_pce_approximation

include("compute_integrals.jl")
export compute_integrals

end #module