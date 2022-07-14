push!(LOAD_PATH,"../src/")
using PCEUncertainty
using Pkg; Pkg.add("Documenter")
using Documenter
makedocs(
         sitename = "PCEUncertainty.jl",
         modules  = [PCEUncertainty],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/liamfdoherty/PCEUncertainty.jl",
)