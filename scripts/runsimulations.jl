using Revise
using CairoMakie
includet(joinpath("..", "src", "simulations.jl"))

#df = test_galilean_invariance(; Ri=100, 𝕂ᵘ=100.0, 𝕂ᵇ=100.0, grid_t=:TriA, hmt_scheme=:standard, le=6.25e3)

df = let 
    θU = π/6
    β  = 0.0
    𝕂ᵘ = 100.0
    𝕂ᵇ = 100.0
    le = 6.25e3 
    Nz = 8
    test_baroclinic(; θU, β, 𝕂ᵘ, 𝕂ᵇ, le, Nz)
end

