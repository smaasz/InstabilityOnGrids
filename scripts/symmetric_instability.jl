using DataFrames
using Revise
using TensorKit

includet("/Users/stmaas001/Projects/OceanFlows/SemiAnalyticInstabilityAnalysis/src/grids.jl")

N  = 1.0e-3 #"s^-1"
Ri = 1/2
f₀ = -1.0e-4 #"s^-1"
M2 = √(N^2 * f₀^2 / Ri)

ms = zeros(100,5);

df = DataFrame(
    "ks"   => [
        range(start=1e-9, stop=π/3.125e3, length=100),
        range(start=1e-9, stop=π/3.125e3, length=100),
        range(start=1e-9, stop=π/6.25e3, length=100),
        range(start=1e-9, stop=π/12.5e3, length=100),
    ],
    "a"    => [floatmin(Float64), 3.125, 6.25, 12.5],
    "vs"   => [zeros(100), zeros(100), zeros(100), zeros(100)],
    "type" => ["true", "adv. form", "adv. form", "adv. form"],
    "grid" => ["B", "B", "B", "B"],
    "N"    => [N, N, N, N],
    "Ri"   => [Ri, Ri, Ri, Ri],
    "f₀"   => [f₀, f₀, f₀, f₀],
);

for row in eachrow(df)
    a = row[:a] * 1e3 # to m
    for (i,k) in enumerate(row[:ks])
        A = FourierSymbols.build_system_free_surface(Tensor([0.; k], ℂ^2); a, Ri, N=N, f₀=f₀, Nz=64, g=1e8)
        A_ = convert(Array, convert(TensorMap, A))
        F = eigen(A_)
        ms[i,:] .= real.(F.values[end-4:end])
    end
    row[:vs] = ms[:,end]
end

open("/Users/stmaas001/Projects/OceanFlows/SemiAnalyticInstabilityAnalysis/data/symmetric_instability.json", "w") do f
    write(f, objecttable(df))
end
