using Revise
using TensorKit

includet("/Users/stmaas001/Projects/OceanFlows/SemiAnalyticInstabilityAnalysis/src/grids.jl")

N  = 1.0e-3 #"s^-1"
f₀ = -1.0e-4 #"s^-1"
Nz = 32

df = DataFrame(
    "ks"   => [
        range(start=1e-9, stop=π/6.25e3, length=100),
        range(start=1e-9, stop=π/6.25e3, length=100),
        range(start=1e-9, stop=π/12.5e3, length=100),
        range(start=1e-9, stop=π/25.e3, length=100),
        range(start=1e-9, stop=π/3.125e3, length=100),
        range(start=1e-9, stop=π/3.125e3, length=100),
        range(start=1e-9, stop=π/6.25e3, length=100),
        range(start=1e-9, stop=π/12.5e3, length=100),
    ],
    "a"    => [floatmin(Float64), 6.25, 12.5, 25.0, floatmin(Float64), 3.125, 6.25, 12.5],
    "vs"   => fill(zeros(100), 8),
    "vecs" => fill(zeros(ComplexF64,100,Nz*5+1), 8),
    "type" => ["true", "adv. form", "adv. form", "adv. form", "true", "adv. form", "adv. form", "adv. form"],
    "grid" => fill("B", 8),
    "N"    => fill(N, 8),
    "Ri"   => [100., 100., 100., 100., 0.5, 0.5, 0.5, 0.5],
    "f₀"   => fill(f₀, 8),
);

ms = zeros(100,5);
vecs = zeros(ComplexF64, 100, Nz * 5 + 1);
for row in eachrow(df)
    a = row[:a] * 1e3 # to m
    Ri = row[:Ri]
    for (i,k) in enumerate(row[:ks])
        A = FourierSymbols.build_system_free_surface(Tensor(Ri < 1 ? [0; k] : [k; 0.], ℂ^2); a, Ri, N, f₀, Nz, g=1e8)
        A_ = convert(Array, convert(TensorMap, A))
        F = eigen(A_)
        ms[i,:] .= real.(F.values[end-4:end])
        vecs[i,:] .= F.vectors[:,end]
    end
    row[:vs] = ms[:,end]
    row[:vecs] = vecs
end

open("/Users/stmaas001/Projects/OceanFlows/SemiAnalyticInstabilityAnalysis/data/baroclinic_instability.json", "w") do f
    write(f, objecttable(df))
end
