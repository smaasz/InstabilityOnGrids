using Revise
using TensorKit
using DataFrames
using MultiplesOfPi
import Dates
import JSON3

includet("/Users/stmaas001/Projects/OceanFlows/InstabilityOnGrids/src/grids.jl")
using .FourierSymbols: rot, e₁

nk = 100 # number of wavenumbers

f0 = -1.0e-4;
N  = 1.0e-3;
Ri = 0.5;
H  = 4.0e3;

M2 = √(N^2 * f0^2 / Ri);

u0 = -M2 / f0 * H/2;
θU = Pi // 6;
experiment = Dict(
    :a => 3.125,
    :N  => N,
    :f0 => f0,
    :Nz => 64,
    :type => "u0=$u0,direct-advective",
    :grid => "tri-B",
    :Ri => Ri,
    :theta => (Ri > 1.0 ? 0.0 : π/2) + θU,
    :thetaU => θU,
    :u0 => u0,
    :date => Dates.format(Dates.now(), "yy-mm-dd-HH:MM"),
    :ks => range(start=1e-9, stop=π/3.125e3, length=nk),
)
statespacedim = experiment[:Nz] * 5 + 1;

ms = zeros(nk);
vecs = zeros(ComplexF64, nk, statespacedim);

(; f0, N, Ri, a, Nz, theta, thetaU) = NamedTuple(experiment)
a *= 1e3 # conversion to meters
for (i,k) in enumerate(experiment[:ks])
    A = FourierSymbols.build_system_free_surface(e₁*k; a, Ri, N, f₀=f0, Nz, g=1e8, θ=theta, θU=thetaU,u0=u0)
    A_ = convert(Array, convert(TensorMap, A))
    F = eigen(A_)
    ms[i] = real.(F.values[end])
    vecs[i,:] .= F.vectors[:,end]
end
experiment[:vs] = ms
#experiment[:vecs] = vecs
open("/Users/stmaas001/Projects/OceanFlows/InstabilityOnGrids/data/data.jsonl", "a") do f
    jso = JSON3.write(experiment)
    write(f, @sprintf("%s\n", jso))
end

#----------------------------------------should be deprecated----------------------------------------

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
    #"vecs" => fill(zeros(ComplexF64,100,Nz*5+1), 8),
    "type" => ["true", "adv. form", "adv. form", "adv. form", "true", "adv. form", "adv. form", "adv. form"],
    "grid" => fill("tri-B", 8),
    "N"    => fill(N, 8),
    "Ri"   => [100., 100., 100., 100., 0.5, 0.5, 0.5, 0.5],
    "f0"   => fill(f₀, 8),
    "date" => fill(Dates.format(Dates.now(), "yy-mm-dd-HH:MM"), 8),
    "theta" => fill(θ, 8),
    "thetaU" => fill(θU, 8),
    "Nz"   => fill(Nz, 8),
);


for row in eachrow(df)
    a = row[:a] * 1e3 # to m
    Ri = row[:Ri]
    for (i,k) in enumerate(row[:ks])
        A = FourierSymbols.build_system_free_surface(Tensor(Ri < 1 ? [0; k] : [k; 0.], ℂ^2); a, Ri, N, f₀, Nz, g=1e8, θ=θ, θU=θU)
        A_ = convert(Array, convert(TensorMap, A))
        F = eigen(A_)
        ms[i] = real.(F.values[end])
        vecs[i,:] .= F.vectors[:,end]
    end
    row[:vs] = ms
    #row[:vecs] = vecs
    open("/Users/stmaas001/Projects/OceanFlows/InstabilityOnGrids/data/data.jsonl", "a") do f
        jso = JSON3.write(convert(NamedTuple, row))
        write(f, @sprintf("%s\n", jso))
    end
end
