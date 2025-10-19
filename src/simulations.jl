using Revise
import GridOperatorAnalysis: eady_jacobian, bb, compute_phases
import GridOperatorAnalysis
import LinearAlgebra: eigen, norm
using DataFrames
import DrWatson: produce_or_load, @unpack, @dict, @strdict
import Symbolics: substitute, build_function, unwrap, @variables, simplify
import Symbolics
import SymbolicUtils: Postwalk, PassThrough
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# flags
const tofloat = true

# Constants
const f₀ = -1e-4
const g  = 1e9
const N² = 1e-6
const Ri = only(@variables(Ri))
const U  = 1.0
const H  = 4000
const k  = only(@variables(k))
const l  = only(@variables(l))
const le = only(@variables(le))
const 𝕂ᵘ = only(@variables(𝕂ᵘ))
const 𝕂ᵇ = only(@variables(𝕂ᵇ))
const θU = only(@variables(θU))
const β  = only(@variables(β))

const ϕ = compute_phases(k, l, le)

const hmt_schemes = Dict(
    :TriA => [
	:standard => "Standard",
    ],
    :TriB => [
	:asc => "advective form, streamline derivative on cells",
	:avi => "advective form, vector-invariant", 
	:fdv => "flux form, divergence on vertices", 
	:fdcre => "flux form, diverence on cells with reconstruction on edges"
    ],
    :TriC => [
	:ICON => "ICON"
    ],
    :HexC => [
	:MPAS => "MPAS"
    ]
)

const hst_schemes = Dict(
    :TriA => [
        :low => "low",
        :high => "high",
    ],
    :TriB => [
        :low => "low",
        :high => "high",
    ],
    :TriC => [
        :low => "low",
    ],
    :HexC => [
        :low => "low",
        :high => "high",
    ],
)

function analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, hst_scheme, dissip_scheme, le, Nz)
    θ = (Ri > 1 ? 0 : π/2) + θU
    Kmax = min(1e-2, 2/√3*π/le)

    Ks  = range(1e-10, Kmax*1.1, 400)
    iωs = Complex{Float64}[]
    for K in Ks
        k = K * cos(θ)
	l = K * sin(θ)

	jac = -ComplexF64[unwrap(eady_jac[i,j](k, l, Ri, le, f₀, g, N², 𝕂ᵘ, 𝕂ᵇ, θU, β)) for i=1:size(eady_jac,1), j=1:size(eady_jac,2)]

        vals, vecs = eigen(jac)
        push!(iωs, vals[end])
    end

    Dict(
        :Ri            => Ri,
        :θU            => θU,
        :β             => β,
        :𝕂ᵘ            => 𝕂ᵘ,
        :𝕂ᵇ            => 𝕂ᵇ,
        :grid_t        => grid_t,
        :hmt_scheme    => hmt_scheme,
        :hst_scheme    => hst_scheme,
        :dissip_scheme => dissip_scheme,
        :le            => le,
        :Nz            => Nz,
        :Ks            => Ks,
        :iωs           => iωs
    )
end

@inline function f(x::Int)
    nothing
end
@inline function f(x::Real)
    convert(Float64, x)
end
@inline function f(x)
    nothing
end

const rewriter = Postwalk(PassThrough(f))

function runeady(config)
    @unpack grid_t, hmt_scheme, hst_scheme, dissip_scheme, Nz, H = config
    # overwrite constant values of f₀, N², g
    @variables f₀ N² g
    eady_jac = eady_jacobian(Val(grid_t), k, l, le; ϕ, g, f₀, N², Ri, 𝕂ᵘ, 𝕂ᵇ, H, Nz, U, θU, β, hmt_scheme, hst_scheme, dissip_scheme)
    if tofloat
        eady_jac = simplify.(eady_jac; rewriter)
    end
    eady_jac_ex = build_function.(eady_jac, k, l, Ri, le, f₀, g, N², 𝕂ᵘ, 𝕂ᵇ, θU, β; expression=Val{true})
    @strdict(eady_jac_ex)
end

function initialdf()
    df = DataFrame(
        Ri = Float64[],
        θU = Float64[],
        β  = Float64[],
        𝕂ᵘ = Float64[],
        𝕂ᵇ = Float64[],
        grid_t = Symbol[],
        hmt_scheme = Symbol[],
        hst_scheme = Symbol[],
        dissip_scheme = Symbol[],
        le         = Float64[],
        Nz         = Int[],
        Ks         = Vector{Float64}[],
        iωs        = Vector{Complex{Float64}}[],
    )
end

# Run Experiments
function test_galilean_invariance(; Ri, θU=π/6, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz=8)
    @assert grid_t in keys(hmt_schemes)  && hmt_scheme in first.(hmt_schemes[grid_t])

    config = @dict(grid_t, hmt_scheme, Nz, H)
    path   = joinpath(@__DIR__, "..", "data")
    data, file = produce_or_load(runeady, config, path)
    @unpack eady_jac_ex, eady_sys = data
    eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]

    df = initialdf()
    for β in -0.5:0.1:0.5
        eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
        push!(df, eadyinstance)
    end
    df
end

function test_flow_angle(; Ri, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz=8)
    @assert grid_t in keys(hmt_schemes)  && hmt_scheme in first.(hmt_schemes[grid_t])

    config = @dict(grid_t, hmt_scheme, Nz, H)
    path   = joinpath(@__DIR__, "..", "data")
    data, file = produce_or_load(runeady, config, path)
    @unpack eady_jac_ex, eady_sys = data
    eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
    
    df = initialdf()
    for θU in 0:π/24:π/6
        eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
        push!(df, eadyinstance)
    end
    df
end

# grid_t and (or?) hmt_scheme are variable
function test_baroclinic(; θU=π/6, β=0, 𝕂ᵘ, 𝕂ᵇ, le, Nz=8)
    df = initialdf()
    Ri = 100

    for grid_t in [:TriA, :TriB, :TriC, :HexC]
        for hmt_scheme in first.(hmt_schemes[grid_t])

            config = @dict(grid_t, hmt_scheme, Nz, H)
            path   = joinpath(@__DIR__, "..", "data")
            data, file = produce_or_load(runeady, config, path)
            @unpack eady_jac_ex, eady_sys = data
            eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]

            eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
            push!(df, eadyinstance)
        end
    end
    df
end

function test_symmetric(; θU=π/6, β=0, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz=8)
    df = initialdf()
    Ri = 1/2

    for grid_t in [:TriA, :TriB, :TriC, :HexC]
        for hmt_scheme in hmt_schemes[grid_t]

            config = @dict(grid_t, hmt_scheme, Nz, H)
            path   = joinpath(@__DIR__, "..", "data")
            data, file = produce_or_load(runeady, config, path)
            @unpack eady_jac_ex, eady_sys = data
            eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
            
            eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
            push!(df, eadyinstance)
        end
    end
    df
end

# function testall(; θUs, βs, 𝕂ᵘs, 𝕂ᵇs, les, Nz=8)
#     df = initialdf()
#     for grid_t in [:TriA, :TriB, :TriC, :HexC]
#         for hmt_scheme in first.(hmt_schemes[grid_t])

#             config = @dict(grid_t, hmt_scheme, Nz, H)
#             path   = joinpath(@__DIR__, "..", "data")
#             data, file = produce_or_load(runeady, config, path)
#             @unpack eady_jac_ex, eady_sys = data
#             eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]

#             for Ri in [1/2, 100]
#                 for θU in θUs
#                     for β in βs
#                         for 𝕂ᵘ in 𝕂ᵘs
#                             for 𝕂ᵇ in 𝕂ᵇs
#                                 for le in les
#                                     eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
#                                     push!(df, eadyinstance)
#                                 end
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     df
# end

function testall(; nθUs, nβs, 𝕍ᵘs, 𝕍ᵇs, les, Nz=16)
    df = initialdf()
    for grid_t in [:TriA, :TriB, :TriC, :HexC]
        for hmt_scheme in first.(hmt_schemes[grid_t])
            for hst_scheme in first.(hst_schemes[grid_t])
                for dissip_scheme in [:biharmonic]
                    config = @dict(grid_t, hmt_scheme, hst_scheme, dissip_scheme, Nz, H)
                    path   = joinpath(@__DIR__, "..", "data")
                    data, file = produce_or_load(runeady, config, path)
                    @unpack eady_jac_ex, = data
                    eady_jac = [@RuntimeGeneratedFunction(eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
                    for nθU in nθUs
                        θU = nθU * π/12
                        for nβ in nβs
                            β  = nβ  * 0.1
                            for Ri in [1/2, 100]
                                for 𝕍ᵘ in 𝕍ᵘs
                                    for 𝕍ᵇ in 𝕍ᵇs
                                        for le in les
                                            𝕂ᵘ = dissip_scheme == :harmonic ? 𝕍ᵘ * le : 𝕍ᵘ * le^3
                                            𝕂ᵇ = dissip_scheme == :harmonic ? 𝕍ᵇ * le : 𝕍ᵇ * le^3
                                            eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, hst_scheme, dissip_scheme, le, Nz)
                                            push!(df, eadyinstance)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    df
end

