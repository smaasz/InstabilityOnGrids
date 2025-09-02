using Revise
import GridOperatorAnalysis: eady_jacobian, bb
import GridOperatorAnalysis
import LinearAlgebra: eigen, norm
using DataFrames
import DrWatson: produce_or_load, @unpack, @dict, @strdict
import Symbolics: substitute, build_function, unwrap
import Symbolics
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Constants
const f₀ = -1e-4
const g  = 1e9
const N² = 1e-6
const U  = 1.0
const H  = 4000

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


# function analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
#     θ = (Ri > 1 ? 0 : π/2) + θU
#     Kmax = let  # check if correct
# 	_bb = Symbolics.unwrap.(substitute.(bb, Ref(Dict(GridOperatorAnalysis.sqrt3 => √3, GridOperatorAnalysis.le => le))))
# 	π / norm(_bb * [cos(θ); sin(θ)], Inf)
#     end

#     Ks  = range(1e-10, Kmax, 400)
#     iωs = Complex{Float64}[]
#     for K in Ks
#         k = K * cos(θ)
# 	l = K * sin(θ)

# 	jac = [Complex{Float64}(eady_jac[i,j](k, l, le, Ri, N², g, f₀, 𝕂ᵘ, 𝕂ᵇ, θU, β)) for i=1:size(eady_jac,1), j=1:size(eady_jac,2)]

#         vals, vecs = eigen(jac)
#         push!(iωs, vals[end])
#     end

#     Dict(
#         :Ri         => Ri,
#         :θU         => θU,
#         :β          => β,
#         :𝕂ᵘ         => 𝕂ᵘ,
#         :𝕂ᵇ         => 𝕂ᵇ,
#         :grid_t     => grid_t,
#         :hmt_scheme => hmt_scheme,
#         :le         => le,
#         :Nz         => Nz,
#         :Ks         => Ks,
#         :iωs        => iωs
#     )
# end

function analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, dissip_scheme, le, Nz)
    θ = (Ri > 1 ? 0 : π/2) + θU
    Kmax = let  # check if correct
	_bb = unwrap.(substitute.(bb, Ref(Dict(GridOperatorAnalysis.sqrt3 => √3, GridOperatorAnalysis.le => le))))
	π / norm(_bb * [cos(θ); sin(θ)], Inf)
    end

    Ks  = range(1e-10, Kmax, 400)
    iωs = Complex{Float64}[]
    for K in Ks
        k = K * cos(θ)
	l = K * sin(θ)

	jac = ComplexF64[unwrap(eady_jac[i,j](k, l, Ri, le, 𝕂ᵘ, 𝕂ᵇ)) for i=1:size(eady_jac,1), j=1:size(eady_jac,2)]

        vals, vecs = eigen(jac)
        push!(iωs, vals[end])
    end

    Dict(
        :Ri         => Ri,
        :θU         => θU,
        :β          => β,
        :𝕂ᵘ         => 𝕂ᵘ,
        :𝕂ᵇ         => 𝕂ᵇ,
        :grid_t     => grid_t,
        :hmt_scheme => hmt_scheme,
        :le         => le,
        :Nz         => Nz,
        :Ks         => Ks,
        :iωs        => iωs
    )
end

# function runeady(config)
#     @unpack grid_t, hmt_scheme, Nz, H = config
#     eady_jac_ex, eady_sys = eady_jacobian(Val(grid_t); Nz, H, U, hmt_scheme)
#     @strdict(eady_jac_ex, eady_sys)    
# end

function runeady(config)
    @unpack grid_t, hmt_scheme, dissip_scheme, nθU, nβ, Nz, H
    θU = nθU * π/12
    nβ = nβ  * 0.1
    @parameters Ri le 𝕂ᵘ 𝕂ᵇ
    eady_jac = eady_jacobian_wo_mtk(Val(grid_t), le; g, f₀, Ri, N², 𝕂ᵘ, 𝕂ᵇ, θU, β, Nz, H, U, hmt_scheme, dissip_scheme)
    eady_jac_ex = build_function.(eady_jac, GridOperatorAnalysis.k, GridOperatorAnalysis.l, Ri, le, 𝕂ᵘ, 𝕂ᵇ; expression=Val{false})
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

function testall(; nθUs, nβs, 𝕂ᵘs, 𝕂ᵇs, les, Nz=16)
    df = initialdf()
    for grid_t in [:TriA, :TriB, :TriC, :HexC]
        for hmt_scheme in first.(hmt_schemes[grid_t])
            for nθU in nθUs
                θU = nθU * π/12
                for nβ in nβs
                    β  = nβ  * 0.1
                    config = @dict(grid_t, hmt_scheme, Nz, H, nθU, nβ)
                    path   = joinpath(@__DIR__, "..", "data")
                    data, file = produce_or_load(runeady, config, path)
                    @unpack eady_jac_ex = data
                    eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
                    for Ri in [1/2, 100]
                        for 𝕂ᵘ in 𝕂ᵘs
                            for 𝕂ᵇ in 𝕂ᵇs
                                for le in les
                                    eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, dissip_scheme, le, Nz)
                                    push!(df, eadyinstance)
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

