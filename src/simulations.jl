using Base.Threads
using DataFrames
import DrWatson: produce_or_load, @unpack, @dict, @strdict
using Revise
import GridOperatorAnalysis: eady_background_flow, bb, e, c, v, nS, dims, t, z, sqrt3
import GridOperatorAnalysis: TriAFlow, TriAFlowFT, TriBFlow, TriBFlowFT, TriCFlow, TriCFlowFT, HexCFlow, HexCFlowFT
import GridOperatorAnalysis: colpt_type, colptidx, compute_phases, sqrt3_subs
import GridOperatorAnalysis: construct_lcc_sys, fourier_transform_expression, fourier_transform_sys, fouriersymbols
import GridOperatorAnalysis.TriA
import GridOperatorAnalysis.TriB
import GridOperatorAnalysis.TriC
import GridOperatorAnalysis.HexC
import GridOperatorAnalysis
import LinearAlgebra: eigen, I, Diagonal
import ProgressLogging: @progress
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
import Symbolics: substitute, @variables, taylor_coeff, simplify, coeff, taylor, expand, det, Num, expand_derivatives, Differential, build_function
import Symbolics
import SymbolicUtils: Postwalk, PassThrough


# Parameters of the experiment
# @variables f₀ g N² Ri M² β θU 𝕂ᵘ 𝕂ᵇ h a H Nz

# const f₀ = only(@variables(f₀))
# const g  = only(@variables(g))
# const N² = only(@variables(N²))
# const Ri = only(@variables(Ri))
# const M² = only(@variables(M²))
# const β  = only(@variables(β))
# const θU = only(@variables(θU))
# const 𝕂ᵘ = only(@variables(𝕂ᵘ))
# const 𝕂ᵇ = only(@variables(𝕂ᵇ))

const a  = only(@variables(a))
const h  = only(@variables(h))

# Wavenumbers
#@variables k l
const k = only(@variables(k))
const l = only(@variables(l))

const ϕ = compute_phases(k, l, a)

# vertical variable
# const z = only(@variables(z))

# Pertubation parameter
#const ϵ = only(@variables(ϵ))


function computesymbols(config)
    @unpack grid_t, hmt_scheme, hst_scheme, dissip_scheme = config
    @variables f₀ g N² Ri M² β θU 𝕂ᵘ 𝕂ᵇ H Nz

    flow_t   = getflow(grid_t)
    flowft_t = getflowft(grid_t)
    
    bflow = eady_background_flow(Val(grid_t), a; f₀, N², Ri, θU, β);
    dflow = let
        if colpt_type(flow_t, :u⃗) == :edge
	    @variables (du⃗(t,z))[-nS:nS,-nS:nS,1:dims(colpt_type(flow_t, :u⃗))]
	else
	    @variables (du⃗(t,z))[1:2,-nS:nS,-nS:nS,1:dims(colpt_type(flow_t, :u⃗))]
	end
	@variables begin
            (db(t,z))[-nS:nS,-nS:nS, 1:dims(colpt_type(flow_t, :b))]
            (dη(t))[-nS:nS, -nS:nS, 1:dims(colpt_type(flow_t, :η))]
            (∫∇ᵀdu⃗dz(t,z))[-nS:nS,-nS:nS, 1:dims(colpt_type(flow_t, :∫∇ᵀu⃗dz))]
            (dw(t,z))[-nS:nS, -nS:nS, 1:dims(colpt_type(flow_t, :w))]
            (dp(t,z))[-nS:nS, -nS:nS, 1:dims(colpt_type(flow_t, :p))]
        end
        flow_t{nS}(; u⃗ = du⃗, w = dw, ∫∇ᵀu⃗dz = ∫∇ᵀdu⃗dz, b = db, p = dp, η = dη)
    end

    # construct the linearized constant coefficient system
    sys = construct_lcc_sys(bflow, dflow, a; f₀, g, 𝕂ᵘ, 𝕂ᵇ, hmt_scheme, hst_scheme, dissip_scheme)

    fflow = let
        if colpt_type(flow_t, :u⃗) == :edge
	    @variables (u⃗̂(t,z))[1:dims(colpt_type(flowft_t, :u⃗))]
	else
	    @variables (u⃗̂(t,z))[1:2, 1:dims(colpt_type(flowft_t, :u⃗))]
	end
	@variables begin
            (b̂(t,z))[1:dims(colpt_type(flowft_t, :b))]
            (η̂(t))[1:dims(colpt_type(flowft_t, :η))]
            (∫∇ᵀu⃗̂dz(t))[1:dims(colpt_type(flowft_t, :∫∇ᵀu⃗dz))]
            (ŵ(t,z))[1:dims(colpt_type(flowft_t, :w))]
            (p̂(t,z))[1:dims(colpt_type(flowft_t, :p))]
        end
        fflow = flowft_t{Num}(; u⃗=u⃗̂, w=ŵ, ∫∇ᵀu⃗dz=∫∇ᵀu⃗̂dz, b=b̂, p=p̂, η=η̂)
    end

    ftsys = fourier_transform_sys(Val(grid_t), sys; dflow, fflow, ϕ)

    # fourier symbols
    fsymbols = fouriersymbols(Val(grid_t), ftsys, fflow; k, l, f₀, g, N², M²=√(f₀^2*N²/Ri), θU, β, 𝕂ᵘ, 𝕂ᵇ, a, h, doapprox=false, dophasesubs=false)

    # functions to build fourier symbols
    fsyms = Dict()
    for (name, fsymbol) in pairs(fsymbols)
	fsyms[name] = build_function(fsymbol, z, f₀, N², Ri, θU, β, k, l, a, h; expression=Val{true})
    end 
    @strdict(fsyms)
end


# compute vertical multiplication operators
const vfops = let
    @variables f₀ N² M² Ri θU β H
    Ū  = (z + H * (1//2 + β)) * -M²/f₀
    Ūz = expand_derivatives(Differential(z)(Ū))

    vops = Dict(
        :U  => Ū * cos(θU),
        :V  => Ū * sin(θU),
        :Uz => Ūz * cos(θU),
        :Vz => Ūz * sin(θU),
        :Bx => M² * -sin(θU),
        :By => M² * cos(θU),
    )

    @assert isequal(vops[:Bx], f₀ * vops[:Vz])
    @assert isequal(vops[:By], -f₀ * vops[:Uz])

    (;
        [Symbol("$(name)f") => build_function(substitute(vop, Dict(M² => √(f₀^2 * N²/ Ri))), z, f₀, N², Ri, θU, β, H; expression=Val{false}) for (name, vop) in pairs(vops)]...
    )    
end


# assembling of system matrix

## initialize buffer for the symbols
# const btriA = let
#     du = dims(:vertex)
#     ds = dims(:vertex)
#     Dict(
#         :Gx    => zeros(ComplexF64, 2, du, 2, du),
#         :Gy    => zeros(ComplexF64, 2, du, 2, du),
#         :M     => zeros(ComplexF64, 2, du, 2, du),
#         :A⁽ˣ⁾  => zeros(ComplexF64, 2, du, ds),
#         :A⁽ʸ⁾  => zeros(ComplexF64, 2, du, ds),
#         :G     => zeros(ComplexF64, 2, du, ds),
#         :Dᵘ    => zeros(ComplexF64, 2, du, 2, du),
#         :I     => zeros(ComplexF64, ds, ds),
#         :Av⁽ˣ⁾ => zeros(ComplexF64, ds, 2, du),
#         :Av⁽ʸ⁾ => zeros(ComplexF64, ds, 2, du),
#         :Γx    => zeros(ComplexF64, ds, ds),
#         :Γy    => zeros(ComplexF64, ds, ds),
#         :Dᵇ    => zeros(ComplexF64, ds, ds),
#         :D     => zeros(ComplexF64, ds, 2, du),
#     )
# end
# const btriB = let
#     du = dims(:cell)
#     ds = dims(:vertex)
#     Dict(
#         :Gx    => zeros(ComplexF64, 2, du, 2, du),
#         :Gy    => zeros(ComplexF64, 2, du, 2, du),
#         :M     => zeros(ComplexF64, 2, du, 2, du),
#         :A⁽ˣ⁾  => zeros(ComplexF64, 2, du, ds),
#         :A⁽ʸ⁾  => zeros(ComplexF64, 2, du, ds),
#         :G     => zeros(ComplexF64, 2, du, ds),
#         :Dᵘ    => zeros(ComplexF64, 2, du, 2, du),
#         :I     => zeros(ComplexF64, ds, ds),
#         :Av⁽ˣ⁾ => zeros(ComplexF64, ds, 2, du),
#         :Av⁽ʸ⁾ => zeros(ComplexF64, ds, 2, du),
#         :Γx    => zeros(ComplexF64, ds, ds),
#         :Γy    => zeros(ComplexF64, ds, ds),
#         :Dᵇ    => zeros(ComplexF64, ds, ds),
#         :D     => zeros(ComplexF64, ds, 2, du),
#     )
# end
# const btriC = let
#     du = dims(:edge)
#     ds = dims(:cell)
#     Dict(
#         :Gx    => zeros(ComplexF64, du, du),
#         :Gy    => zeros(ComplexF64, du, du),
#         :M     => zeros(ComplexF64, du, du),
#         :A⁽ˣ⁾  => zeros(ComplexF64, du, ds),
#         :A⁽ʸ⁾  => zeros(ComplexF64, du, ds),
#         :G     => zeros(ComplexF64, du, ds),
#         :Dᵘ    => zeros(ComplexF64, du, du),
#         :I     => zeros(ComplexF64, ds, ds),
#         :Av⁽ˣ⁾ => zeros(ComplexF64, ds, du),
#         :Av⁽ʸ⁾ => zeros(ComplexF64, ds, du),
#         :Γx    => zeros(ComplexF64, ds, ds),
#         :Γy    => zeros(ComplexF64, ds, ds),
#         :Dᵇ    => zeros(ComplexF64, ds, ds),
#         :D     => zeros(ComplexF64, ds, du),
#     )
# end
# const bhexC = let
#     du = dims(:edge)
#     ds = dims(:vertex)
#     Dict(
#         :Gx    => zeros(ComplexF64, du, du),
#         :Gy    => zeros(ComplexF64, du, du),
#         :M     => zeros(ComplexF64, du, du),
#         :A⁽ˣ⁾  => zeros(ComplexF64, du, ds),
#         :A⁽ʸ⁾  => zeros(ComplexF64, du, ds),
#         :G     => zeros(ComplexF64, du, ds),
#         :Dᵘ    => zeros(ComplexF64, du, du),
#         :I     => zeros(ComplexF64, ds, ds),
#         :Av⁽ˣ⁾ => zeros(ComplexF64, ds, du),
#         :Av⁽ʸ⁾ => zeros(ComplexF64, ds, du),
#         :Γx    => zeros(ComplexF64, ds, ds),
#         :Γy    => zeros(ComplexF64, ds, ds),
#         :Dᵇ    => zeros(ComplexF64, ds, ds),
#         :D     => zeros(ComplexF64, ds, du),
#     )
# end


function createbuffer(::Val{grid_t}) where {grid_t}
    if grid_t == :TriA
        du = dims(:vertex)
        ds = dims(:vertex)
        Dict(
            :Gx    => zeros(ComplexF64, 2, du, 2, du),
            :Gy    => zeros(ComplexF64, 2, du, 2, du),
            :M     => zeros(ComplexF64, 2, du, 2, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, 2, du, ds),
            :A⁽ʸ⁾  => zeros(ComplexF64, 2, du, ds),
            :G     => zeros(ComplexF64, 2, du, ds),
            :Dᵘ    => zeros(ComplexF64, 2, du, 2, du),
            :I     => zeros(ComplexF64, ds, ds),
            :Av⁽ˣ⁾ => zeros(ComplexF64, ds, 2, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, ds, 2, du),
            :Γx    => zeros(ComplexF64, ds, ds),
            :Γy    => zeros(ComplexF64, ds, ds),
            :Dᵇ    => zeros(ComplexF64, ds, ds),
            :D     => zeros(ComplexF64, ds, 2, du),
        )
    elseif grid_t == :TriB
        du = dims(:cell)
        ds = dims(:vertex)
        Dict(
            :Gx    => zeros(ComplexF64, 2, du, 2, du),
            :Gy    => zeros(ComplexF64, 2, du, 2, du),
            :M     => zeros(ComplexF64, 2, du, 2, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, 2, du, ds),
            :A⁽ʸ⁾  => zeros(ComplexF64, 2, du, ds),
            :G     => zeros(ComplexF64, 2, du, ds),
            :Dᵘ    => zeros(ComplexF64, 2, du, 2, du),
            :I     => zeros(ComplexF64, ds, ds),
            :Av⁽ˣ⁾ => zeros(ComplexF64, ds, 2, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, ds, 2, du),
            :Γx    => zeros(ComplexF64, ds, ds),
            :Γy    => zeros(ComplexF64, ds, ds),
            :Dᵇ    => zeros(ComplexF64, ds, ds),
            :D     => zeros(ComplexF64, ds, 2, du),
        )
    elseif grid_t == :TriC
        du = dims(:edge)
        ds = dims(:cell)
        Dict(
            :Gx    => zeros(ComplexF64, du, du),
            :Gy    => zeros(ComplexF64, du, du),
            :M     => zeros(ComplexF64, du, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, du, ds),
            :A⁽ʸ⁾  => zeros(ComplexF64, du, ds),
            :G     => zeros(ComplexF64, du, ds),
            :Dᵘ    => zeros(ComplexF64, du, du),
            :I     => zeros(ComplexF64, ds, ds),
            :Av⁽ˣ⁾ => zeros(ComplexF64, ds, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, ds, du),
            :Γx    => zeros(ComplexF64, ds, ds),
            :Γy    => zeros(ComplexF64, ds, ds),
            :Dᵇ    => zeros(ComplexF64, ds, ds),
            :D     => zeros(ComplexF64, ds, du),
        )
    else
        du = dims(:edge)
        ds = dims(:vertex)
        Dict(
            :Gx    => zeros(ComplexF64, du, du),
            :Gy    => zeros(ComplexF64, du, du),
            :M     => zeros(ComplexF64, du, du),
            :A⁽ˣ⁾  => zeros(ComplexF64, du, ds),
            :A⁽ʸ⁾  => zeros(ComplexF64, du, ds),
            :G     => zeros(ComplexF64, du, ds),
            :Dᵘ    => zeros(ComplexF64, du, du),
            :I     => zeros(ComplexF64, ds, ds),
            :Av⁽ˣ⁾ => zeros(ComplexF64, ds, du),
            :Av⁽ʸ⁾ => zeros(ComplexF64, ds, du),
            :Γx    => zeros(ComplexF64, ds, ds),
            :Γy    => zeros(ComplexF64, ds, ds),
            :Dᵇ    => zeros(ComplexF64, ds, ds),
            :D     => zeros(ComplexF64, ds, du),
        )
    end
end

## methods for assembling
function systemmat(grid_t::Union{Val{:TriA}, Val{:TriB}}, fsyms, b, k, l; g, f₀, N², H, Nz, Ri, θU, β, a, Vᵘ, Vᵇ, dissip_scheme, useidealized=Dict())
    Δz  = H / Nz
    h   = a * √3/2

    # conversion to dissipation parameters
    if dissip_scheme == :biharmonic
	𝕂ᵘ = Vᵘ * a^3
	𝕂ᵇ = Vᵇ * a^3
    else
	𝕂ᵘ = Vᵘ * a
	𝕂ᵇ = Vᵇ * a
    end

    # compute fourier symbols
    for (name, fsym) in pairs(fsyms)
        useideal = get(useidealized, name, false)
	fsym[2](b[name], z, f₀, N², Ri, θU, β, k, l, useideal ? 1e-20 : a, useideal ? √3/2*1e-20 : h)
    end

    # vertical operators
    @unpack Uf, Vf, Uzf, Vzf, Bxf, Byf = vfops
    U̲ = Diagonal([Uf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲ = Diagonal([Vf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗ = (U̲, V̲)
    U̲z = Diagonal([Uzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲z = Diagonal([Vzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗z = (U̲z, V̲z)
    B̲x = Diagonal([Bxf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲y = Diagonal([Byf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲⃗ = (B̲x, B̲y)
    
    flow_t = getflow(grid_t)
    du = dims(colpt_type(flow_t, :u⃗))
    ds = dims(colpt_type(flow_t, :b))
    W̲ = let
	M = Δz * 1/2 * [iV ≤ iVi ≤ iV+1 for iV=1:Nz, iVi=1:Nz+1] * [iV < iVi ? 1 : 0 for iVi=1:Nz+1, iV=1:Nz]
	[kron(-b[:D][:,jTH,:], M) for jTH=1:2]
    end
    
    P̲ = let
	M = Δz * 1/2 * [iV ≤ iVi ? 1 : 0 for iV=1:Nz, iVi=1:Nz-1] * [iVi ≤ iV ≤ iVi+1 for iVi=1:Nz-1, iV=1:Nz]
	kron(I(ds), -M)
    end

    # assembling	
    S̲ = zeros(ComplexF64, (2*du+ds)*Nz+ds, (2*du+ds)*Nz+ds)

    ru⃗ = [(iTH-1)*du*Nz+1:iTH*du*Nz for iTH=1:2]
    rb = 2*du*Nz+1:(2*du+ds)*Nz
    rη = (2*du+ds)*Nz+1:(2*du+ds)*Nz+ds
    
    # U⃗
    @views for iTH = 1:2
	for jTH = 1:2
	    kron!(S̲[ru⃗[iTH], ru⃗[jTH]], b[:Gx][iTH,:,jTH,:], U̲)
	    S̲[ru⃗[iTH], ru⃗[jTH]] += kron(b[:Gy][iTH,:,jTH,:], V̲)
	    S̲[ru⃗[iTH], ru⃗[jTH]] += kron(b[:A⁽ˣ⁾][iTH,:,:], U̲z) * W̲[jTH]
	    S̲[ru⃗[iTH], ru⃗[jTH]] += kron(b[:A⁽ʸ⁾][iTH,:,:], V̲z) * W̲[jTH]
	    S̲[ru⃗[iTH], ru⃗[jTH]] += f₀ * kron(b[:M][iTH,:,jTH,:], I(Nz))
	    S̲[ru⃗[iTH], ru⃗[jTH]] += 𝕂ᵘ * kron(b[:Dᵘ][iTH,:,jTH,:], I(Nz))
	end

	kron!(S̲[ru⃗[iTH], rb], b[:G][iTH,:,:], I(Nz))
	S̲[ru⃗[iTH], rb] *= P̲

	kron!(S̲[ru⃗[iTH], rη], b[:G][iTH,:,:], g*ones(Nz,1))
    end

    # b
    @views for jTH = 1:2
	S̲[rb, ru⃗[jTH]] += N² * kron(b[:I], I(Nz)) * W̲[jTH]
	S̲[rb, ru⃗[jTH]] += kron(b[:Av⁽ˣ⁾][:,jTH,:], B̲x)
	S̲[rb, ru⃗[jTH]] += kron(b[:Av⁽ʸ⁾][:,jTH,:], B̲y)
    end
    @views kron!(S̲[rb, rb], b[:Γx], U̲)
    @views S̲[rb, rb] += kron(b[:Γy], V̲)
    @views S̲[rb, rb] += 𝕂ᵇ * kron(b[:Dᵇ], I(Nz))

    # η
    @views for jTH = 1:2
	kron!(S̲[rη, ru⃗[jTH]], b[:D][:,jTH,:], Δz*ones(1,Nz))
    end

    S̲[rη, rη] .= U̲[end] * b[:Γx] + V̲[end] * b[:Γy]
    
    S̲
end

function systemmat(grid_t::Union{Val{:TriC}, Val{:HexC}}, fsyms, b, k, l; g, f₀, N², H, Nz, Ri, θU, β, a, Vᵘ, Vᵇ, dissip_scheme, useidealized=Dict())
    Δz  = H / Nz
    h   = a * √3/2

    # conversion to dissipation parameters
    if dissip_scheme == :biharmonic
	𝕂ᵘ = Vᵘ * a^3
	𝕂ᵇ = Vᵇ * a^3
    else
	𝕂ᵘ = Vᵘ * a
	𝕂ᵇ = Vᵇ * a
    end

    # compute fourier symbols
    for (name, fsym) in pairs(fsyms)
        useideal = get(useidealized, name, false)
	fsym[2](b[name], z, f₀, N², Ri, θU, β, k, l, useideal ? 1e-20 : a, useideal ? √3/2*1e-20 : h)
    end

    # vertical operators
    @unpack Uf, Vf, Uzf, Vzf, Bxf, Byf = vfops
    U̲ = Diagonal([Uf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲ = Diagonal([Vf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗ = (U̲, V̲)
    U̲z = Diagonal([Uzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    V̲z = Diagonal([Vzf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    U̲⃗z = (U̲z, V̲z)
    B̲x = Diagonal([Bxf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲y = Diagonal([Byf(((iV-1/2)-Nz) * Δz, f₀, N², Ri, θU, β, H) for iV=1:Nz])
    B̲⃗ = (B̲x, B̲y)

    flow_t = getflow(grid_t)
    du = dims(colpt_type(flow_t, :u⃗))
    ds = dims(colpt_type(flow_t, :b))
    W̲ = let
	M = Δz * 1/2 * [iV ≤ iVi ≤ iV+1 for iV=1:Nz, iVi=1:Nz+1] * [iV < iVi ? 1 : 0 for iVi=1:Nz+1, iV=1:Nz]
	kron(-b[:D], M)
    end
    
    P̲ = let
	M = Δz * 1/2 * [iV ≤ iVi ? 1 : 0 for iV=1:Nz, iVi=1:Nz-1] * [iVi ≤ iV ≤ iVi+1 for iVi=1:Nz-1, iV=1:Nz]
	kron(I(ds), -M)
    end

    # assembling	
    S̲ = zeros(ComplexF64, (du+ds)*Nz+ds, (du+ds)*Nz+ds)

    ru⃗, rb, rη = (1:du*Nz, du*Nz+1:(du+ds)*Nz, (du+ds)*Nz+1:(du+ds)*Nz+ds)
    # U⃗

    @views begin
	kron!(S̲[ru⃗, ru⃗], b[:Gx], U̲)
	S̲[ru⃗, ru⃗] += kron(b[:Gy], V̲)
	S̲[ru⃗, ru⃗] += kron(b[:A⁽ˣ⁾], U̲z) * W̲
	S̲[ru⃗, ru⃗] += kron(b[:A⁽ʸ⁾], V̲z) * W̲
	S̲[ru⃗, ru⃗] += f₀ * kron(b[:M], I(Nz))
	S̲[ru⃗, ru⃗] += 𝕂ᵘ * kron(b[:Dᵘ], I(Nz))
	
	kron!(S̲[ru⃗, rb], b[:G], I(Nz))
	S̲[ru⃗, rb] *= P̲
	
	kron!(S̲[ru⃗, rη], b[:G], g * ones(Nz,1))
    end

    # b
    @views begin
	S̲[rb, ru⃗] += N² * kron(b[:I], I(Nz)) * W̲
	S̲[rb, ru⃗] += kron(b[:Av⁽ˣ⁾], B̲x)
	S̲[rb, ru⃗] += kron(b[:Av⁽ʸ⁾], B̲y)

	kron!(S̲[rb, rb], b[:Γx], U̲)
	S̲[rb, rb] += kron(b[:Γy], V̲)
	S̲[rb, rb] += 𝕂ᵇ * kron(b[:Dᵇ], I(Nz))
    end

    # η
    @views begin
	kron!(S̲[rη, ru⃗], b[:D], Δz*ones(1,Nz))
	S̲[rη, rη] .+= U̲[end] * b[:Γx] + V̲[end] * b[:Γy]
    end
    
    S̲
end


function analyzeinstability(config, fsyms; kwargs...)
    @unpack grid_t, hmt_scheme, hst_scheme, dissip_scheme, g, f₀, N², Ri, θU, β, Vᵘ, Vᵇ, a, H, Nz = config
    θ = (Ri > 1 ? 0 : π/2) + θU
    Kmax = min(2/√3*π/6.25e3, 2/√3*π/a)

    nK = 500
    Ks  = range(1e-10, Kmax*1.1, nK)
    iωs = zeros(ComplexF64, nK)
    bfsyms = createbuffer(Val(grid_t))
    @progress for (iK, K) in enumerate(Ks)
        k = K * cos(θ)
	l = K * sin(θ)

        S̲ = systemmat(Val(grid_t), fsyms, bfsyms, k, l; g, f₀, N², H, Nz, Ri, θU, β, a, Vᵘ, Vᵇ, dissip_scheme, kwargs...)

        F = eigen(-S̲)
        iωs[iK] = F.values[end]
    end

    instance = copy(config)
    
    instance[:Ks]  = Ks
    instance[:iωs] = iωs
    instance
end

function eadyexperiments(; nθUs, nβs, Vᵘs, Vᵇs, as, Nz=16, H=4000.0, f₀=-1e-4, g=1e9, N²=1e-6)
    nt = Threads.nthreads()
    dfs = [initialdf() for i in 1:nt]
    #df = initialdf()
    Threads.@threads for grid_t in [:TriA, :TriB, :TriC, :HexC]
        Threads.@threads for hmt_scheme in first.(hmt_schemes[grid_t])
            Threads.@threads for hst_scheme in first.(hst_schemes[grid_t])
                Threads.@threads for dissip_scheme in [:biharmonic]
                    config = @dict(grid_t, hmt_scheme, hst_scheme, dissip_scheme)
                    path   = joinpath(@__DIR__, "..", "data")
                    data, file = produce_or_load(computesymbols, config, path)
                    @unpack fsyms, = data
                    fsyms_generated = Dict([name => (@RuntimeGeneratedFunction(fsym[1]), @RuntimeGeneratedFunction(fsym[2])) for (name, fsym) in pairs(fsyms)])
                    Threads.@threads for nθU in nθUs
                        θU = nθU * π/12
                        Threads.@threads for nβ in nβs
                            β  = nβ  * 0.1
                            Threads.@threads for Ri in [0.5, 100.0]
                                Threads.@threads for Vᵘ in Vᵘs
                                    Threads.@threads for Vᵇ in Vᵇs
                                        Threads.@threads for a in as
                                            fullconfig = @dict(grid_t, hmt_scheme, hst_scheme, dissip_scheme, g, f₀, N², Ri, θU, β, Vᵘ, Vᵇ, a, H, Nz)
                                             eadyinstance = analyzeinstability(fullconfig, fsyms_generated)
                                            df = dfs[Threads.threadid()]
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
    vcat(dfs...)
end

function symbolicsymbols(grid_t, hmt_scheme, hst_scheme, dissip_scheme, doapprox=false, dophasesubs=false)
    config = @dict(grid_t, hmt_scheme, hst_scheme, dissip_scheme)
    path   = joinpath(@__DIR__, "..", "data")
    data, file = produce_or_load(computesymbols, config, path)
    @unpack fsyms, = data
    fsyms_generated = Dict([name => (@RuntimeGeneratedFunction(fsym[1]), @RuntimeGeneratedFunction(fsym[2])) for (name, fsym) in pairs(fsyms)])

    @variables f₀ g N² Ri M² β θU 𝕂ᵘ 𝕂ᵇ H Nz
    symbolicsyms = Dict([
        name => fsym[1](z, f₀, N², Ri, θU, β, k, l, a, h) for (name, fsym) in pairs(fsyms_generated)# if name ≠ :Γy
    ])
    symbolicsyms
end


# some helper functions (not optimal)
function initialdf()
    df = DataFrame(
        g             = Float64[],
        f₀            = Float64[],
        N²            = Float64[],
        Ri            = Float64[],
        θU            = Float64[],
        β             = Float64[],
        Vᵘ            = Float64[],
        Vᵇ            = Float64[],
        grid_t        = Symbol[],
        hmt_scheme    = Symbol[],
        hst_scheme    = Symbol[],
        dissip_scheme = Symbol[],
        a             = Float64[],
        H             = Float64[],
        Nz            = Int[],
        Ks            = Vector{Float64}[],
        iωs           = Vector{Complex{Float64}}[],
    )
end

function getflow(grid_t::Symbol)
    if grid_t == :TriA
	TriAFlow
    elseif grid_t == :TriB
	TriBFlow
    elseif grid_t == :TriC
	TriCFlow
    else
	HexCFlow
    end
end

function getflow(::Val{grid_t}) where {grid_t}
    if grid_t == :TriA
	TriAFlow
    elseif grid_t == :TriB
	TriBFlow
    elseif grid_t == :TriC
	TriCFlow
    else
	HexCFlow
    end
end

function getflowft(grid_t)
    if grid_t == :TriA
	TriAFlowFT
    elseif grid_t == :TriB
	TriBFlowFT
    elseif grid_t == :TriC
	TriCFlowFT
    else
	HexCFlowFT
    end
end

const hmt_schemes = Dict(
    :TriA => [
	:standard => "Standard",
    ],
    :TriB => [
	:asc => "advective form, streamline derivative on cells",
	:avi => "advective form, vector-invariant", 
	:fdv => "flux form, divergence on vertices", 
	#:fdcre => "flux form, diverence on cells with reconstruction on edges"
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

# function analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, hst_scheme, dissip_scheme, le, Nz)
#     θ = (Ri > 1 ? 0 : π/2) + θU
#     Kmax = min(1e-2, 2/√3*π/le)

#     Ks  = range(1e-10, Kmax*1.1, 400)
#     iωs = Complex{Float64}[]
#     for K in Ks
#         k = K * cos(θ)
# 	l = K * sin(θ)

# 	jac = -ComplexF64[unwrap(eady_jac[i,j](k, l, Ri, le, f₀, g, N², 𝕂ᵘ, 𝕂ᵇ, θU, β)) for i=1:size(eady_jac,1), j=1:size(eady_jac,2)]

#         vals, vecs = eigen(jac)
#         push!(iωs, vals[end])
#     end

#     Dict(
#         :Ri            => Ri,
#         :θU            => θU,
#         :β             => β,
#         :𝕂ᵘ            => 𝕂ᵘ,
#         :𝕂ᵇ            => 𝕂ᵇ,
#         :grid_t        => grid_t,
#         :hmt_scheme    => hmt_scheme,
#         :hst_scheme    => hst_scheme,
#         :dissip_scheme => dissip_scheme,
#         :le            => le,
#         :Nz            => Nz,
#         :Ks            => Ks,
#         :iωs           => iωs
#     )
# end

# @inline function f(x::Int)
#     nothing
# end
# @inline function f(x::Real)
#     convert(Float64, x)
# end
# @inline function f(x)
#     nothing
# end

# const rewriter = Postwalk(PassThrough(f))

# function runeady(config)
#     @unpack grid_t, hmt_scheme, hst_scheme, dissip_scheme, Nz, H = config
#     # overwrite constant values of f₀, N², g
#     @variables f₀ N² g
#     eady_jac = eady_jacobian(Val(grid_t), k, l, le; ϕ, g, f₀, N², Ri, 𝕂ᵘ, 𝕂ᵇ, H, Nz, U, θU, β, hmt_scheme, hst_scheme, dissip_scheme)
#     if tofloat
#         eady_jac = simplify.(eady_jac; rewriter)
#     end
#     eady_jac_ex = build_function.(eady_jac, k, l, Ri, le, f₀, g, N², 𝕂ᵘ, 𝕂ᵇ, θU, β; expression=Val{true})
#     @strdict(eady_jac_ex)
# end

# # Run Experiments
# function test_galilean_invariance(; Ri, θU=π/6, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz=8)
#     @assert grid_t in keys(hmt_schemes)  && hmt_scheme in first.(hmt_schemes[grid_t])

#     config = @dict(grid_t, hmt_scheme, Nz, H)
#     path   = joinpath(@__DIR__, "..", "data")
#     data, file = produce_or_load(runeady, config, path)
#     @unpack eady_jac_ex, eady_sys = data
#     eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]

#     df = initialdf()
#     for β in -0.5:0.1:0.5
#         eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
#         push!(df, eadyinstance)
#     end
#     df
# end

# function test_flow_angle(; Ri, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz=8)
#     @assert grid_t in keys(hmt_schemes)  && hmt_scheme in first.(hmt_schemes[grid_t])

#     config = @dict(grid_t, hmt_scheme, Nz, H)
#     path   = joinpath(@__DIR__, "..", "data")
#     data, file = produce_or_load(runeady, config, path)
#     @unpack eady_jac_ex, eady_sys = data
#     eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
    
#     df = initialdf()
#     for θU in 0:π/24:π/6
#         eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
#         push!(df, eadyinstance)
#     end
#     df
# end

# # grid_t and (or?) hmt_scheme are variable
# function test_baroclinic(; θU=π/6, β=0, 𝕂ᵘ, 𝕂ᵇ, le, Nz=8)
#     df = initialdf()
#     Ri = 100

#     for grid_t in [:TriA, :TriB, :TriC, :HexC]
#         for hmt_scheme in first.(hmt_schemes[grid_t])

#             config = @dict(grid_t, hmt_scheme, Nz, H)
#             path   = joinpath(@__DIR__, "..", "data")
#             data, file = produce_or_load(runeady, config, path)
#             @unpack eady_jac_ex, eady_sys = data
#             eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]

#             eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
#             push!(df, eadyinstance)
#         end
#     end
#     df
# end

# function test_symmetric(; θU=π/6, β=0, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz=8)
#     df = initialdf()
#     Ri = 1/2

#     for grid_t in [:TriA, :TriB, :TriC, :HexC]
#         for hmt_scheme in hmt_schemes[grid_t]

#             config = @dict(grid_t, hmt_scheme, Nz, H)
#             path   = joinpath(@__DIR__, "..", "data")
#             data, file = produce_or_load(runeady, config, path)
#             @unpack eady_jac_ex, eady_sys = data
#             eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
            
#             eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
#             push!(df, eadyinstance)
#         end
#     end
#     df
# end

# # function testall(; θUs, βs, 𝕂ᵘs, 𝕂ᵇs, les, Nz=8)
# #     df = initialdf()
# #     for grid_t in [:TriA, :TriB, :TriC, :HexC]
# #         for hmt_scheme in first.(hmt_schemes[grid_t])

# #             config = @dict(grid_t, hmt_scheme, Nz, H)
# #             path   = joinpath(@__DIR__, "..", "data")
# #             data, file = produce_or_load(runeady, config, path)
# #             @unpack eady_jac_ex, eady_sys = data
# #             eady_jac = [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]

# #             for Ri in [1/2, 100]
# #                 for θU in θUs
# #                     for β in βs
# #                         for 𝕂ᵘ in 𝕂ᵘs
# #                             for 𝕂ᵇ in 𝕂ᵇs
# #                                 for le in les
# #                                     eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, le, Nz)
# #                                     push!(df, eadyinstance)
# #                                 end
# #                             end
# #                         end
# #                     end
# #                 end
# #             end
# #         end
# #     end
# #     df
# # end

# function testall(; nθUs, nβs, 𝕍ᵘs, 𝕍ᵇs, les, Nz=16)
#     df = initialdf()
#     for grid_t in [:TriA]#[:TriA, :TriB, :TriC, :HexC]
#         for hmt_scheme in first.(hmt_schemes[grid_t])
#             for hst_scheme in first.(hst_schemes[grid_t])
#                 for dissip_scheme in [:biharmonic]
#                     config = @dict(grid_t, hmt_scheme, hst_scheme, dissip_scheme, Nz, H)
#                     path   = joinpath(@__DIR__, "..", "data")
#                     data, file = produce_or_load(runeady, config, path)
#                     @unpack eady_jac_ex, = data
#                     eady_jac = [@RuntimeGeneratedFunction(eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
#                     for nθU in nθUs
#                         θU = nθU * π/12
#                         for nβ in nβs
#                             β  = nβ  * 0.1
#                             for Ri in [1/2, 100]
#                                 for 𝕍ᵘ in 𝕍ᵘs
#                                     for 𝕍ᵇ in 𝕍ᵇs
#                                         for le in les
#                                             𝕂ᵘ = dissip_scheme == :harmonic ? 𝕍ᵘ * le : 𝕍ᵘ * le^3
#                                             𝕂ᵇ = dissip_scheme == :harmonic ? 𝕍ᵇ * le : 𝕍ᵇ * le^3
#                                             eadyinstance = analyzeinstability(; eady_jac, Ri, θU, β, 𝕂ᵘ, 𝕂ᵇ, grid_t, hmt_scheme, hst_scheme, dissip_scheme, le, Nz)
#                                             push!(df, eadyinstance)
#                                         end
#                                     end
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

