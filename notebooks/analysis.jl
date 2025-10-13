### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 500f352c-6e16-11f0-215d-4f5a3075cb33
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	using Revise
	using CairoMakie
	using DataFrames
	import DrWatson: produce_or_load, @unpack, @dict, @strdict
	import GridOperatorAnalysis: bb, eady_jacobian, eady_background_flow, colpt_type
	import GridOperatorAnalysis: construct_flowΔxz, construct_flowkz
	import GridOperatorAnalysis: construct_lcc_sys, fourier_transform_sys, lowering_sys
	import GridOperatorAnalysis: VertexIndex, EdgeIndex, CellIndex
	import GridOperatorAnalysis: t, z, nS, evalat, dims, v, c, e
	import GridOperatorAnalysis: TriAFlow, TriBFlow, TriCFlow, HexCFlow
	import GridOperatorAnalysis: TriA, TriB, TriC, HexC
	import GridOperatorAnalysis
	using HypertextLiteral: @htl, @html_str
	using LaTeXStrings
	import LinearAlgebra: eigen, norm
	using PlutoUI
	using RuntimeGeneratedFunctions
	import Symbolics: substitute, @syms, Num, Differential, expand, taylor_coeff
	import Symbolics
	import SymbolicUtils
	using Unicode
	using UUIDs: uuid1
end;

# ╔═╡ e653beef-bd70-4f94-867c-b802f416408e
begin
	include(joinpath("..", "src", "simulations.jl"))
end;

# ╔═╡ 534dfd2f-f7e5-4253-83dc-0be9e94cba01
md"""
## Constants
"""

# ╔═╡ b02f29a5-7f20-49dd-aa0a-9b30d8b4a1ed
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	const f₀ = -1e-4
	const g  = 1e9
	const N² = 1e-6
	const U  = 1.0
	const H  = 4000
	const Nz = 4
end;
  ╠═╡ =#

# ╔═╡ 99dab097-3490-46cb-bf26-62b60b51d3c2
const Nz = 4

# ╔═╡ 2e485b46-7aa1-483d-9c25-95cd98439ccc
md"""
## Discretization
"""

# ╔═╡ 5870aa5c-cda6-496f-8836-c9e749e0c629
begin
	@syms vin₁::Int vin₂::Int vin₃::Int
	vin = VertexIndex((vin₁, vin₂, vin₃))
	@syms cin₁::Int cin₂::Int cin₃::Int
	cin = CellIndex((cin₁, cin₂, cin₃))
	@syms ein₁::Int ein₂::Int ein₃::Int
	ein = EdgeIndex((ein₁, ein₂, ein₃))
end;

# ╔═╡ a6e1eda0-c7b6-46be-b435-b999d08d83ba
md"""
__Choice of the grid__: $(@bind _grid_t PlutoUI.Select([:TriA => "triangular A-Grid", :TriC => "triangular C-Grid", :TriB => "triangular B-Grid", :HexC => "hexagonal C-Grid"]; default=:TriA))
"""

# ╔═╡ a3eecd80-eb24-47f6-ab77-1bb4b6e109bb
grid_t = Val(_grid_t)

# ╔═╡ b9fa1a9b-7989-4f1c-abe8-de153cf8085f
dflow = let
	if colpt_type(grid_t, :u⃗) == :edge
        @variables (du⃗(t,z))[-nS:nS,-nS:nS,1:dims(colpt_type(grid_t, :u⃗))]
    else 
        @variables (du⃗(t,z))[1:2,-nS:nS,-nS:nS,1:dims(colpt_type(grid_t, :u⃗))]
    end
    @variables begin
        (db(t,z))[-nS:nS,-nS:nS, 1:dims(colpt_type(grid_t, :b))]
        (dη(t))[-nS:nS, -nS:nS, 1:dims(colpt_type(grid_t, :η))]
        (∫∇ᵀdu⃗dz(t,z))[-nS:nS,-nS:nS, 1:dims(colpt_type(grid_t, :∫∇ᵀu⃗dz))]
        (dw(t,z))[-nS:nS, -nS:nS, 1:dims(colpt_type(grid_t, :w))]
        (dp(t,z))[-nS:nS, -nS:nS, 1:dims(colpt_type(grid_t, :p))]
    end
    construct_flowΔxz(grid_t, nS; u⃗ = du⃗, w = dw, ∫∇ᵀu⃗dz = ∫∇ᵀdu⃗dz, b = db, p = dp, η = dη)
end;

# ╔═╡ e9c5a7ef-6412-4b59-9aa7-d240277550e0
begin
	∂ₜ = Differential(t)
	∂₃ = Differential(z)
end

# ╔═╡ 6e5eb153-8d92-47f9-bee5-8f823abd11e9
function select_hmt_scheme(grid_t)
	schemes =
	if grid_t == :TriA
		[
			:standard => "Standard",
		]
	elseif grid_t == :TriB
		[
			:asc => "advective form, streamline derivative on cells",
			:avi => "advective form, vector-invariant", 
			:fdv => "flux form, divergence on vertices", 
			:fdcre => "flux form, diverence on cells with reconstruction on edges"
		]
	elseif grid_t == :TriC
		[
			:ICON => "ICON"
		]
	elseif grid_t == :HexC
		[
			:MPAS => "MPAS"
		]
	end
	PlutoUI.Select(schemes)
end

# ╔═╡ 826aa7f7-c0f1-4dd4-978f-dc093f970d84
function select_hst_scheme(grid_t)
	schemes =
	if grid_t == :TriA
		[
			:low => "low-order accurate",
		]
	elseif grid_t == :TriB
		[
			:low => "low-order accurate",
			:high => "high-order accurate",
		]
	elseif grid_t == :TriC
		[
			:low => "low-order accurate",
		]
	elseif grid_t == :HexC
		[
			:low => "low-order accurate",
		]
	end
	PlutoUI.Select(schemes)
end

# ╔═╡ b1a33624-23c1-4102-85b7-1e63980f3bf2
function EqSchemes(schemes)
	PlutoUI.combine() do Child
	@htl("""
	$([
		@htl("$(Child(string(i), Select(s))) +")	
		for (i, s) in enumerate(schemes[1:end-1])
	 ])
	$(Child(string(length(schemes)), Select(schemes[end]))) = 0 
	""")
	end
end

# ╔═╡ 29bef65a-0414-4304-a718-189b6885d2e3
schemes = Dict(
	:TriA => [
		[
			("∂ₜ(evalat(vout, vin, u⃗[iTH]))", 0) => "∂ₜ(evalat(vout, vin, u⃗[iTH]))"],
		[
			("TriA.u⃗ᵀ∇(vout, vin, u⃗, u⃗)[iTH]", 0) => "TriA.u⃗ᵀ∇(vout, vin, u⃗, u⃗)[iTH]", 
			("TriA.u⃗ᵀ∇(vout, vin, u⃗, u⃗)[iTH]", 1) => "evalat(vout, vin, ū⃗ᵀ∇u⃗[iTH])"
		],
		[
			("f₀ * evalat(vout, vin, u⃗⊥[iTH])", 0) => "f₀ * evalat(vout, vin, u⃗⊥[iTH])"
		],
		[
			("evalat(vout, vin, w) * ∂₃(evalat(vout, vin, u⃗[iTH]))", 0) => "evalat(vout, vin, w) * ∂₃(evalat(vout, vin, u⃗[iTH]))"
		],
		[
			("TriA.∇vv(vout, vin, p)[iTH]", 0) => "TriA.∇vv(vout, vin, p)[iTH]", 
			("TriA.∇vv(vout, vin, p)[iTH]", 1) => "evalat(vout, vin, ∇p[iTH])"
		],
		[
			("g * TriA.∇vv(vout, vin, η)[iTH]", 0) => "g * TriA.∇vv(vout, vin, η)[iTH]", 
			("g * TriA.∇vv(vout, vin, η)[iTH]", 1) => "g * evalat(vout, vin, ∇η[iTH])"
		],
		[
			("𝕂ᵘ * -TriA.Δ(vout, vin, u⃗[iTH])", 0) => "𝕂ᵘ * -TriA.Δ(vout, vin, u⃗[iTH])", 
			("𝕂ᵘ * TriA.Δ(vout, v, TriA.Δ(v, vin, u⃗[iTH]))", 0) => "𝕂ᵘ * TriA.Δ(vout, v, TriA.Δ(v, vin, u⃗[iTH]))", 
			("𝕂ᵘ * -TriA.Δ(vout, vin, u⃗[iTH])", 1) => "𝕂ᵘ * -evalat(vout, vin, Δu⃗[iTH])", 
			("𝕂ᵘ * TriA.Δ(vout, v, TriA.Δ(v, vin, u⃗[iTH]))", 1) => "𝕂ᵘ * evalat(vout, vin, Δ²u⃗[iTH])"
		],
	],
	:TriB => [
		[
			("∂ₜ(evalat(cout, cin, u⃗[iTH]))", 0) => "∂ₜ(evalat(cout, cin, u⃗[iTH]))"
		],
		[
			("TriB.u⃗ᵀ∇_asc(cout, cin, u⃗, u⃗)[iTH]", 0) => "TriB.u⃗ᵀ∇_asc(cout, cin, u⃗, u⃗)[iTH]", 
			("TriB.u⃗ᵀ∇_avi(cout, cin, u⃗, u⃗)[iTH]", 0) => "TriB.u⃗ᵀ∇_avi(cout, cin, u⃗, u⃗)[iTH]", 
			("TriB.u⃗ᵀ∇_fdv(cout, cin, u⃗, u⃗)[iTH]", 0) => "TriB.u⃗ᵀ∇_fdv(cout, cin, u⃗, u⃗)[iTH]", 
			("TriB.u⃗ᵀ∇_fdcre(cout, cin, u⃗, u⃗)[iTH]", 0) => "TriB.u⃗ᵀ∇_fdcre(cout, cin, u⃗, u⃗)[iTH]", 
			("TriB.u⃗ᵀ∇_asc(cout, cin, u⃗, u⃗)[iTH]", 1) => "evalat(vout, vin, u⃗ᵀ∇u⃗[iTH])"
		],
		[
			("f₀ * evalat(cout, cin, u⃗⊥[iTH])", 0) => "f₀ * evalat(cout, cin, u⃗⊥[iTH])"
		],
		[
			("TriB.av_cv(cout, vin, w) * ∂₃(evalat(cout, cin, u⃗[iTH]))", 0) => "TriB.av_cv(cout, vin, w) * ∂₃(evalat(cout, cin, u⃗[iTH]))",
			("TriB.av_cv(cout, vin, w) * ∂₃(evalat(cout, cin, u⃗[iTH]))", 1) => "evalat(cout, cin, w̄) * ∂₃(evalat(cout, cin, u⃗[iTH]))"
		],
		[
			("TriB.∇cv(cout, vin, p)[iTH]", 0) => "TriB.∇cv(cout, vin, p)[iTH]", 
			("TriB.∇cv(cout, vin, p)[iTH]", 1) => "evalat(cout, vin, ∇p[iTH])"
		],
		[
			("g * TriB.∇cv(cout, vin, η)[iTH]", 0) => "g * TriB.∇cv(cout, vin, η)[iTH]", 
			("g * TriB.∇cv(cout, vin, η)[iTH]", 1) => "g * evalat(cout, vin, ∇η[iTH])"
		],
		[
			("𝕂ᵘ * -TriB.Δ⃗(cout, cin, u⃗)[iTH]", 0) => "𝕂ᵘ * -TriB.Δ⃗(cout, cin, u⃗)[iTH]",
			("𝕂ᵘ * -TriB.Δ⃗(cout, cin, u⃗)[iTH]", 1) => "𝕂ᵘ * -evalat(cout, cin, Δu⃗[iTH])", 
			("𝕂ᵘ * TriB.Δ⃗(cout, c, TriB.Δ⃗(c, cin, u⃗))[iTH]", 0) => "𝕂ᵘ * TriB.Δ⃗(cout, c, TriB.Δ⃗(c, cin, u⃗))[iTH]",
			("𝕂ᵘ * TriB.Δ⃗(cout, c, TriB.Δ⃗(c, cin, u⃗))[iTH]", 1) => "𝕂ᵘ * evalat(cout, cin, Δ²u⃗)",
			
		],
	],
	:TriC => [
		[
			("∂ₜ(evalat(eout, ein, u⃗))", 0) => "∂ₜ(evalat(eout, ein, u⃗))"
		],
		[
			("TriC.u⃗ᵀ∇(eout, ein, u⃗, u⃗)", 0) => "TriC.u⃗ᵀ∇(eout, ein, u⃗, u⃗)",
			("TriC.u⃗ᵀ∇(eout, ein, u⃗, u⃗)", 1) => "evalat(eout, ein, u⃗ᵀ∇u⃗)"
		],
		[
			("TriC.ℳ̃(eout, ein, v, u⃗, f₀)", 0) => "TriC.ℳ̃(eout, ein, v, u⃗, f₀)"
		],
		[
			("TriC.Pᵀec(eout, cin, w * TriC.Pce(cin, ein, ∂₃(u⃗)))", 0) => "TriC.Pᵀec(eout, cin, w * TriC.Pce(cin, ein, ∂₃(u⃗)))", 
			("TriC.Pᵀec(eout, cin, w * TriC.Pce(cin, ein, ∂₃(u⃗)))", 1) => "evalat(eout, ein, w̄) * ∂₃(evalat(eout, ein, ū))"
		],
		[
			("TriC.ℳ(eout,  e, TriC.∇ec(e, cin, p))", 0) => "TriC.ℳ(eout,  e, TriC.∇ec(e, cin, p))"
		],
		[
			("g * TriC.ℳ(eout,  e, TriC.∇ec(e, cin, η))", 0) => "g * TriC.ℳ(eout,  e, TriC.∇ec(e, cin, η))"
		],
		[
			("𝕂ᵘ *- TriC.Δ⃗(eout, ein, u⃗)", 0) => "𝕂ᵘ *- TriC.Δ⃗(eout, ein, u⃗)",
			("𝕂ᵘ * TriC.Δ⃗(eout, e, TriC.Δ⃗(e, ein, u⃗))", 0) => "𝕂ᵘ * TriC.Δ⃗(eout, e, TriC.Δ⃗(e, ein, u⃗))", 
		],
	]
)

# ╔═╡ 92853482-fd36-4fe1-b515-9c7696b71c95
bschemes = Dict(
	:TriA => [
		[
			("∂ₜ(evalat(vout, vin, b))", 0) => "∂ₜ(evalat(vout, vin, b))"
		],
		[
			("TriA.u⃗∇ᵀ(vout, vin, u⃗, b)", 0) => "TriA.u⃗∇ᵀ(vout, vin, u⃗, b)",
			("TriA.u⃗∇ᵀ(vout, vin, u⃗, b)", 1) => "evalat(vout, vin, u⃗ᵀ∇b)",
			("TriA.u⃗∇ᵀ(vout, vin, u⃗, b-_le^2/8*TriA.Δ(vin,vin,b))", 0) => "TriA.u⃗∇ᵀ_high(vout, vin, u⃗, b)",
			("evalat(vout, vin, ū⃗ᵀ∇b)+evalat(vout, vin, u⃗ᵀ∇b̄)", 0) => "evalat(vout, vin, ū⃗ᵀ∇b)+evalat(vout, vin, u⃗ᵀ∇b̄)"
		],
		[
			("evalat(vout, vin, w) * ∂₃(evalat(vout, vin, b))", 0) => "evalat(vout, vin, w) * ∂₃(evalat(vout, vin, b))"
		],
		[
			("𝕂ᵇ * -TriA.Δ(vout, vin, b)", 0) => "𝕂ᵇ * -TriA.Δ(vout, vin, b)", 
			("𝕂ᵇ * -TriA.Δ(vout, vin, b)", 1) => "𝕂ᵇ * -evalat(vout, vin, Δb)",
			("𝕂ᵇ * TriA.Δ(vout, v, TriA.Δ(v, vin, b))", 0) => "𝕂ᵇ * TriA.Δ(vout, v, TriA.Δ(v, vin, b))",
			("𝕂ᵇ * TriA.Δ(vout, v, TriA.Δ(v, vin, b))", 1) => "𝕂ᵇ * evalat(vout, vin, Δ²b)"
		],
	],
	:TriB => [
		[
			("∂ₜ(evalat(vout, vin, b))", 0) => "∂ₜ(evalat(vout, vin, b))"
		],
		[
			("TriB.u⃗∇ᵀ(vout, cin, vin, u⃗, b; γ=3//4)", 0) => "TriB.u⃗∇ᵀ(vout, cin, vin, u⃗, b; γ=3//4)", 
		 	("TriB.u⃗∇ᵀ_low(vout, cin, vin, u⃗, b)", 0) => "TriB.u⃗∇ᵀ_low(vout, cin, vin, u⃗, b)", 
			("TriB.u⃗∇ᵀ_low(vout, cin, vin, u⃗, b)", 1) => "evalat(vout, vin, u⃗ᵀ∇b)"
		],
		[
			("evalat(vout, vin, w) * ∂₃(evalat(vout, vin, b))", 0) => "evalat(vout, vin, w) * ∂₃(evalat(vout, vin, b))"
		],
		[
			("𝕂ᵇ * -TriB.Δ(vout, vin, b)", 0) => "𝕂ᵇ * -TriB.Δ(vout, vin, b)",
			("𝕂ᵇ * -TriB.Δ(vout, vin, b)", 1) => "𝕂ᵇ * -evalat(vout, vin, Δb)",
			("𝕂ᵇ * TriB.Δ(vout, v, TriB.Δ(v, vin, b))", 0) => "𝕂ᵇ * TriB.Δ(vout, v, TriB.Δ(v, vin, b))",
			("𝕂ᵇ * TriB.Δ(vout, v, TriB.Δ(v, vin, b))", 1) => "𝕂ᵇ * evalat(vout, vin, Δ²b)", 
		],
	],
	#, , 
	:TriC => [
		[
			("∂ₜ(evalat(cout, cin, b))", 0) => "∂ₜ(evalat(cout, cin, b))"
		],
		[
			("TriC.u⃗∇ᵀ(cout, ein, cin, u⃗, b)", 0) => "TriC.u⃗∇ᵀ(cout, ein, cin, u⃗, b)",
			("TriC.u⃗∇ᵀ(cout, ein, cin, u⃗, b-_le^2/24*TriC.Δ(cin,cin,b))", 0) => "TriC.u⃗∇ᵀ_high(cout, ein, cin, u⃗, b)"
		],
		[
			("evalat(cout, cin, w) * ∂₃(evalat(cout, cin, b))", 0) => "evalat(cout, cin, w) * ∂₃(evalat(cout, cin, b))"
		],
		[
			("𝕂ᵇ * -TriC.Δ(cout, c, TriC.Δ(c, cin, b))", 0) => "𝕂ᵇ * -TriC.Δ(cout, c, TriC.Δ(c, cin, b))", 
			("𝕂ᵇ * TriC.Δ(cout, cin, b)", 0) => "𝕂ᵇ * TriC.Δ(cout, cin, b)"
		]
	],
)

# ╔═╡ 86ab6145-af2f-4e0d-90f9-4f5492f05ee1
ηschemes = Dict(
	:TriA => [
		[
			("∂ₜ(evalat(vout, vin, η))", 0) => "∂ₜ(evalat(vout, vin, η))"
		],
		[
			("evalat(vout, vin, ∫∇ᵀu⃗dz)", 0) => "evalat(vout, vin, ∫∇ᵀu⃗dz)"
		]
	],
	:TriB => [
		[
			("∂ₜ(evalat(vout, vin, η))", 0) => "∂ₜ(evalat(vout, vin, η))"
		],
		[
			("evalat(vout, vin, ∫∇ᵀu⃗dz)", 0) => "evalat(vout, vin, ∫∇ᵀu⃗dz)"
		]
	],
	:TriC => [
		[
			("∂ₜ(evalat(cout, cin, η))", 0) => "∂ₜ(evalat(cout, cin, η))"
		],
		[
			("evalat(cout, cin, ∫∇ᵀu⃗dz)", 0) => "evalat(cout, cin, ∫∇ᵀu⃗dz)"
		]
	],
)

# ╔═╡ 6c58a2e9-7b61-415e-b020-20b7170082e2
cschemes = Dict(
	:TriA => [
		[
			("TriA.∇ᵀvv(vout, vin, u⃗)", 0) => "TriA.∇ᵀvv(vout, vin, u⃗)", 
			("TriA.∇ᵀvv(vout, vin, u⃗)", 1) => "evalat(vout, vin, ∇ᵀu⃗)"
		],
		[
			("∂₃(evalat(vout, vin, w))", 0) => "∂₃(evalat(vout, vin, w))"
		]
	],
	:TriB => [
		[
			("TriB.∇ᵀvc(vout, cin, u⃗)", 0) => "TriB.∇ᵀvc(vout, cin, u⃗)", 
			("TriB.∇ᵀvc(vout, cin, u⃗)", 1) => "evalat(vout, vin, ∇ᵀu⃗)"
		],
		[
			("∂₃(evalat(vout, vin, w))", 0) => "∂₃(evalat(vout, vin, w))"
		]
	],
	:TriC => [
		[
			("TriC.∇ᵀce(cout, e, TriC.ℳ(e, ein, u⃗))", 0) => "TriC.∇ᵀce(cout, e, TriC.ℳ(e, ein, u⃗))", 
			("TriC.∇ᵀce(cout, e, TriC.ℳ(e, ein, u⃗))", 1) => "evalat(cout, cin, ∇ᵀu⃗)"
		],
		[
			("∂₃(evalat(cout, cin, w))", 0) => "∂₃(evalat(cout, cin, w))"
		]
	],
)

# ╔═╡ f4fe6284-3ea3-4e57-af2b-664fda7066a3
pschemes = Dict(
	:TriA => [
		[
			("∂₃(evalat(vout, vin, p))", 0) => "∂₃(evalat(vout, vin, p))"
		],
		[
			("-evalat(vout, vin, b)", 0) => "-evalat(vout, vin, b)"
		]
	],
	:TriB => [
		[
			("∂₃(evalat(vout, vin, p))", 0) => "∂₃(evalat(vout, vin, p))"
		],
		[
			("-evalat(vout, vin, b)", 0) => "-evalat(vout, vin, b)"
		]
	],
	:TriC => [
		[
			("∂₃(evalat(cout, cin, p))", 0) => "∂₃(evalat(cout, cin, p))"
		],
		[
			("-evalat(cout, cin, b)", 0) => "-evalat(cout, cin, b)"
		]
	]
)

# ╔═╡ 51d642b5-b044-4c2a-b7e0-aab048231544
WideCell(md"""
##### Horizontal Momentum Transport Equation
$(@bind vterms EqSchemes(schemes[_grid_t]))
"""; max_width=1500)

# ╔═╡ bbf8aa4b-0b66-40b0-b5ab-56f2b4baa641
WideCell(
md"""
__Buoyancy transport balance__:

$(@bind bterms EqSchemes(bschemes[_grid_t]))
"""
; max_width=1500)

# ╔═╡ cffd4061-d02d-4e42-bbca-cc990927824b
WideCell(
md"""
__Surface elevation equation__:
	
$(@bind ηterms EqSchemes(ηschemes[_grid_t]))
"""
; max_width=1500)

# ╔═╡ c690b0a4-6130-4320-869f-2764fb681f3f
WideCell(
md"""
__Continuity Equation__:
	
$(@bind wterms EqSchemes(cschemes[_grid_t]))
"""
; max_width=1500)

# ╔═╡ b4024a03-2a30-4756-aebe-3f2e776cd8d6
WideCell(
md"""
__Hydrostatic Balance__:
	
$(@bind pterms EqSchemes(pschemes[_grid_t]))
"""
; max_width=1500)

# ╔═╡ 42b320c0-d9a3-4a84-bcc6-3498d71ceec8
md"""
__Discretization scheme for the momentum transport__: 
$(@bind hmt_scheme select_hmt_scheme(_grid_t))
"""

# ╔═╡ 07e1fa2c-7ef0-4a22-9c75-6640a626fe2c
md"""
__Discretization scheme for the momentum transport__: 
$(@bind hst_scheme select_hst_scheme(_grid_t))
"""

# ╔═╡ 694e01ee-5b10-44e4-a477-e8593ededa8c
md"""
__Dissipation scheme__: $(@bind dissip_scheme PlutoUI.Select([:harmonic => "harmonic", :biharmonic => "biharmonic"]; default=:biharmonic))
"""

# ╔═╡ 51c00889-ad8c-400d-b764-a4261e5a76e1
md"""
## Fourier Transform
"""

# ╔═╡ 4867203f-402e-4d3a-8ad0-ef0ba8478ed1
fflow = let
	if colpt_type(grid_t, :u⃗) == :edge
        @variables (u⃗̂(t,z))[1:dims(colpt_type(grid_t, :u⃗))]
    else
        @variables (u⃗̂(t,z))[1:2, 1:dims(colpt_type(grid_t, :u⃗))]
    end
    @variables begin
        (b̂(t,z))[1:dims(colpt_type(grid_t, :b))]
        (η̂(t))[1:dims(colpt_type(grid_t, :η))]
        (∫∇ᵀu⃗̂dz(t))[1:dims(colpt_type(grid_t, :∫∇ᵀu⃗dz))]
        (ŵ(t,z))[1:dims(colpt_type(grid_t, :w))]
        (p̂(t,z))[1:dims(colpt_type(grid_t, :p))]
    end
    construct_flowkz(grid_t, Num; u⃗=u⃗̂, w=ŵ, ∫∇ᵀu⃗dz=∫∇ᵀu⃗̂dz, b=b̂, p=p̂, η=η̂)
end;

# ╔═╡ cf11e62e-ee93-4b2f-9cdc-2cb228ca040c
md"""
## Instability Analysis
"""

# ╔═╡ 28470719-44cf-4b6a-a515-29f7f6a4acb2
config = @dict(grid_t, hmt_scheme, dissip_scheme, Nz, H);

# ╔═╡ 7c1d69bf-8612-4072-95b5-27b57e59b2ac
path = joinpath(@__DIR__, "..", "data")

# ╔═╡ 944eb1fc-3207-48fe-bec3-97f73e7a523a
# ╠═╡ disabled = true
#=╠═╡
begin
	data, file = produce_or_load(config, path) do config
		@unpack grid_t, hmt_scheme, Nz, H = config
		eady_jac_ex, eady_sys = eady_jacobian(Val(grid_t); Nz, H, U, hmt_scheme, dissip_scheme)
		@strdict(eady_jac_ex, eady_sys)
	end
	@unpack eady_jac_ex, eady_sys = data
end;
  ╠═╡ =#

# ╔═╡ 1aa38ca9-ab99-4d2b-ad5d-c5a18d055229
Symbolics.@variables _Ri _le _g _f₀ _N² _𝕂ᵘ _𝕂ᵇ _θU _β k l

# ╔═╡ a832b824-afdb-458c-9f91-ea023c52852e
bflow = eady_background_flow(grid_t, _le; f₀=_f₀, N²=_N², Ri=_Ri, H, U=1, θU=_θU, β=0);

# ╔═╡ 0ef5eb5f-6c6c-4ecc-87bd-16247f9056e4
begin
	@variables ϵ
	pflow = let
		state = bflow.state .+ ϵ .* dflow.state
		pflow = typeof(dflow)(state)
	end
end;

# ╔═╡ 0d0a0e92-2d63-4745-9114-8a10c5be58de
begin
	@variables Im
	upin = grid_t == Val(:TriA) ? vin : (grid_t == Val(:TriB) ? cin : ein)
	spin = grid_t == Val(:TriA) ? vin : (grid_t == Val(:TriB) ? vin : cin)
	u⃗ = if grid_t == Val(:TriC)
		pflow.u⃗[upin[1], upin[2], upin[3]]
	else
		[pflow.u⃗[iTH, upin[1], upin[2], upin[3]] for iTH=1:2]
	end
	b = pflow.b[spin[1], spin[2], spin[3]]
	η = pflow.η[spin[1], spin[2], spin[3]]
	w = pflow.w[spin[1], spin[2], spin[3]]
	p = pflow.p[spin[1], spin[2], spin[3]]
	∫∇ᵀu⃗dz = pflow.∫∇ᵀu⃗dz[spin[1], spin[2], spin[3]]
	u⃗ᵀ∇b̄ = let
		_u⃗  = ϵ * [dflow.u⃗[iTH, vin[1], vin[2], vin[3]] for iTH=1:2]
		∇b̄ = TriA.∇vv(VertexIndex((0,0,1)), vin, bflow.b[vin[1], vin[2], vin[3]])
		∇b̄ = substitute.(∇b̄, Ref(Dict(GridOperatorAnalysis.le => _le)))
		∇b̄ = Symbolics.simplify.(∇b̄; expand=true)
		@show ∇b̄
		_u⃗' * ∇b̄
	end
	ū⃗ᵀ∇b = let
		ū⃗  = [bflow.u⃗[iTH, vin[1], vin[2], vin[3]] for iTH=1:2]
		∇b = ϵ * Im * [k; l] * fflow.b[1]
		ū⃗' * ∇b
	end
end;

# ╔═╡ 61150eb4-c9ce-4fd8-8314-780b9c1c29b6
u⃗⊥ = [-u⃗[2]; u⃗[1]];

# ╔═╡ 0ad8eba1-f992-4f24-aa90-38ca629f8c72
ϕ = let
	ϕ = GridOperatorAnalysis.compute_phases(k, l, _le)
	#for (k, v) in pairs(ϕ)
	#	newv = Symbolics.simplify.(v; rewriter=rtrig)
	#	newv = Symbolics.expand.(newv)
	#	newv = Symbolics.substitute.(newv, Ref(sqrt3subs))
	#	ϕ[k] = Symbolics.simplify.(newv)
	#end
	ϕ
end;

# ╔═╡ 3588e3e2-193d-4c1c-9dbf-c046803cface
rtrig = let
	function p(x)
		 !isequal(x, _θU)
	end
	rcos = let
		x = Symbolics.variable(:x)
		cosx = Symbolics.taylor(cos(x),x,0:10; rationalize=false)
		Symbolics.@rule cos(~x::p) => substitute(cosx, Dict([x=>~x]))
	end
	rsin = let
		x = Symbolics.variable(:x)
		sinx = Symbolics.taylor(sin(x),x,0,0:10; rationalize=false)
		Symbolics.@rule sin(~x::p)=>substitute(sinx, Dict([x=>~x]))
	end
	SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(SymbolicUtils.RestartedChain([rcos, rsin])))
end

# ╔═╡ 9bc2c550-4008-42fc-9ac0-607d49f8f319
# ╠═╡ disabled = true
#=╠═╡
jac = let
	#jac = eady_jacobian_wo_mtk(Val(grid_t), _le; g, f₀, N², Ri=_Ri, Nz=32, 𝕂ᵘ=_𝕂ᵘ, 𝕂ᵇ=_𝕂ᵇ, θU=θU, β=β, hmt_scheme, hst_scheme, dissip_scheme)
	@showtime jac = GridOperatorAnalysis.eady_jacobian(Val(grid_t), k, l, _le; ϕ, g=_g, f₀=_f₀, Ri=_Ri, N²=_N², 𝕂ᵘ=_𝕂ᵘ, 𝕂ᵇ=_𝕂ᵇ, Nz=16, θU, β, hmt_scheme, hst_scheme, dissip_scheme)
	#jac = substitute.(jac, Ref(Dict(_g => g, _f₀ => f₀, _N² => N²)))
	#@showtime jac = Symbolics.expand.(jac)
	#@showtime Symbolics.simplify_fractions.(jac)
	#@showtime jac = Symbolics.simplify.(jac; rewriter=rtrig)
	#@showtime jac = substitute.(jac, Ref(sqrt3subs))
	#@showtime jac = Symbolics.simplify_fractions.(jac)	
end;
  ╠═╡ =#

# ╔═╡ e3bc1fcc-fb96-404f-8204-675174b3afe1
md"""
## Plotting
"""

# ╔═╡ 70f46342-4b5c-45ce-9bf9-13a44fd780fc
import CSV

# ╔═╡ be945a73-9afb-4b9a-ac78-4514fbb0e33e
df = let
	nβs  = [0, 5]
	nθUs = [0, 2]
	les = [6.25e3, 12.5e3]
	𝕍ᵘs = [0.0, 0.005]
	𝕍ᵇs = [0.0, 0.005]
	savepath = joinpath(@__DIR__, "..", "data", "allsims.csv")
	if isfile(savepath)
		df = CSV.read(savepath, DataFrame)
		select(df, Not(:Ks, :iωs), :Ks => ByRow(x->eval(Meta.parse(x))) => :Ks, :iωs => ByRow(x-> eval(Meta.parse(x))) => :iωs)
	else
		df = testall(; nθUs, nβs, les, 𝕍ᵘs, 𝕍ᵇs)
		CSV.write(savepath, df)
		df
	end
end

# ╔═╡ 3adadede-5704-441f-9756-08f0f820c723
md"""
### Baroclinic Axis
"""

# ╔═╡ 2bfe2cd8-cf70-4efa-a0d2-5a3fe447e9e8
md"""
#### Without Dissipation
"""

# ╔═╡ a371a682-d689-4147-8c8f-49a65a29cd8e
md"""
#### With Biharmonic Dissipation
"""

# ╔═╡ ffc466bf-7d98-4de1-9c12-0e8e2d559ed5
md"""
### Symmetric Axis
"""

# ╔═╡ 9790e852-5292-4ef1-9c11-a49344d423af
md"""
#### Without Dissipation
"""

# ╔═╡ 1effe5fd-e458-4bc3-8318-88f54396d32e
md"""
#### With Dissipation
"""

# ╔═╡ 6a193b8c-8eac-422a-a6b0-cb3962a44c39
md"""
#### Alignment of the Flow with the Mesh
"""

# ╔═╡ 96d447fc-0e4e-4c2d-8c98-fe50547906e8
md"""
##### Baroclinic Axis
"""

# ╔═╡ 72cdb1b1-0ef2-468b-9e25-af3732e00ee8
md"""
##### Symmetric Axis
"""

# ╔═╡ b861da6e-16f4-4799-bf5c-fca588b523b6
md"""
#### B-Grid on Triangular Mesh
"""

# ╔═╡ 5f414b9f-342d-4e80-8b72-ebc306b0cc78
html"""<hr>"""

# ╔═╡ 6e93f103-3988-4a68-8804-eabc14f12f5d
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ 5d57177f-c6c5-496f-942f-755b25fa957d
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;

# ╔═╡ c29f6ae0-e69b-4ad8-8615-8dd9a4aee408
floataside(md"""
Flow direction: $(@bind θU PlutoUI.Slider([0.0, π/12, π/6]; show_value=true))

Flow shift (in -HM²/f₀): $(@bind β PlutoUI.Slider(-0.5:0.1:0.5; default=0.0, show_value=true))

le: $(@bind le PlutoUI.Slider([1e-28, 1.575e3, 3.125e3, 6.25e3, 12.5e3, 25e3]; default=6.25e3, show_value=true))

Ri: $(@bind Ri Select([100//1, 1//2]; default=100//1))

Vᵘ: $(@bind Vᵘ PlutoUI.Slider([0, 1e-3, 5e-3, 1e-2]; show_value=true))

Vᵇ: $(@bind Vᵇ PlutoUI.Slider([0, 1e-3, 5e-3, 1e-2]; show_value=true))
""")

# ╔═╡ 4767f71a-fd13-4c26-8ec4-bf16ac6028a4
if dissip_scheme == :biharmonic
	𝕂ᵘ = Vᵘ * le^3
	𝕂ᵇ = Vᵇ * le^3
else
	𝕂ᵘ = Vᵘ * le
	𝕂ᵇ = Vᵇ * le
end

# ╔═╡ 69f7116f-5062-43b3-a004-acd5de35ed1e
let
	Ri     = 100
	size   = Ri > 1 ? (1500, 550) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/le)
		#min(1e-2, π/le * 2/√3 / norm(inv([1;0;;0.5;√3/2]) * [cos(θU); sin(θU)], 1))
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 3.0 : 2.5
	limits = Ri > 1 ? (0.0, 1.0, -0.02, 0.38) : (0.0, 1.1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :β=>x->x.==β, :θU=>x->x.==θU, :le=>x->x.==le, :𝕂ᵘ=>x->x.≈0.0, :𝕂ᵇ=>x->x.≈0.0, :Ri=>x->x.==100, :dissip_scheme=>x->Symbol.(x).==:harmonic)
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme) = row
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t))"* (grid_t == "TriA" ? "" : ":$(String(hmt_scheme))"), linewidth=3)
	end
	#axislegend()
	axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=36);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 953ca291-388d-4689-a891-fcdd62b659f7
let
	Ri     = 100
	size   = Ri > 1 ? (1400, 700) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/le)
		#min(1e-2, π/le * 2/√3 / norm(inv([1;0;;0.5;√3/2]) * [cos(θU); sin(θU)], 1))
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, 1, -0.1, 0.4) : (0.0, 1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :β=>x->x.==β, :θU=>x->x.==θU, :le=>x->x.==le, :𝕂ᵘ=>x->x.≈Vᵘ*le^3, :𝕂ᵇ=>x->x.≈Vᵇ*le^3, :Ri=>x->x.==100, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme) = row
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t)):$(String(hmt_scheme))", linewidth=3)
	end
	#axislegend()
	f[1,2] = Legend(f, ax; merge=true, valign=:top);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 93f798b4-b999-48e2-9cc4-2371d238ee9d
let
	Ri = 1/2
	size   = Ri > 1 ? (1400, 700) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/le)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, 1, -0.1, 0.4) : (0.0, 1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :β=>x->x.==β, :θU=>x->x.==θU, :le=>x->x.==le, :𝕂ᵘ=>x->x.≈0.0, :𝕂ᵇ=>x->x.≈0.0, :Ri=>x->x.==Ri, :dissip_scheme=>x->Symbol.(x).==:harmonic)
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme) = row
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t)):$(String(hmt_scheme))", linewidth=3)
	end
	#axislegend()
	f[1,2] = Legend(f, ax; merge=true, valign=:top);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ c68fde99-fd3d-454e-b8d3-c510b35ac02b
let
	Ri = 1/2
	size   = Ri > 1 ? (1400, 700) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/le)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, 1, -0.1, 0.4) : (0.0, 1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :β=>x->x.==β, :θU=>x->x.==θU, :le=>x->x.==le, :𝕂ᵘ=>x->x.≈Vᵘ*le^3, :𝕂ᵇ=>x->x.≈Vᵇ*le^3, :Ri=>x->x.==Ri, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme) = row
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t)):$(String(hmt_scheme))", linewidth=3)
	end
	#axislegend()
	f[1,2] = Legend(f, ax; merge=true, valign=:top);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 8d013161-0743-4907-b133-1f2bacbcc1d1
let
	Ri = 100
	size   = Ri > 1 ? (1400, 2*700) : (1400, 2*500)
	fₛ     = min(1e-2, 2/√3*π/le)
		#min(1e-2, π/le * 2/√3 / norm(inv([1;0;;0.5;√3/2]) * [cos(θU); sin(θU)], 1))
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, 1, -0.1, 0.4) : (0.0, 1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	for (i, β) in enumerate([0, 0.5])
		ax = Axis(f[i,1];
				 xlabel = "wavenumber / fₛ",
				 ylabel = "growth rate / N M⁻²",
				 xticks,
				 aspect,
				 limits,
				 )
		subdf = subset(df, :β=>x->x.≈β, :θU=>x->x.≈θU, :le=>x->x.==le, :𝕂ᵘ=>x->x.≈0.0, :𝕂ᵇ=>x->x.≈0.0, :Ri=>x->x.==Ri, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
		for row in eachrow(subdf)
			(; Ks, iωs, grid_t, hmt_scheme) = row
			lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t)):$(String(hmt_scheme))", linewidth=3)
		end
		f[i,2] = Legend(f, ax; merge=true, valign=:top)
		#Label(f[i,1, TopLeft()], "($i)", )
	end
	#axislegend()
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 10f7c553-6875-4d22-be26-58947e8a381d
let
	Ri = 1/2
	size   = Ri > 1 ? (1400, 2*700) : (1400, 2*500)
	fₛ     = min(1e-2, 2/√3*π/le)
		#min(1e-2, π/le * 2/√3 / norm(inv([1;0;;0.5;√3/2]) * [cos(θU); sin(θU)], 1))
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, 1, -0.1, 0.4) : (0.0, 1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	for (i, θU) in enumerate([0, π/6])
		ax = Axis(f[i,1];
				 xlabel = "wavenumber / fₛ",
				 ylabel = "growth rate / N M⁻²",
				 xticks,
				 aspect,
				 limits,
				 )
		subdf = subset(df, :β=>x->x.==β, :θU=>x->x.≈θU, :le=>x->x.==le, :𝕂ᵘ=>x->x.≈0.0, :𝕂ᵇ=>x->x.≈0.0, :Ri=>x->x.==Ri, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
		for row in eachrow(subdf)
			(; Ks, iωs, grid_t, hmt_scheme) = row
			lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t)):$(String(hmt_scheme))", linewidth=3)
		end
		f[i,2] = Legend(f, ax; merge=true, valign=:top)
		#Label(f[i,1, TopLeft()], "($i)", )
	end
	#axislegend()
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 2cdf19d9-870e-495b-bba9-2c1e29f28ba6
let
	size   = Ri > 1 ? (1400, 700) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/le)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, 1, -0.1, 0.4) : (0.0, 1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	if dissip_scheme == :biharmonic
		𝕂ᵘ = Vᵘ * le^3
		𝕂ᵇ = Vᵇ * le^3
	else
		𝕂ᵘ = Vᵘ * le
		𝕂ᵇ = Vᵇ * le
	end
	subdf = subset(df, :β=>x->x.==β, :θU=>x->x.==θU, :le=>x->x.==le, :𝕂ᵘ=>x->x.==𝕂ᵘ, :𝕂ᵇ=>x->x.==𝕂ᵇ, :Ri=>x->x.==Ri, :grid_t=>x->x.=="TriB")
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme) = row
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t)):$(String(hmt_scheme))", linewidth=3)
	end
	#axislegend()
	f[1,2] = Legend(f, ax; merge=true, valign=:top);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ abcf94dc-cfac-47e6-b903-7fbf3a63a7aa
sqrt3subs = Dict(
	GridOperatorAnalysis.sqrt3^2=>3, GridOperatorAnalysis.sqrt3^3=>3*GridOperatorAnalysis.sqrt3, GridOperatorAnalysis.sqrt3^4=>9, GridOperatorAnalysis.sqrt3^5=>9*GridOperatorAnalysis.sqrt3, GridOperatorAnalysis.sqrt3^6=>27,
	GridOperatorAnalysis.sqrt3^7=>27*GridOperatorAnalysis.sqrt3,
	GridOperatorAnalysis.sqrt3^8=>81,
);

# ╔═╡ c9ac2630-ea0e-4048-a161-8c9f925d1d20
let
	mdt = md"""
```math
[∂ₜ(evalat(vout, vin, u⃗[iTH])), hmt(vout, vin, u⃗, u⃗)[iTH], f₀ * evalat(vout, vin, u⃗⊥[iTH]), evalat(vout, vin, w) * ∂₃(evalat(vout, vin, u⃗[iTH])), TriA.∇vv(vout, vin, p)[iTH], g * TriA.∇vv(vout, vin, η)[iTH], -Dᵘₕ[iTH]]
```
"""
	for (i,c) in enumerate(mdt.content)
		if typeof(c) == Markdown.LaTeX
			c.formula = replace(c.formula, r"(.)⃗"=>s"\\vec{\1}")
			c.formula = replace(c.formula, r"vout"=>s"v_{\\text{out}}")
			c.formula = replace(c.formula, r"vin"=>s"v_{\\text{in}}")
			c.formula = replace(c.formula, r"iTH"=>s"i_{\\text{TH}}")
			c.formula = replace(c.formula, r"evalat"=>s"\\mathrm{evalat}")
			c.formula = replace(c.formula, r"hmt"=>s"\\mathrm{hmt}")
			c.formula = replace(c.formula, r"\*"=>s"\\cdot")
		end
	end
	mdt
end

# ╔═╡ d35d5d13-d7a2-44ee-8855-be4aa31edf4c
function assemble_sys!(sys, vterms, bterms, ηterms, wterms, pterms; 𝕂ᵘ, 𝕂ᵇ)
	(; state, momentum_transport_eq, hydrostatic_balance_eq, continuity_eq, buoyancy_transport_eq, surface_elevation_eq) = sys
	# horizontal momentum transport equation
	for iH=1:dims(colpt_type(grid_t, :momentum_transport_eq))
		if colpt_type(grid_t, :u⃗) == :edge
			ts = []
			for t in collect(first.(values(vterms)))
				t = replace(t, "vout" => "VertexIndex((0,0,$iH))", "cout" => "CellIndex((0,0,$iH))", "eout" => "EdgeIndex((0,0,$iH))", "𝕂ᵘ" => "$𝕂ᵘ")
				et = @eval $(Meta.parse(t))
				push!(ts, et)
			end
			push!(momentum_transport_eq[iH], ts...)
		else
			for iTH=1:2
				ts = []
				for t in collect(first.(values(vterms)))
					t = replace(t, "iTH" => "$iTH", "vout" => "VertexIndex((0,0,$iH))", "cout" => "CellIndex((0,0,$iH))", "eout" => "EdgeIndex((0,0,$iH))", "𝕂ᵘ" => "$𝕂ᵘ")
					et = @eval $(Meta.parse(t))
					push!(ts, et)
				end
				push!(momentum_transport_eq[iTH, iH], ts...)
			end
		end
	end
	# horizontal buoyancy transport equation
	for iH=1:dims(colpt_type(grid_t, :buoyancy_transport_eq))
		ts = []
		for t in collect(first.(values(bterms)))
			t = replace(t, "vout" => "VertexIndex((0,0,$iH))", "cout" => "CellIndex((0,0,$iH))", "eout" => "EdgeIndex((0,0,$iH))", "𝕂ᵇ" => "$𝕂ᵇ")
			et = @eval $(Meta.parse(t))
			push!(ts, et)
		end
		push!(buoyancy_transport_eq[iH], ts...)
	end
	# surface elevation equation
	for iH=1:dims(colpt_type(grid_t, :surface_elevation_eq))
		ts = []
		for t in collect(first.(values(ηterms)))
			t = replace(t, "vout" => "VertexIndex((0,0,$iH))", "cout" => "CellIndex((0,0,$iH))", "eout" => "EdgeIndex((0,0,$iH))")
			et = @eval $(Meta.parse(t))
			push!(ts, et)
		end
		push!(surface_elevation_eq[iH], ts...)
	end
	# continuity equation
	for iH=1:dims(colpt_type(grid_t, :continuity_eq))
		ts = []
		for t in collect(first.(values(wterms)))
			t = replace(t, "vout" => "VertexIndex((0,0,$iH))", "cout" => "CellIndex((0,0,$iH))", "eout" => "EdgeIndex((0,0,$iH))")
			et = @eval $(Meta.parse(t))
			push!(ts, et)
		end
		push!(sys.continuity_eq[iH], ts...)
	end
	# hydrostatic balance equation
	for iH=1:dims(colpt_type(grid_t, :hydrostatic_balance_eq))
		ts = []
		for t in collect(first.(values(pterms)))
			t = replace(t, "vout" => "VertexIndex((0,0,$iH))", "cout" => "CellIndex((0,0,$iH))", "eout" => "EdgeIndex((0,0,$iH))")
			et = @eval $(Meta.parse(t))
			push!(ts, et)
		end
		push!(sys.hydrostatic_balance_eq[iH], ts...)
	end

	for eq in sys.state
		eq .= substitute.(eq, Ref(Dict(GridOperatorAnalysis.sqrt3^2=>3//1, GridOperatorAnalysis.le=>_le)))
		eq .= [Symbolics.expand(taylor_coeff(Symbolics.expand(expr), ϵ, 1)) for expr in eq]
	end
end

# ╔═╡ f0f59948-54d8-477e-b2cb-595f68ca13b2
begin
	sys = construct_lcc_sys(grid_t, Vector{Num})
	assemble_sys!(sys, vterms, bterms, ηterms, wterms, pterms; 𝕂ᵘ=_𝕂ᵘ, 𝕂ᵇ=_𝕂ᵇ)
end;

# ╔═╡ 3e557bf0-c8aa-49d5-8860-03465dc06731
function lowestorder(expr)
	expr = expand(expr)
	lowexpr = 0
	ps = [(0,0),(1,0), (0,1), (2,0), (1,1), (0,2)]
	for (pk, pl) in ps
		kexpr  = taylor_coeff(expr, k, pk)
		klexpr = taylor_coeff(kexpr, l, pl)
		lowexpr += klexpr * k^pk * l^pl
	end
	lowexpr
end

# ╔═╡ 3ff311ec-aaaf-42f1-b11a-25be4632192c
function exactop(iH, cp_t_out, expr)
	expr = substitute(expr, Dict(GridOperatorAnalysis.sqrt3^2=>3//1, GridOperatorAnalysis.le=>_le))
	expr = Symbolics.expand(taylor_coeff(Symbolics.expand(expr), ϵ, 1))
	fexpr = GridOperatorAnalysis.fourier_transform_expression(iH, cp_t_out, expr; dflow, fflow, ϕ)
	fexpr = Symbolics.simplify(fexpr; expand=true, rewriter=rtrig)
	fexpr = substitute(fexpr, Dict(_le => 1e-20, GridOperatorAnalysis.sqrt3 => √3))
	Symbolics.simplify(fexpr; expand=true)
end

# ╔═╡ 4efa8fb3-a6a9-4e73-8d61-8c7c45f0d35f
function exactop(fexpr)
	fexpr = Symbolics.simplify(fexpr; expand=true, rewriter=rtrig)
	fexpr = Symbolics.simplify(fexpr)
	fexpr = lowestorder(fexpr)
	fexpr = substitute(fexpr, Dict(_le => 1e-20, GridOperatorAnalysis.sqrt3 => √3))
	Symbolics.simplify(fexpr; expand=true)
end

# ╔═╡ 5ad386c4-1dfa-4c20-b7ab-e6a2d13ae1ce
begin
	fsys = fourier_transform_sys(grid_t, sys; dflow, fflow, ϕ, subs=Dict{Any, Any}(Im=>im))
	# horizontal momentum transport
	for (i, islimit) in enumerate(collect(last.(values(vterms))))
		momentum_transport_eq = fsys.momentum_transport_eq
		for iH=1:dims(colpt_type(grid_t, :momentum_transport_eq))
			if grid_t == Val(:TriC)
				if islimit == 1
					expr = momentum_transport_eq[iH][i]
					momentum_transport_eq[iH][i] = exactop(expr)
				end
			else
				for iTH=1:2
					if islimit == 1
						expr = momentum_transport_eq[iTH, iH][i]
						momentum_transport_eq[iTH, iH][i] = exactop(expr)
					end
				end
			end
		end
	end
	# horizontal buoyancy transport
	for (i, islimit) in enumerate(collect(last.(values(bterms))))
		buoyancy_transport_eq = fsys.buoyancy_transport_eq
		for iH=1:dims(colpt_type(grid_t, :buoyancy_transport_eq))
			if islimit == 1
				expr = buoyancy_transport_eq[iH][i]
				buoyancy_transport_eq[iH][i] = exactop(expr)
			end
		end
	end
	# surface elevation
	for (i, islimit) in enumerate(collect(last.(values(ηterms))))
		surface_elevation_eq = fsys.surface_elevation_eq
		for iH=1:dims(colpt_type(grid_t, :surface_elevation_eq))
			if islimit == 1
				expr = surface_elevation_eq[iH][i]
				surface_elevation_eq[iH][i] = exactop(expr)
			end
		end
	end
	# continuity equation
	for (i, islimit) in enumerate(collect(last.(values(wterms))))
		continuity_eq = fsys.continuity_eq
		for iH=1:dims(colpt_type(grid_t, :continuity_eq))
			if islimit == 1
				expr = continuity_eq[iH][i]
				continuity_eq[iH][i] = exactop(expr)
			end
		end
	end
	# hydrostatic balance equation
	for (i, islimit) in enumerate(collect(last.(values(pterms))))
		hydrostatic_balance_eq = fsys.hydrostatic_balance_eq
		for iH=1:dims(colpt_type(grid_t, :hydrostatic_balance_eq))
			if islimit == 1
				expr = hydrostatic_balance_eq[iH][i]
				hydrostatic_balance_eq[iH][i] = exactop(expr)
			end
		end
	end
end;

# ╔═╡ a77404c3-511d-4ae6-9e6c-2da0cf490ee9
# ╠═╡ disabled = true
#=╠═╡
fsys.momentum_transport_eq[2,1][7] = (k^2 + l^2) * fflow.u⃗[2,1] * _𝕂ᵘ
  ╠═╡ =#

# ╔═╡ a2687cd8-4e59-414a-8560-e8293935e6b0
jac = lowering_sys(grid_t, fsys, fflow; Nz=16);

# ╔═╡ c6332760-298c-42f0-ba13-0936cab3fdcd
fun = let
	Symbolics.build_function.(substitute.(jac, Ref(Dict(GridOperatorAnalysis.sqrt3=>√3))), k, l, _Ri, _le, _f₀, _g, _N², _𝕂ᵘ, _𝕂ᵇ, _θU; expression=Val{false})
end;

# ╔═╡ ab60e360-e826-496d-b382-e868f640d85c
Ks, iωs = let
	θ = (Ri > 1 ? 0 : π/2) + θU
    Kmax = min(1e-2, 2/√3*π/le)
	#min(1e-2, π/le * 2 / norm(inv([1;0;;0.5;√3/2]) * [cos(mod(θ, π/3)); sin(mod(θ, π/3))], 1))
	if dissip_scheme == :biharmonic
		𝕂ᵘ = Vᵘ * le^3
		𝕂ᵇ = Vᵇ * le^3
	else
		𝕂ᵘ = Vᵘ * le
		𝕂ᵇ = Vᵇ * le
	end
    Ks  = range(1e-10, Kmax, 400)
    iωs = Complex{Float64}[]
    for K in Ks
        k = K * cos(θ)
		l = K * sin(θ)

		#jac = [ComplexF64(eady_jac[i,j](k, l, le, Ri, N², g, f₀, 𝕂ᵘ, 𝕂ᵇ, θU, β)) for i=1:size(eady_jac,1), j=1:size(eady_jac,2)]
		jac = [-ComplexF64(fun[i,j](k,l, Ri, le, f₀, g, N², 𝕂ᵘ, 𝕂ᵇ, θU)) for i=1:size(fun,1), j=1:size(fun,2)]

        vals, vecs = eigen(jac)
        push!(iωs, vals[end])
    end
	(Ks, iωs)
end

# ╔═╡ b15d7752-cf88-4e49-95c2-69935c08f448
let
	size   = Ri > 1 ? (1400, 700) : (1400, 500)
	θ = (Ri > 1 ? 0 : π/2) + θU
	fₛ     = min(1e-2, 2/√3*π/le)
	fₛ = 2/√3*π/6.25e3
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, 1.1, -0.1, 0.4) : (0.0, 1.1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	lines!(ax, Ks./ fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(_grid_t)):$(String(hmt_scheme))", linewidth=3)
	axislegend()
	#f[1,2] = Legend(f, ax; merge=true, valign=:top);
	#colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 2ced61a3-8306-423f-b0f7-71ab9906186a
function exactops(grid_t)
	if grid_t == Val(:TriA)
		∇ᵀu⃗ = let
			a = exactop(1, :vertex, TriA.∇ᵀvv(VertexIndex((0,0,1)), vin, u⃗))
			ϵ * (real(a) + Im * imag(a))
		end
		Δu⃗ = let
			as = []
			for iTH=1:2
				a = exactop(1, :vertex, TriA.Δ(VertexIndex((0,0,1)), vin, u⃗[iTH]))
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		Δ²u⃗ = let
			as = []
			for iTH=1:2
				a = exactop(1, :vertex, TriA.Δ(VertexIndex((0,0,1)), v, TriA.Δ(v, vin, u⃗[iTH])))
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
		end
		Δb = let
			a = exactop(1, :vertex, TriA.Δ(VertexIndex((0,0,1)), vin, b))
			ϵ * (real(a) + Im * imag(a))
		end
		Δ²b = let
			a = exactop(1, :vertex, TriA.Δ(VertexIndex((0,0,1)), v, TriA.Δ(v, vin, b)))
			ϵ * (real(a) + Im * imag(a))
		end
		∇p = let
			expr = TriA.∇vv(VertexIndex((0,0,1)), vin, p)
			as = []
			for iTH=1:2
				a = exactop(1, :vertex, expr[iTH])
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		∇η = let
			expr = TriA.∇vv(VertexIndex((0,0,1)), vin, η)
			as = []
			for iTH=1:2
				a = exactop(1, :vertex, expr[iTH])
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		u⃗ᵀ∇u⃗ = let
			as = []
			expr = TriA.u⃗ᵀ∇(VertexIndex((0,0,1)), vin, u⃗, u⃗)
			for iTH=1:2
				a = exactop(1, :vertex, expr[iTH])
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		u⃗ᵀ∇b = let
			a = exactop(1, :vertex, TriA.u⃗∇ᵀ(VertexIndex((0,0,1)), vin, u⃗, b))
			ϵ * (real(a) + Im * imag(a))
		end
		(; ∇ᵀu⃗, Δu⃗, Δ²u⃗, Δb, Δ²b, ∇p, ∇η, u⃗ᵀ∇u⃗, u⃗ᵀ∇b)
	elseif grid_t == Val(:TriB)
		∇ᵀu⃗ = let
			a = exactop(1, :vertex, TriA.∇ᵀvc(VertexIndex((0,0,1)), cin, u⃗))
			ϵ * (real(a) + Im * imag(a))
		end
		Δu⃗ = let
			as = []
			for iTH=1:2
				a = exactop(1, :vertex, TriA.Δ(VertexIndex((0,0,1)), vin, u⃗[iTH]))
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		Δ²u⃗ = let
			as = []
			for iTH=1:2
				bs = []
				for iH=1:dims(:cell)
					b = exactop(1, :vertex, TriA.Δ(CellIndex((0,0,iH)), c, TriA.Δ(c, vin, u⃗[iTH])))
					push!(bs, ϵ * (real(a) + Im * imag(a)))
				end
			end
		end
		Δb = let
			a = exactop(1, :vertex, TriA.Δ(VertexIndex((0,0,1)), vin, b))
			ϵ * (real(a) + Im * imag(a))
		end
		Δ²b = let
			a = exactop(1, :vertex, TriA.Δ(VertexIndex((0,0,1)), v, TriA.Δ(v, vin, b)))
			ϵ * (real(a) + Im * imag(a))
		end
		∇p = let
			expr = TriA.∇vv(VertexIndex((0,0,1)), vin, p)
			as = []
			for iTH=1:2
				a = exactop(1, :vertex, expr[iTH])
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		∇η = let
			expr = TriA.∇vv(VertexIndex((0,0,1)), vin, η)
			as = []
			for iTH=1:2
				a = exactop(1, :vertex, expr[iTH])
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		u⃗ᵀ∇u⃗ = let
			as = []
			expr = TriA.u⃗ᵀ∇(VertexIndex((0,0,1)), vin, u⃗, u⃗)
			for iTH=1:2
				a = exactop(1, :vertex, expr[iTH])
				push!(as, ϵ * (real(a) + Im * imag(a)))
			end
			as
		end
		u⃗ᵀ∇b = let
			a = exactop(1, :vertex, TriA.u⃗∇ᵀ(VertexIndex((0,0,1)), vin, u⃗, b))
			ϵ * (real(a) + Im * imag(a))
		end
		(; ∇ᵀu⃗, Δu⃗, Δ²u⃗, Δb, Δ²b, ∇p, ∇η, u⃗ᵀ∇u⃗, u⃗ᵀ∇b)
	elseif grid_t == Val(:TriC)
		(;)
	else
		(;)
	end
end

# ╔═╡ ee93d7e8-14a0-4021-8f3d-ddf14a79d608
# ╠═╡ disabled = true
#=╠═╡
(; ∇ᵀu⃗, Δu⃗, Δ²u⃗, Δb, Δ²b, ∇p, ∇η, u⃗ᵀ∇u⃗, u⃗ᵀ∇b) = exactops(grid_t)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═500f352c-6e16-11f0-215d-4f5a3075cb33
# ╟─534dfd2f-f7e5-4253-83dc-0be9e94cba01
# ╠═b02f29a5-7f20-49dd-aa0a-9b30d8b4a1ed
# ╠═e653beef-bd70-4f94-867c-b802f416408e
# ╠═99dab097-3490-46cb-bf26-62b60b51d3c2
# ╟─2e485b46-7aa1-483d-9c25-95cd98439ccc
# ╠═5870aa5c-cda6-496f-8836-c9e749e0c629
# ╟─a6e1eda0-c7b6-46be-b435-b999d08d83ba
# ╠═a3eecd80-eb24-47f6-ab77-1bb4b6e109bb
# ╠═a832b824-afdb-458c-9f91-ea023c52852e
# ╠═b9fa1a9b-7989-4f1c-abe8-de153cf8085f
# ╠═0ef5eb5f-6c6c-4ecc-87bd-16247f9056e4
# ╠═0d0a0e92-2d63-4745-9114-8a10c5be58de
# ╠═ee93d7e8-14a0-4021-8f3d-ddf14a79d608
# ╠═4767f71a-fd13-4c26-8ec4-bf16ac6028a4
# ╠═f0f59948-54d8-477e-b2cb-595f68ca13b2
# ╠═e9c5a7ef-6412-4b59-9aa7-d240277550e0
# ╟─6e5eb153-8d92-47f9-bee5-8f823abd11e9
# ╟─826aa7f7-c0f1-4dd4-978f-dc093f970d84
# ╟─b1a33624-23c1-4102-85b7-1e63980f3bf2
# ╠═61150eb4-c9ce-4fd8-8314-780b9c1c29b6
# ╟─29bef65a-0414-4304-a718-189b6885d2e3
# ╟─92853482-fd36-4fe1-b515-9c7696b71c95
# ╟─86ab6145-af2f-4e0d-90f9-4f5492f05ee1
# ╟─6c58a2e9-7b61-415e-b020-20b7170082e2
# ╟─f4fe6284-3ea3-4e57-af2b-664fda7066a3
# ╟─51d642b5-b044-4c2a-b7e0-aab048231544
# ╟─bbf8aa4b-0b66-40b0-b5ab-56f2b4baa641
# ╟─cffd4061-d02d-4e42-bbca-cc990927824b
# ╟─c690b0a4-6130-4320-869f-2764fb681f3f
# ╟─b4024a03-2a30-4756-aebe-3f2e776cd8d6
# ╟─42b320c0-d9a3-4a84-bcc6-3498d71ceec8
# ╟─07e1fa2c-7ef0-4a22-9c75-6640a626fe2c
# ╟─694e01ee-5b10-44e4-a477-e8593ededa8c
# ╟─51c00889-ad8c-400d-b764-a4261e5a76e1
# ╠═4867203f-402e-4d3a-8ad0-ef0ba8478ed1
# ╠═5ad386c4-1dfa-4c20-b7ab-e6a2d13ae1ce
# ╠═a77404c3-511d-4ae6-9e6c-2da0cf490ee9
# ╠═a2687cd8-4e59-414a-8560-e8293935e6b0
# ╟─cf11e62e-ee93-4b2f-9cdc-2cb228ca040c
# ╟─c29f6ae0-e69b-4ad8-8615-8dd9a4aee408
# ╠═28470719-44cf-4b6a-a515-29f7f6a4acb2
# ╠═7c1d69bf-8612-4072-95b5-27b57e59b2ac
# ╠═944eb1fc-3207-48fe-bec3-97f73e7a523a
# ╠═ab60e360-e826-496d-b382-e868f640d85c
# ╠═1aa38ca9-ab99-4d2b-ad5d-c5a18d055229
# ╠═0ad8eba1-f992-4f24-aa90-38ca629f8c72
# ╟─3588e3e2-193d-4c1c-9dbf-c046803cface
# ╠═9bc2c550-4008-42fc-9ac0-607d49f8f319
# ╠═c6332760-298c-42f0-ba13-0936cab3fdcd
# ╟─e3bc1fcc-fb96-404f-8204-675174b3afe1
# ╠═b15d7752-cf88-4e49-95c2-69935c08f448
# ╠═be945a73-9afb-4b9a-ac78-4514fbb0e33e
# ╠═70f46342-4b5c-45ce-9bf9-13a44fd780fc
# ╟─3adadede-5704-441f-9756-08f0f820c723
# ╟─2bfe2cd8-cf70-4efa-a0d2-5a3fe447e9e8
# ╟─69f7116f-5062-43b3-a004-acd5de35ed1e
# ╟─a371a682-d689-4147-8c8f-49a65a29cd8e
# ╟─953ca291-388d-4689-a891-fcdd62b659f7
# ╟─ffc466bf-7d98-4de1-9c12-0e8e2d559ed5
# ╟─9790e852-5292-4ef1-9c11-a49344d423af
# ╟─93f798b4-b999-48e2-9cc4-2371d238ee9d
# ╟─1effe5fd-e458-4bc3-8318-88f54396d32e
# ╟─c68fde99-fd3d-454e-b8d3-c510b35ac02b
# ╟─6a193b8c-8eac-422a-a6b0-cb3962a44c39
# ╟─96d447fc-0e4e-4c2d-8c98-fe50547906e8
# ╟─8d013161-0743-4907-b133-1f2bacbcc1d1
# ╟─72cdb1b1-0ef2-468b-9e25-af3732e00ee8
# ╟─10f7c553-6875-4d22-be26-58947e8a381d
# ╟─b861da6e-16f4-4799-bf5c-fca588b523b6
# ╟─2cdf19d9-870e-495b-bba9-2c1e29f28ba6
# ╟─5f414b9f-342d-4e80-8b72-ebc306b0cc78
# ╟─6e93f103-3988-4a68-8804-eabc14f12f5d
# ╟─5d57177f-c6c5-496f-942f-755b25fa957d
# ╠═abcf94dc-cfac-47e6-b903-7fbf3a63a7aa
# ╟─c9ac2630-ea0e-4048-a161-8c9f925d1d20
# ╠═d35d5d13-d7a2-44ee-8855-be4aa31edf4c
# ╠═3e557bf0-c8aa-49d5-8860-03465dc06731
# ╠═3ff311ec-aaaf-42f1-b11a-25be4632192c
# ╠═4efa8fb3-a6a9-4e73-8d61-8c7c45f0d35f
# ╠═2ced61a3-8306-423f-b0f7-71ab9906186a
