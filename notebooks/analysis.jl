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
	import GridOperatorAnalysis: bb, eady_jacobian#, eady_jacobian_wo_mtk, 
	import GridOperatorAnalysis
	using HypertextLiteral: @htl, @html_str
	using LaTeXStrings
	import LinearAlgebra: eigen, norm
	using PlutoUI
	using RuntimeGeneratedFunctions
	import Symbolics: substitute
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

# ╔═╡ a6e1eda0-c7b6-46be-b435-b999d08d83ba
md"""
__Choice of the grid__: $(@bind grid_t PlutoUI.Select([:TriA => "triangular A-Grid", :TriC => "triangular C-Grid", :TriB => "triangular B-Grid", :HexC => "hexagonal C-Grid"]; default=:TriA))
"""

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

# ╔═╡ 0389932f-20e8-4346-95a2-051249ac3025
WideCell(
md"""
__Horizontal momentum transport balance__:
	
   $(@bind test1 Select([:a => "∂ₜ(evalat(vout, vin, u⃗[iTH]))"])) 
\+ $(@bind test2 Select([:a => "TriA.u⃗ᵀ∇(vout, vin, u⃗, u⃗)[iTH]", :b => "evalat(vout, vin, u⃗ᵀ∇u⃗[iTH])"]))
\+ $(@bind test2 Select([:a => "f₀ * evalat(vout, vin, u⃗⊥[iTH])"]))
\+ $(@bind test4 Select([:a => "evalat(vout, vin, w) * ∂₃(evalat(vout, vin, u⃗[iTH]))"]))
\+ $(@bind test5 Select([:a => "TriA.∇vv(vout, vin, p)[iTH]",     :b => "evalat(vout, vin, ∇p[iTH])"]))
\+ $(@bind test6 Select([:a => "g * TriA.∇vv(vout, vin, η)[iTH]", :b => "g * evalat(vout, vin, ∇η[iTH])"]))
\+ $(@bind test7 Select([:a => "-Dᵘₕ[iTH]", :b => "evalat(vout, vin, Δ⃗u⃗[iTH])"]))
= 0
"""
; max_width=1500)

# ╔═╡ bbf8aa4b-0b66-40b0-b5ab-56f2b4baa641
WideCell(
md"""
__Buoyancy transport balance__:
	
   $(@bind b1 Select([:a => "∂ₜ(evalat(vout, vin, b))"])) 
\+ $(@bind b2 Select([:a => "TriA.u⃗∇ᵀ(vout, vin, u⃗, b)", :b => "evalat(vout, vin, u⃗ᵀ∇u⃗[iTH])"]))
\+ $(@bind b3 Select([:a => "evalat(vout, vin, w) * ∂₃(evalat(vout, vin, b))"]))
\+ $(@bind b4 Select([:a => "-Dᵇₕ", :b => "evalat(vout, vin, Δb)"]))
= 0
"""
; max_width=1500)

# ╔═╡ cffd4061-d02d-4e42-bbca-cc990927824b
WideCell(
md"""
__Surface elevation equation__:
	
   $(@bind η1 Select([:a => "∂ₜ(evalat(vout, vin, η))"])) 
\+ $(@bind η2 Select([:a => "evalat(vout, vin, ∫∇ᵀu⃗dz)"]))
= 0
"""
; max_width=1500)

# ╔═╡ c690b0a4-6130-4320-869f-2764fb681f3f
WideCell(
md"""
__Continuity Equation__:
	
   $(@bind c1 Select([:a => "TriA.∇ᵀvv(vout, vin, u⃗)", :b => "evalat(vout, vin, ∇ᵀu⃗)"])) 
\+ $(@bind c2 Select([:a => "∂₃(evalat(vout, vin, w))"]))
= 0
"""
; max_width=1500)

# ╔═╡ b4024a03-2a30-4756-aebe-3f2e776cd8d6
WideCell(
md"""
__Hydrostatic Balance__:
	
   $(@bind p1 Select([:a => "∂₃(evalat(vout, vin, p))"])) 
\+ $(@bind p2 Select([:a => "-evalat(vout, vin, b)"]))
= 0
"""
; max_width=1500)

# ╔═╡ 42b320c0-d9a3-4a84-bcc6-3498d71ceec8
md"""
__Discretization scheme for the momentum transport__: 
$(@bind hmt_scheme select_hmt_scheme(grid_t))
"""

# ╔═╡ 07e1fa2c-7ef0-4a22-9c75-6640a626fe2c
md"""
__Discretization scheme for the momentum transport__: 
$(@bind hst_scheme select_hst_scheme(grid_t))
"""

# ╔═╡ 694e01ee-5b10-44e4-a477-e8593ededa8c
md"""
__Dissipation scheme__: $(@bind dissip_scheme PlutoUI.Select([:harmonic => "harmonic", :biharmonic => "biharmonic"]; default=:biharmonic))
"""

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

# ╔═╡ 8ef7642b-a0f8-4948-a505-4da71e1dc2ef
begin
	RuntimeGeneratedFunctions.init(@__MODULE__)
	eady_jac =  [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
end;

# ╔═╡ 1aa38ca9-ab99-4d2b-ad5d-c5a18d055229
Symbolics.@variables _Ri _le _g _f₀ _N² _𝕂ᵘ _𝕂ᵇ _θU _β k l

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

le: $(@bind le PlutoUI.Slider([1e-20, 1.575e3, 3.125e3, 6.25e3, 12.5e3, 25e3]; default=6.25e3, show_value=true))

Ri: $(@bind Ri Select([100//1, 1//2]; default=100//1))

Vᵘ: $(@bind Vᵘ PlutoUI.Slider([0, 1e-3, 5e-3, 1e-2]; show_value=true))

Vᵇ: $(@bind Vᵇ PlutoUI.Slider([0, 1e-3, 5e-3, 1e-2]; show_value=true))
""")

# ╔═╡ 9bc2c550-4008-42fc-9ac0-607d49f8f319
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

# ╔═╡ c6332760-298c-42f0-ba13-0936cab3fdcd
fun = let
	Symbolics.build_function.(substitute.(jac, Ref(Dict(GridOperatorAnalysis.sqrt3=>√3))), k, l, _Ri, _le, _f₀, _g, _N², _𝕂ᵘ, _𝕂ᵇ; expression=Val{false})
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
		jac = [-Symbolics.unwrap(fun[i,j](k,l, Ri, le, f₀, g, N², 𝕂ᵘ, 𝕂ᵇ)) for i=1:size(fun,1), j=1:size(fun,2)]

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
	lines!(ax, Ks./ fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t)):$(String(hmt_scheme))", linewidth=3)
	axislegend()
	#f[1,2] = Legend(f, ax; merge=true, valign=:top);
	#colsize!(f.layout, 1, Aspect(1, aspect))
	f
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

# ╔═╡ 2cdf19d9-870e-495b-bba9-2c1e29f28ba6
let
	size   = Ri > 1 ? (1400, 700) : (1400, 500)
	fₛ     = Ks[end]
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "fₛ/8", "fₛ/4","3fₛ/8", "fₛ/2","5fₛ/8", "3fₛ/4","7fₛ/8", "fₛ"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "fₛ/4", "fₛ/2", "3fₛ/4", "fₛ"]
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

# ╔═╡ Cell order:
# ╠═500f352c-6e16-11f0-215d-4f5a3075cb33
# ╟─534dfd2f-f7e5-4253-83dc-0be9e94cba01
# ╠═b02f29a5-7f20-49dd-aa0a-9b30d8b4a1ed
# ╠═e653beef-bd70-4f94-867c-b802f416408e
# ╠═99dab097-3490-46cb-bf26-62b60b51d3c2
# ╟─2e485b46-7aa1-483d-9c25-95cd98439ccc
# ╟─a6e1eda0-c7b6-46be-b435-b999d08d83ba
# ╟─6e5eb153-8d92-47f9-bee5-8f823abd11e9
# ╟─826aa7f7-c0f1-4dd4-978f-dc093f970d84
# ╟─0389932f-20e8-4346-95a2-051249ac3025
# ╟─bbf8aa4b-0b66-40b0-b5ab-56f2b4baa641
# ╟─cffd4061-d02d-4e42-bbca-cc990927824b
# ╟─c690b0a4-6130-4320-869f-2764fb681f3f
# ╟─b4024a03-2a30-4756-aebe-3f2e776cd8d6
# ╟─42b320c0-d9a3-4a84-bcc6-3498d71ceec8
# ╟─07e1fa2c-7ef0-4a22-9c75-6640a626fe2c
# ╟─694e01ee-5b10-44e4-a477-e8593ededa8c
# ╟─cf11e62e-ee93-4b2f-9cdc-2cb228ca040c
# ╟─c29f6ae0-e69b-4ad8-8615-8dd9a4aee408
# ╠═28470719-44cf-4b6a-a515-29f7f6a4acb2
# ╠═7c1d69bf-8612-4072-95b5-27b57e59b2ac
# ╠═944eb1fc-3207-48fe-bec3-97f73e7a523a
# ╠═8ef7642b-a0f8-4948-a505-4da71e1dc2ef
# ╠═ab60e360-e826-496d-b382-e868f640d85c
# ╠═1aa38ca9-ab99-4d2b-ad5d-c5a18d055229
# ╠═0ad8eba1-f992-4f24-aa90-38ca629f8c72
# ╟─3588e3e2-193d-4c1c-9dbf-c046803cface
# ╠═9bc2c550-4008-42fc-9ac0-607d49f8f319
# ╠═c6332760-298c-42f0-ba13-0936cab3fdcd
# ╟─e3bc1fcc-fb96-404f-8204-675174b3afe1
# ╟─b15d7752-cf88-4e49-95c2-69935c08f448
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
# ╟─10f7c553-6875-4d22-be26-58947e8a381d
# ╠═96d447fc-0e4e-4c2d-8c98-fe50547906e8
# ╟─8d013161-0743-4907-b133-1f2bacbcc1d1
# ╠═2cdf19d9-870e-495b-bba9-2c1e29f28ba6
# ╟─5f414b9f-342d-4e80-8b72-ebc306b0cc78
# ╟─6e93f103-3988-4a68-8804-eabc14f12f5d
# ╟─5d57177f-c6c5-496f-942f-755b25fa957d
# ╠═abcf94dc-cfac-47e6-b903-7fbf3a63a7aa
# ╟─c9ac2630-ea0e-4048-a161-8c9f925d1d20
