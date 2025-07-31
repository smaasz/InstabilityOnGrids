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
	import GridOperatorAnalysis: eady_jacobian, bb
	import GridOperatorAnalysis
	using HypertextLiteral: @htl, @html_str
	import LinearAlgebra: eigen, norm
	using PlutoUI
	using RuntimeGeneratedFunctions
	import Symbolics: substitute
	import Symbolics
	using UUIDs: uuid1
end;

# ╔═╡ 534dfd2f-f7e5-4253-83dc-0be9e94cba01
md"""
## Constants
"""

# ╔═╡ b02f29a5-7f20-49dd-aa0a-9b30d8b4a1ed
begin
	const f₀ = -1e-4
	const g  = 1e9
	const N² = 1e-6
	const U  = 1.0
	const H  = 4000
	const Nz = 8
end;

# ╔═╡ 2e485b46-7aa1-483d-9c25-95cd98439ccc
md"""
## Discretization
"""

# ╔═╡ a6e1eda0-c7b6-46be-b435-b999d08d83ba
md"""
__Choice of the grid__: $(@bind grid_t PlutoUI.Select([:TriA => "triangular A-Grid", :TriC => "triangular C-Grid", :TriB => "triangular B-Grid", :HexC => "hexagonal C-Grid"]; default=:TriA))
"""

# ╔═╡ 6e5eb153-8d92-47f9-bee5-8f823abd11e9
function select_scheme(grid_t)
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

# ╔═╡ 42b320c0-d9a3-4a84-bcc6-3498d71ceec8
md"""
__Discretization scheme for the momentum transport__: 
$(@bind hmt_scheme select_scheme(grid_t))
"""

# ╔═╡ cf11e62e-ee93-4b2f-9cdc-2cb228ca040c
md"""
## Instability Analysis
"""

# ╔═╡ 28470719-44cf-4b6a-a515-29f7f6a4acb2
config = @dict(grid_t, hmt_scheme, Nz, H);

# ╔═╡ 7c1d69bf-8612-4072-95b5-27b57e59b2ac
path = joinpath(@__DIR__, "..", "data")

# ╔═╡ 944eb1fc-3207-48fe-bec3-97f73e7a523a
begin
	data, file = produce_or_load(config, path) do config
		@unpack grid_t, hmt_scheme, Nz, H = config
		eady_jac_ex, eady_sys = eady_jacobian(Val(grid_t); Nz, H, U, hmt_scheme)
		@strdict(eady_jac_ex, eady_sys)
	end
	@unpack eady_jac_ex, eady_sys = data
end;

# ╔═╡ 8ef7642b-a0f8-4948-a505-4da71e1dc2ef
begin
	RuntimeGeneratedFunctions.init(@__MODULE__)
	eady_jac =  [@RuntimeGeneratedFunction(GridOperatorAnalysis, eady_jac_ex[i,j]) for i=1:size(eady_jac_ex,1), j=1:size(eady_jac_ex,2)]
end;

# ╔═╡ e3bc1fcc-fb96-404f-8204-675174b3afe1
md"""
## Plotting
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

le: $(@bind le PlutoUI.Slider([1e-9, 1.575e3, 3.125e3, 6.25e3, 12.5e3]; default=6.25e3, show_value=true))

Ri: $(@bind Ri Select([100//1, 1//2]; default=100//1))

𝕂ᵘ: $(@bind 𝕂ᵘ PlutoUI.Slider([0, 1e0, 1e1, 1e2, 1e3]; show_value=true))

𝕂ᵇ: $(@bind 𝕂ᵇ PlutoUI.Slider([0, 1e0, 1e1, 1e2, 1e3]; show_value=true))
""")

# ╔═╡ ab60e360-e826-496d-b382-e868f640d85c
Ks, iωs = let
	θ = (Ri > 1 ? 0 : π/2) + θU
    Kmax = let  # check if correct
	_bb = Symbolics.unwrap.(substitute.(bb, Ref(Dict(GridOperatorAnalysis.sqrt3 => √3, GridOperatorAnalysis.le => le))))
	π / norm(_bb * [cos(θ); sin(θ)], Inf)
    end

    Ks  = range(1e-10, Kmax, 400)
    iωs = Complex{Float64}[]
    for K in Ks
        k = K * cos(θ)
		l = K * sin(θ)

		jac = [Complex{Float64}(eady_jac[i,j](k, l, le, Ri, N², g, f₀, 𝕂ᵘ, 𝕂ᵇ, θU, β)) for i=1:size(eady_jac,1), j=1:size(eady_jac,2)]

        vals, vecs = eigen(jac)
        push!(iωs, vals[end])
    end
	(Ks, iωs)
end

# ╔═╡ b15d7752-cf88-4e49-95c2-69935c08f448
let
	size   = Ri > 1 ? (1400, 700) : (1300, 800)
	fₛ     = Ks[end]
	xticks = if Ri > 1
		xs = collect(0.0:fₛ/8:(fₛ*1.1))
    	ls = ["0.0", "fₛ/8", "fₛ/4","3fₛ/8", "fₛ/2","5fₛ/8", "3fₛ/4","7fₛ/8", "fₛ"]
    	(xs, ls)
	else
		xs = collect(0.0:fₛ/4:(fₛ*1.1))
    	ls = ["0.0", "fₛ/4", "fₛ/2", "3fₛ/4", "fₛ"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 1.5 : 2.5
	limits = Ri > 1 ? (0.0, xs[end], -0.1, 0.4) : (0.0, xs[end], -0.1, 1.0)
	M²     = √(N² * f₀^2 / Ri)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	lines!(ax, Ks, real.(iωs) .* (sqrt(N²) / abs(M²)), label="growth rate", linewidth=3)
	#f[1,2] = Legend(f, ax; merge=true, valign=:top);
	#colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ Cell order:
# ╠═500f352c-6e16-11f0-215d-4f5a3075cb33
# ╟─534dfd2f-f7e5-4253-83dc-0be9e94cba01
# ╠═b02f29a5-7f20-49dd-aa0a-9b30d8b4a1ed
# ╟─2e485b46-7aa1-483d-9c25-95cd98439ccc
# ╟─a6e1eda0-c7b6-46be-b435-b999d08d83ba
# ╟─6e5eb153-8d92-47f9-bee5-8f823abd11e9
# ╟─42b320c0-d9a3-4a84-bcc6-3498d71ceec8
# ╟─cf11e62e-ee93-4b2f-9cdc-2cb228ca040c
# ╟─c29f6ae0-e69b-4ad8-8615-8dd9a4aee408
# ╠═28470719-44cf-4b6a-a515-29f7f6a4acb2
# ╠═7c1d69bf-8612-4072-95b5-27b57e59b2ac
# ╠═944eb1fc-3207-48fe-bec3-97f73e7a523a
# ╠═8ef7642b-a0f8-4948-a505-4da71e1dc2ef
# ╟─ab60e360-e826-496d-b382-e868f640d85c
# ╟─e3bc1fcc-fb96-404f-8204-675174b3afe1
# ╟─b15d7752-cf88-4e49-95c2-69935c08f448
# ╟─5f414b9f-342d-4e80-8b72-ebc306b0cc78
# ╟─6e93f103-3988-4a68-8804-eabc14f12f5d
# ╟─5d57177f-c6c5-496f-942f-755b25fa957d
