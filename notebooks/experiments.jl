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

# ╔═╡ 92b5008e-bfcc-11f0-0d84-7baec76549aa
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	using CairoMakie
	import DataFrames: subset
	using HypertextLiteral: @htl, @html_str
	using PlutoUI
	using UUIDs: uuid1
end

# ╔═╡ d93466d1-f0da-474c-8e16-2f297806ac80
include(joinpath(@__DIR__, "..", "src", "simulations.jl"));

# ╔═╡ 1e186fab-0a9d-414d-af11-818f9e109087
TableOfContents(; depth=4)

# ╔═╡ 3b35fe39-912a-412d-9406-ecf9d2458bb1
ndf = let
	nθUs = [0]
	nβs  = [0]
	Vᵘs  = [0.0]
	Vᵇs  = [0.0]
	as   = [6.25e3, 12.5e3]
	eadyexperiments(; nθUs, nβs, Vᵘs, Vᵇs, as, Nz=32, f₀=-1e-4, g=1e9, N²=1e-6)
end

# ╔═╡ b8cce2f9-150a-49fc-beba-5c96126650a4
# ╠═╡ disabled = true
#=╠═╡
begin
	using CSV
	CSV.write("../data/sims32t.csv", df)
end
  ╠═╡ =#

# ╔═╡ 5983fc73-61b9-4f1f-aaf6-1a528e04ba88
md"""
Grid type: $(@bind sgrid_t PlutoUI.Select([:TriA => "TriA", :TriB => "TriB", :TriC => "TriC", :HexC => "HexC"]))
"""

# ╔═╡ 38b5ec65-ab21-4a03-8e02-a4811855da5a
fsyms = symbolicsymbols(sgrid_t, :fdcre, :low, :biharmonic);

# ╔═╡ 9e04b993-5103-4f55-b9d7-56f14cc28afd
html"""<hr>"""

# ╔═╡ 7ce186b6-f5e7-48bb-8f15-e2e5708f7e93
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ d83ec40d-0151-490a-a4b1-84459432896f
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
            	right: 80rem;
            	top: $(top)px;
            	width: 230px;
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

# ╔═╡ 9e0ed423-7d1e-4928-a9ab-e89a7d0f9975
begin
	floataside(md"""
Number of layers: $(@bind sNz PlutoUI.Slider([16, 32, 64]; show_value=true))
			   
Flow direction: $(@bind sθU PlutoUI.Slider([0.0, π/12, π/6]; show_value=true))

Flow shift (in -HM²/f₀): $(@bind sβ PlutoUI.Slider(0.0:0.5:0.5; default=0.0, show_value=true))

a: $(@bind sa PlutoUI.Slider([6.25e3, 12.5e3]; default=6.25e3, show_value=true))

Ri: $(@bind sRi Select([100.0, 0.5]; default=100//1))

hmt-scheme: $(@bind shmt_scheme Select(hmt_schemes[sgrid_t]))
			   
hst-scheme: $(@bind shst_scheme Select([:low => "low-order accurate", :high => "high-order accurate"]))

Vᵘ: $(@bind sVᵘ PlutoUI.Slider([0, 1e-3, 1e-2]; show_value=true))

Vᵇ: $(@bind sVᵇ PlutoUI.Slider([0, 1e-3, 1e-2]; show_value=true))
""")
end

# ╔═╡ 6ca04b7b-34f8-4072-bccd-487965a4cee2
let
	subdf = subset(ndf, :grid_t => x->x.==sgrid_t, :Ri => x->x.==sRi, :hst_scheme => x->x.==shst_scheme, :θU => x->x.==sθU, :β => x->x.==sβ, :a => x->x.==sa, :Vᵘ => x->x.==sVᵘ, :Vᵇ => x->x.==sVᵇ, :hmt_scheme => x->x.==shmt_scheme)
	(; Ks, iωs, f₀, N²) = first(subdf)
	
	size   = sRi > 1 ? (1500, 550) : (1400, 500)
	θ      = (sRi > 1 ? 0 : π/2) + sθU
	fₛ     = min(1e-2, 2/√3*π/sa)
	#fₛ = 2/√3*π/6.25e3
	xticks = if sRi > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = sRi > 1 ? 3.0 : 2.5
	limits = sRi > 1 ? (0.0, 1.0, -0.1, 0.4) : (0.0, 1.1, -0.1, 1.0)
	M²     = √(N² * f₀^2 / sRi)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 #limits,
			 )
	lines!(ax, Ks./ fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(sgrid_t))", linewidth=3)
	f
end

# ╔═╡ Cell order:
# ╠═92b5008e-bfcc-11f0-0d84-7baec76549aa
# ╟─1e186fab-0a9d-414d-af11-818f9e109087
# ╠═d93466d1-f0da-474c-8e16-2f297806ac80
# ╠═3b35fe39-912a-412d-9406-ecf9d2458bb1
# ╠═b8cce2f9-150a-49fc-beba-5c96126650a4
# ╟─5983fc73-61b9-4f1f-aaf6-1a528e04ba88
# ╟─9e0ed423-7d1e-4928-a9ab-e89a7d0f9975
# ╠═6ca04b7b-34f8-4072-bccd-487965a4cee2
# ╠═38b5ec65-ab21-4a03-8e02-a4811855da5a
# ╟─9e04b993-5103-4f55-b9d7-56f14cc28afd
# ╟─7ce186b6-f5e7-48bb-8f15-e2e5708f7e93
# ╟─d83ec40d-0151-490a-a4b1-84459432896f
