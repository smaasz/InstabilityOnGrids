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
	using CSV
	import DataFrames: subset, DataFrame
	import DataFrames
	import GridOperatorAnalysis: trigwriter, trigexpand, phasesubs
	using HypertextLiteral: @htl, @html_str
	using PlutoUI
	using UUIDs: uuid1
	using RuntimeGeneratedFunctions
	RuntimeGeneratedFunctions.init(@__MODULE__)
end

# ╔═╡ d93466d1-f0da-474c-8e16-2f297806ac80
include(joinpath(@__DIR__, "..", "src", "simulations.jl"));

# ╔═╡ 1e186fab-0a9d-414d-af11-818f9e109087
TableOfContents(; depth=4)

# ╔═╡ 3b35fe39-912a-412d-9406-ecf9d2458bb1
# ╠═╡ disabled = true
#=╠═╡
edf = let
	nθUs = [0, 2]
	nβs  = [0, 5]
	Vᵘs  = [0.001]
	Vᵇs  = [0.001]
	as   = [12.5e3]
	eadyexperiments(; nθUs, nβs, Vᵘs, Vᵇs, as, Nz=64, f₀=-1e-4, g=1e9, N²=1e-6)
end
  ╠═╡ =#

# ╔═╡ b8cce2f9-150a-49fc-beba-5c96126650a4
# ╠═╡ disabled = true
#=╠═╡
begin
	CSV.write("../data/sims3264.csv", df)
end
  ╠═╡ =#

# ╔═╡ 95871d65-f459-4973-83ab-9c1079c0d84b
# ╠═╡ disabled = true
#=╠═╡
df = let
	df = vcat(
		CSV.read("../data/sims32all.csv", DataFrames.DataFrame),
		CSV.read("../data/sims32exact.csv", DataFrames.DataFrame),
		CSV.read("../data/sims64df.csv", DataFrames.DataFrame),
		CSV.read("../data/sims64dc.csv", DataFrames.DataFrame),
		CSV.read("../data/simsA64dc.csv", DataFrames.DataFrame),
		CSV.read("../data/simsA64df.csv", DataFrames.DataFrame),
		CSV.read("../data/simsC64df.csv", DataFrames.DataFrame),
		CSV.read("../data/simsC64dc.csv", DataFrames.DataFrame),
		CSV.read("../data/simsHexC64df.csv", DataFrames.DataFrame),
		CSV.read("../data/simsHexC64dc.csv", DataFrames.DataFrame),
		CSV.read("../data/sims64ddf.csv", DataFrames.DataFrame),
		CSV.read("../data/sims64ddc.csv", DataFrames.DataFrame),
	)
	df =DataFrames.select(df, DataFrames.Not(:Ks, :iωs), :Ks => DataFrames.ByRow(x->eval(Meta.parse(x))) => :Ks, :iωs => DataFrames.ByRow(x-> eval(Meta.parse(x))) => :iωs)
	DataFrames.sort!(df, [:grid_t, :hmt_scheme])
end
  ╠═╡ =#

# ╔═╡ 65657f8a-696d-49d6-81ec-d400d4a7e69f
df = let
	df = CSV.read("../data/sims3264.csv", DataFrames.DataFrame)
	df = DataFrames.select(df, DataFrames.Not(:Ks, :iωs), :Ks => DataFrames.ByRow(x->eval(Meta.parse(x))) => :Ks, :iωs => DataFrames.ByRow(x-> eval(Meta.parse(x))) => :iωs)
	DataFrames.sort!(df, [:grid_t, :hmt_scheme])
end

# ╔═╡ 5983fc73-61b9-4f1f-aaf6-1a528e04ba88
md"""
Grid type: $(@bind sgrid_t PlutoUI.Select([:TriA => "TriA", :TriB => "TriB", :TriC => "TriC", :HexC => "HexC"]))
"""

# ╔═╡ efda4e1f-3a59-485e-b197-db9c4306d045
md"""
## Plotting
"""

# ╔═╡ 5a0e939d-4f53-4cad-b460-998fa572b499
colordf = let
	colors = RGBAf.(Makie.wong_colors())
	grid_ts = vcat([[grid_t for i in 1:length(hmt_schemes[grid_t])] for grid_t in [:HexC, :TriA, :TriB, :TriC]]...)
	hmt_schemes = vcat([first.(hmt_schemes[grid_t]) for grid_t in [:HexC, :TriA, :TriB, :TriC]]...)
	DataFrame("grid_t" => String.(grid_ts), "hmt_scheme" => String.(hmt_schemes), "color" => colors[1:length(grid_ts)])
end

# ╔═╡ 90ee7041-98a0-4f91-9cc2-8831631ee9ec
md"""
### Baroclinic Axis
"""

# ╔═╡ 64430cc5-d472-4645-be3b-54385e754492
md"""
### Symmetric Axis
"""

# ╔═╡ 594b9693-ff74-4166-ba43-51b284ddd34c
md"""
### Galilean Invariance
"""

# ╔═╡ 3204f04f-e3e0-4a55-b975-b3ce4832d5aa
md"""
__Small wavenumber expansion__: $(@bind doapprox PlutoUI.CheckBox())

__Phase substitutions__: $(@bind dophasesubs PlutoUI.CheckBox())
"""

# ╔═╡ bcb90f3f-5c78-4caf-bac7-bff178811320
const phase_subs = phasesubs(a, h, k, l)

# ╔═╡ 357a4379-b325-4146-b474-9a6be5f47a2f
simplifyexpand(x) = simplify(x; expand=true)

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
            	right: 10rem;
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
Number of layers: $(@bind sNz PlutoUI.Slider([32, 64]; show_value=true))
			   
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

# ╔═╡ a8b4e34c-dcb3-49fc-9cfb-736c1c471139
instance = let
	RuntimeGeneratedFunctions.init(@__MODULE__)
	grid_t = sgrid_t
	hmt_scheme = shmt_scheme
	hst_scheme = shst_scheme
	dissip_scheme = :biharmonic
	config = @dict(grid_t, hmt_scheme, hst_scheme, dissip_scheme)
	path = joinpath(@__DIR__, "..", "data")
	data, file = produce_or_load(computesymbols, config, path)
    @unpack fsyms, = data
    fsyms_generated = Dict([name => (@RuntimeGeneratedFunction(fsym[1]), @RuntimeGeneratedFunction(fsym[2])) for (name, fsym) in pairs(fsyms)])
	f₀ = -1e-4
	g = 1e9
	N² = 1e-6
	Ri = sRi
	θU = sθU
	β = sβ
	Vᵘ = sVᵘ
	Vᵇ = sVᵇ
	a = sa
	H = 4000.0
	Nz = sNz
	fullconfig = @dict(grid_t, hmt_scheme, hst_scheme, dissip_scheme, g, f₀, N², Ri, θU, β, Vᵘ, Vᵇ, a, H, Nz)
	analyzeinstability(fullconfig, fsyms_generated; useidealized=Dict(:A⁽ˣ⁾=>true, :A⁽ʸ⁾=>false, :Av⁽ˣ⁾=>false, :Av⁽ʸ⁾=>true, :Gx=>false, :Gy=>false))
end

# ╔═╡ 7b35b917-876d-44dd-9424-3f5c99ff846c
let
	@unpack Ks, iωs, f₀, N² = instance
	
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
			 limits,
			 )
	lines!(ax, Ks./ fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(sgrid_t))", linewidth=3)
	f
end

# ╔═╡ 6ca04b7b-34f8-4072-bccd-487965a4cee2
let
	subdf = subset(df, :grid_t => x->Symbol.(x).==sgrid_t, :Ri => x->x.==sRi, :hst_scheme => x->Symbol.(x).==shst_scheme, :θU => x->x.==sθU, :β => x->x.==sβ, :a => x->x.==sa, :Vᵘ => x->x.==sVᵘ, :Vᵇ => x->x.==sVᵇ, :hmt_scheme => x->Symbol.(x).==shmt_scheme)
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
			 limits,
			 )
	lines!(ax, Ks./ fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(sgrid_t))", linewidth=3)
	f
end

# ╔═╡ 626ea509-d106-40e9-904f-0ff96912bd94
let
	subdf = subset(df, :Nz=>x->x.==sNz, :β=>x->x.==sβ, :θU=>x->x.==sθU, :a=>x->x.==sa, :Vᵘ=>x->x.≈0.0, :Vᵇ=>x->x.≈0.0, :Ri=>x->x.==100, :hst_scheme=>x->Symbol.(x).==shst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	#DataFrames.innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
end

# ╔═╡ ef767b5f-bd4b-45f9-9600-1cd6f7644765
let
	Ri     = 100
	size   = Ri > 1 ? (1400, 500) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/sa)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 3.4 : 2.5
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :Nz=>x->x.==sNz, :β=>x->x.==sβ, :θU=>x->x.==sθU, :a=>x->x.==sa, :Vᵘ=>x->x.≈0.0, :Vᵇ=>x->x.≈0.0, :Ri=>x->x.==100, :hst_scheme=>x->Symbol.(x).==shst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	subdf = DataFrames.innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
		M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t))"* (grid_t == "TriA" ? "" : ":$(String(hmt_scheme))"), linewidth=3, color=color)
	end
	(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈100))
	M²= √(N² * f₀^2 / Ri)
	lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
	#axislegend()
	axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ ec9534ea-68b5-473c-913a-133b02bed327
let
	Ri     = 100
	size   = Ri > 1 ? (1400, 500) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/sa)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 3.4 : 2.5
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :Nz=>x->x.==sNz, :β=>x->x.==sβ, :θU=>x->x.==sθU, :a=>x->x.==sa, :Vᵘ=>x->x.≈sVᵘ, :Vᵇ=>x->x.≈sVᵇ, :Ri=>x->x.==100, :hst_scheme=>x->Symbol.(x).==shst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	subdf = DataFrames.innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
		M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t))"* (grid_t == "TriA" ? "" : ":$(String(hmt_scheme))"), linewidth=3, color=color)
	end
	(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈100))
	M²= √(N² * f₀^2 / Ri)
	lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
	axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ a9222090-eb04-4086-8fc7-45afe40ba590
let
	Ri     = 0.5
	size   = Ri > 1 ? (1400, 500) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/sa)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 3.4 : 3.4
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :Nz=>x->x.==sNz, :β=>x->x.==sβ, :θU=>x->x.==sθU, :a=>x->x.==sa, :Vᵘ=>x->x.≈0.0, :Vᵇ=>x->x.≈0.0, :Ri=>x->x.==Ri, :hst_scheme=>x->Symbol.(x).==shst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	subdf = DataFrames.innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
		M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t))"* (grid_t == "TriA" ? "" : ":$(String(hmt_scheme))"), linewidth=3, color=color)
	end
	(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈Ri))
	M²= √(N² * f₀^2 / Ri)
	lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
	axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ feeea627-dcca-45dd-b158-a167fbbda814
let
	Ri     = 0.5
	size   = Ri > 1 ? (1400, 500) : (1400, 500)
	fₛ     = min(1e-2, 2/√3*π/sa)
	xticks = if Ri > 1
		xs = collect(0.0:1/8:1.1)
    	ls = ["0.0", "1/8", "1/4","3/8", "1/2","5/8", "3/4","7/8", "1"]
    	(xs, ls)
	else
		xs = collect(0.0:1/4:1.1)
    	ls = ["0.0", "1/4", "1/2", "3/4", "1"]
		(xs, ls)
	end
	aspect = Ri > 1 ? 3.4 : 3.4
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)
	
	f = Figure(; size, fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / fₛ",
			 ylabel = "growth rate / N M⁻²",
			 xticks,
			 aspect,
			 limits,
			 )
	subdf = subset(df, :Nz=>x->x.==sNz, :β=>x->x.==sβ, :θU=>x->x.==sθU, :a=>x->x.==sa, :Vᵘ=>x->x.≈sVᵘ, :Vᵇ=>x->x.≈sVᵇ, :Ri=>x->x.==Ri, :hst_scheme=>x->Symbol.(x).==shst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
	subdf = DataFrames.innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
	for row in eachrow(subdf)
		(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
		M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="$(String(grid_t))"* (grid_t == "TriA" ? "" : ":$(String(hmt_scheme))"), linewidth=3, color=color)
	end
	(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈Ri))
	M²= √(N² * f₀^2 / Ri)
	lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
	axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30);
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 81110f20-424e-4a29-a39e-eb8083615d6b
let
	Ri = 100
	size   = Ri > 1 ? (1400, 2*500) : (1400, 2*500)
	fₛ     = min(1e-2, 2/√3*π/sa)
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
	aspect = Ri > 1 ? 3.4 : 2.5
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)
	
	f = Figure(; size, fontsize=36)
	for (i, θU) in enumerate([0, π/6])
		ax = Axis(f[i,1];
				 xlabel = "wavenumber / fₛ",
				 ylabel = "growth rate / N M⁻²",
				 xticks,
				 aspect,
				 limits,
				 )
		subdf = subset(df, :Nz=>x->x.==sNz, :β=>x->x.≈sβ, :θU=>x->x.≈θU, :a=>x->x.==sa, :Vᵘ=>x->x.≈sVᵘ, :Vᵇ=>x->x.≈sVᵇ, :Ri=>x->x.==Ri, :hst_scheme=>x->Symbol.(x).==shst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
		subdf = DataFrames.innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
		for row in eachrow(subdf)
			(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
			M² = √(N² * f₀^2 / Ri)
			lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label= grid_t == "TriA" ? "$(String(grid_t))" : "$(String(grid_t)):$(String(hmt_scheme))", linewidth=3, color=color)
		end
		(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈Ri))
	M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
		axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30)
		#Label(f[i,1, TopLeft()], "($i)", )
	end
	#axislegend()
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ dd6e6622-53c3-4d7e-bcb0-ddbe7fce25fd
let
	Ri = 0.5
	size   = Ri > 1 ? (1400, 2*500) : (1400, 2*500)
	fₛ     = min(1e-2, 2/√3*π/sa)
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
	aspect = Ri > 1 ? 3.4 : 3.0
	limits = Ri > 1 ? (0.0, 1.01, -0.02, 0.38) : (0.0, 1.01, -0.1, 1.0)
	
	f = Figure(; size, fontsize=36)
	for (i, θU) in enumerate([0, π/6])
		ax = Axis(f[i,1];
				 xlabel = "wavenumber / fₛ",
				 ylabel = "growth rate / N M⁻²",
				 xticks,
				 aspect,
				 limits,
				 )
		subdf = subset(df, :Nz=>x->x.==sNz, :β=>x->x.≈sβ, :θU=>x->x.≈θU, :a=>x->x.==sa, :Vᵘ=>x->x.≈sVᵘ, :Vᵇ=>x->x.≈sVᵇ, :Ri=>x->x.==Ri, :hst_scheme=>x->Symbol.(x).==shst_scheme, :dissip_scheme=>x->Symbol.(x).==:biharmonic)
		subdf = DataFrames.innerjoin(subdf, colordf; on=[:grid_t, :hmt_scheme])
		for row in eachrow(subdf)
			(; Ks, iωs, grid_t, hmt_scheme, f₀, N², color) = row
			M² = √(N² * f₀^2 / Ri)
			lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label= grid_t == "TriA" ? "$(String(grid_t))" : "$(String(grid_t)):$(String(hmt_scheme))", linewidth=3, color=color)
		end
		(; Ks, iωs, grid_t, hmt_scheme, N², f₀) = first(subset(df, :a=>x->x.≈1e-20, :Ri=>x->x.≈0.5))
	M²= √(N² * f₀^2 / Ri)
		lines!(ax, Ks./fₛ, real.(iωs) .* (sqrt(N²) / abs(M²)), label="ideal", linewidth=4, linestyle=:dot, color=:black)
		axislegend(ax; merge=true, valign=:top, orientation=:horizontal, labelsize=30)
		#Label(f[i,1, TopLeft()], "($i)", )
	end
	#axislegend()
	colsize!(f.layout, 1, Aspect(1, aspect))
	f
end

# ╔═╡ 38b5ec65-ab21-4a03-8e02-a4811855da5a
fsymbols = let
	fsymbols = symbolicsymbols(sgrid_t, shmt_scheme, shst_scheme, :biharmonic)
	for (name, fsymbol) in fsymbols
        fsymbol = substitute.(fsymbol, Ref(sqrt3_subs))
        fsymbol = substitute.(fsymbol, Ref(Dict(sqrt3 => h / a * 2))) .|> Symbolics.wrap
        if doapprox
            fsymbol = simplify.(fsymbol; rewriter=trigexpand(4)) .|> simplifyexpand
            fsymbol = substitute.(fsymbol, Ref(Dict(a^2 => 4//3 * h^2, a^4 => 16//9 * h^4))) .|> simplifyexpand
        elseif dophasesubs
            fsymbol = substitute.(fsymbol, Ref(phase_subs)) .|> simplifyexpand
            fsymbol = simplify.(fsymbol; rewriter=trigwriter)
        else
            fsymbol = simplify.(fsymbol; rewriter=trigwriter)
        end
        fsymbols[name] = convert.(Complex{Num}, fsymbol)
    end
	fsymbols
end;

# ╔═╡ f4415ce5-e55f-4f7c-bacb-a304fec58d3e
real.(fsymbols[:A⁽ˣ⁾][:,:,1]) + im*imag.(fsymbols[:A⁽ˣ⁾][:,:,1]) .|> Complex{Num}

# ╔═╡ cd276d71-b1eb-4699-a8f8-0aa0f60e0226
real.(fsymbols[:Γx][:,:,1]) + im*imag.(fsymbols[:Γx][:,:,1]) .|> Complex{Num}

# ╔═╡ Cell order:
# ╠═92b5008e-bfcc-11f0-0d84-7baec76549aa
# ╟─1e186fab-0a9d-414d-af11-818f9e109087
# ╠═d93466d1-f0da-474c-8e16-2f297806ac80
# ╠═a8b4e34c-dcb3-49fc-9cfb-736c1c471139
# ╟─7b35b917-876d-44dd-9424-3f5c99ff846c
# ╠═3b35fe39-912a-412d-9406-ecf9d2458bb1
# ╠═b8cce2f9-150a-49fc-beba-5c96126650a4
# ╠═95871d65-f459-4973-83ab-9c1079c0d84b
# ╠═65657f8a-696d-49d6-81ec-d400d4a7e69f
# ╟─5983fc73-61b9-4f1f-aaf6-1a528e04ba88
# ╟─9e0ed423-7d1e-4928-a9ab-e89a7d0f9975
# ╟─6ca04b7b-34f8-4072-bccd-487965a4cee2
# ╟─efda4e1f-3a59-485e-b197-db9c4306d045
# ╠═5a0e939d-4f53-4cad-b460-998fa572b499
# ╠═626ea509-d106-40e9-904f-0ff96912bd94
# ╟─90ee7041-98a0-4f91-9cc2-8831631ee9ec
# ╟─ef767b5f-bd4b-45f9-9600-1cd6f7644765
# ╟─ec9534ea-68b5-473c-913a-133b02bed327
# ╟─64430cc5-d472-4645-be3b-54385e754492
# ╟─a9222090-eb04-4086-8fc7-45afe40ba590
# ╟─feeea627-dcca-45dd-b158-a167fbbda814
# ╟─594b9693-ff74-4166-ba43-51b284ddd34c
# ╟─81110f20-424e-4a29-a39e-eb8083615d6b
# ╟─dd6e6622-53c3-4d7e-bcb0-ddbe7fce25fd
# ╟─3204f04f-e3e0-4a55-b975-b3ce4832d5aa
# ╠═bcb90f3f-5c78-4caf-bac7-bff178811320
# ╠═357a4379-b325-4146-b474-9a6be5f47a2f
# ╠═38b5ec65-ab21-4a03-8e02-a4811855da5a
# ╠═f4415ce5-e55f-4f7c-bacb-a304fec58d3e
# ╠═cd276d71-b1eb-4699-a8f8-0aa0f60e0226
# ╟─9e04b993-5103-4f55-b9d7-56f14cc28afd
# ╟─7ce186b6-f5e7-48bb-8f15-e2e5708f7e93
# ╟─d83ec40d-0151-490a-a4b1-84459432896f
