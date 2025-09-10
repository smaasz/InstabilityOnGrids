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

# ╔═╡ f59a5438-8d5e-11f0-13fa-d9c703e5f87f
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	using Revise
	import GridOperatorAnalysis: eady_background_flow, bb, e, c, v, nS, dims, t, z 
	import GridOperatorAnalysis: TriAFlow, TriAFlowFT, TriBFlow, TriBFlowFT, TriCFlow, TriCFlowFT
	import GridOperatorAnalysis: colpt_type, colptidx, compute_phases
	import GridOperatorAnalysis: fourier_transform_expression
	import GridOperatorAnalysis.TriA
	import GridOperatorAnalysis.TriB
	import GridOperatorAnalysis.TriC
	import GridOperatorAnalysis
	import Symbolics: substitute, @variables, taylor_coeff, simplify, coeff, taylor, expand
	import Symbolics
	import SymbolicUtils
	using HypertextLiteral: @htl, @html_str
	using PlutoUI
	using UUIDs: uuid1
end;

# ╔═╡ b9cca8c1-0d0f-48f9-b21b-8d6df2cb77aa
md"""
# Calculation of Fourier Symbols
"""

# ╔═╡ 59056484-9c6b-48ad-b741-2bc294d5cc6f
TableOfContents(; depth=4)

# ╔═╡ 0943c748-9fee-401e-851e-de2327d39706
md"""
## Grid
"""

# ╔═╡ d2b72544-78fb-4a0b-b169-aa50ebaeb31d
colpts = (; vertex=v, cell=c, edge=e);

# ╔═╡ 9abfbc35-bd37-4554-949c-27cd2bdfa1a7
md"""
__Grid__: $(@bind grid_t Select([:TriA, :TriB, :TriC]; default=:TriA))
"""

# ╔═╡ cfb78fd5-4435-4c0d-8cf2-1b5fffc2f91e
@variables f₀ g N² Ri le ϵ k l h a M² θU β

# ╔═╡ fb5e7874-c067-4769-aefa-260fd3eca00b
bflow = eady_background_flow(Val(grid_t), a; f₀, N², Ri, θU=0, β=0);

# ╔═╡ bd9551af-d36f-4b8c-99cf-02577b13b9b4
Ū = (z + 4000 * (1//2 + 0)) * -M²/f₀

# ╔═╡ 3b152922-e343-42a1-b863-38814a80b0d7
md"""
## Horizontal Momentum Transport
"""

# ╔═╡ 33e1ab03-9e1a-4014-84ed-1ad393f0445e
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

# ╔═╡ f9adc84b-12c3-4a1c-8400-6ede0108c9e9
md"""
__Scheme for the momentum transport__: $(@bind hmt_scheme select_hmt_scheme(grid_t))
"""

# ╔═╡ f71432b9-cb86-4f5e-a6ba-4ad443140292
ϕ = compute_phases(k, l, a);

# ╔═╡ 583e619a-76da-4f63-8377-b5eac40c96af
md"""
### Fourier Symbols
"""

# ╔═╡ 52def68d-5f55-48f8-8e9c-e2f387800bbc
md"""
#### Small wavenumber approximation
"""

# ╔═╡ 990c73d7-0d2a-457c-81c7-88a1f80cb44a
md"""
## Horizontal Scalar Transport
"""

# ╔═╡ 27ec48bb-d7bb-4b70-9932-08f0e0926504
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

# ╔═╡ e9dd8811-c91e-45e4-8ced-afd7391f248a
md"""
__Discretization scheme for the momentum transport__: 
$(@bind hst_scheme select_hst_scheme(grid_t))
"""

# ╔═╡ 698f9004-f56b-41c0-b6a6-b55c505fb1a8
md"""
### Fourier Symbols
"""

# ╔═╡ 7d55ea61-e0f1-4bd9-a36c-22857ad145e7
html"""<hr>"""

# ╔═╡ 8eea87eb-048e-4436-9d20-3bcc9e76c7b7
md"""
## Appendix
"""

# ╔═╡ f0f8c55d-82fa-4cf4-9678-555c1222e615
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ 9b8a0e33-2422-48e6-93d4-0313c557abe3
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

# ╔═╡ e322d734-1c1a-4ce1-a147-9f8420c64203
getflow(grid_t) =
	if grid_t == :TriA
		TriAFlow
	elseif grid_t == :TriB
		TriBFlow
	elseif grid_t == :TriC
		TriCFlow
	end

# ╔═╡ 8ae523f0-1eeb-4d62-b089-502c52dd1277
begin
	flow_t = getflow(grid_t)
	if grid_t == :TriC
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
    dflow = flow_t{nS}(; u⃗ = du⃗, w = dw, ∫∇ᵀu⃗dz = ∫∇ᵀdu⃗dz, b = db, p = dp, η = dη)
end;

# ╔═╡ e18fbb3e-02cb-4add-ab00-9e0cc7ae935f
floataside(
	@bind inputs PlutoUI.combine() do Child
	u⃗inputs = [
		md"""
		iH: $(Child("iH", Select(1:dims(colpt_type(flow_t, :u⃗)))))
		""",
		md"""
		jH: $(Child("jH", Select(1:dims(colpt_type(flow_t, :u⃗)))))
		"""
	]
	if colpt_type(flow_t, :u⃗) ≠ :edge
		push!(u⃗inputs, md"""
			  iTH: $(Child("iTH", Select(1:2, default=1)))
			  """)
		push!(u⃗inputs, md"""
			  jTH: $(Child("jTH", Select(1:2, default=1)))
			  """)
	end
	binputs = [
		md"""
		iH: $(Child("biH", Select(1:dims(colpt_type(flow_t, :b)))))
		""",
		md"""
		jH: $(Child("bjH", Select(1:dims(colpt_type(flow_t, :b)))))
		"""
	]	
	md"""
	__u input__:
	$(u⃗inputs)
	__b input__:
	$(binputs)
	"""
end; top=500)

# ╔═╡ 3782fb16-26b7-4edd-8591-ef119b1cabef
getflowft(grid_t) =
	if grid_t == :TriA
		TriAFlowFT
	elseif grid_t == :TriB
		TriBFlowFT
	elseif grid_t == :TriC
		TriCFlowFT
	end

# ╔═╡ 780abc6b-d9ab-40d1-9a01-d12b1d3fc3ae
begin
	flowft_t = getflowft(grid_t)
	if grid_t == :TriC
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
    fflow = flowft_t(; u⃗=u⃗̂, w=ŵ, ∫∇ᵀu⃗dz=∫∇ᵀu⃗̂dz, b=b̂, p=p̂, η=η̂)
end;

# ╔═╡ d424025e-3efd-4c4b-8d97-7369ff3618da
function gethmt(grid_t, hmt_scheme)
	if grid_t == :TriA
		if hmt_scheme == :standard
        	TriA.u⃗ᵀ∇
    	else
        	throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
    	end
	elseif grid_t == :TriB
		 if hmt_scheme == :asc
	        TriB.u⃗ᵀ∇_asc
	    elseif hmt_scheme == :avi
	        TriB.u⃗ᵀ∇_avi
	    elseif hmt_scheme == :fdcre
	        TriB.u⃗ᵀ∇_fdcre
	    #elseif hmt_scheme == :fdv
	    #    TriB.∇ᵀcc(cout, cin, u⃗, u⃗) - evalat(cout, cin, u⃗) * TriB.av_cv(cout, v, TriB.∇ᵀvc(v, cin, u⃗))
	    else
	        throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
	    end
	elseif grid_t == :TriC
		if hmt_scheme == :ICON
	        TriC.u⃗ᵀ∇
	    else
	        throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
	    end
	end
end

# ╔═╡ a64854b8-b161-48d2-87da-ebed3e08671c
hmt = let
	(; iH) = inputs
	cp = colpts[colpt_type(flow_t, :u⃗)]
	pu⃗ = grid_t == :TriC ? bflow.u⃗[cp[1], cp[2], cp[3]] .+ ϵ*dflow.u⃗[cp[1], cp[2], cp[3]] : [bflow.u⃗[iTH, cp[1], cp[2], cp[3]] .+ ϵ*dflow.u⃗[iTH, cp[1], cp[2], cp[3]] for iTH=1:2]
	u⃗ᵀ∇ = gethmt(grid_t, hmt_scheme)
	u⃗ᵀ∇(colptidx(0,0,iH,Val(colpt_type(flow_t, :u⃗))), cp, pu⃗, pu⃗)
end;

# ╔═╡ 18a2e213-d9a2-4094-b03a-f57df7c16e82
lhmt = if grid_t == :TriC
	expand(taylor_coeff(hmt, ϵ, 1))
else
	[expand(taylor_coeff(hmt[iTH], ϵ, 1)) for iTH=1:2]
end;

# ╔═╡ ec914f25-22a9-486b-b893-a5484cf109f7
fhmt = 
	let
		(; iH) = inputs
		if grid_t == :TriC
			fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhmt; fflow, dflow, ϕ)
		else
			[fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhmt[iTH]; fflow, dflow, ϕ) for iTH=1:2]
	end
end;

# ╔═╡ 423c5bfc-d6e2-46cc-a29f-05e3d9d3c84b
function gethst(grid_t, hst_scheme)
	if grid_t == :TriA
		if hst_scheme == :low
        	TriA.u⃗∇ᵀ
    	else
        	throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
    	end
	elseif grid_t == :TriB
		 if hst_scheme == :low
	        TriB.u⃗∇ᵀ_low
	    elseif hst_scheme == :high
	        TriB.u⃗∇ᵀ
	    else
	        throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
	    end
	elseif grid_t == :TriC
		if hst_scheme == :low
	        TriC.u⃗∇ᵀ
	    else
	        throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
	    end
	end
end

# ╔═╡ f80621ff-53de-4b21-b314-b8d4bbdbcc61
hst = let
	(; biH) = inputs
	cpu⃗ = colpts[colpt_type(flow_t, :u⃗)]
	cpb = colpts[colpt_type(flow_t, :b)]
	pu⃗ = grid_t == :TriC ? bflow.u⃗[cpu⃗[1], cpu⃗[2], cpu⃗[3]] .+ ϵ*dflow.u⃗[cpu⃗[1], cpu⃗[2], cpu⃗[3]] : [bflow.u⃗[iTH, cpu⃗[1], cpu⃗[2], cpu⃗[3]] .+ ϵ*dflow.u⃗[iTH, cpu⃗[1], cpu⃗[2], cpu⃗[3]] for iTH=1:2]
	pb =  bflow.b[cpb[1], cpb[2], cpb[3]] .+ ϵ*dflow.b[cpb[1], cpb[2], cpb[3]]
	u⃗∇ᵀ = gethst(grid_t, hst_scheme)
	u⃗∇ᵀ(colptidx(0,0,biH,Val(colpt_type(flow_t, :b))), cpu⃗, cpb, pu⃗, pb)
end;

# ╔═╡ 534e0d3c-8453-4096-8aaf-354691a04826
lhst = expand(taylor_coeff(expand(hst), ϵ, 1));

# ╔═╡ a30cee69-cce3-4044-900e-a3a6188b4653
fhst = let
	(; biH) = inputs
	fourier_transform_expression(biH, colpt_type(flow_t, :b), lhst; fflow, dflow, ϕ)
end;

# ╔═╡ ed74cf26-7ae6-43f6-9303-09c46b0fae0a
rtrig = let
	rcos = let
    	x = Symbolics.variable(:x)
    	Symbolics.@rule cos(~x)=>substitute(taylor(cos(x),x,0:5), Dict([x=>~x]))
	end
	rsin = let
    	x = Symbolics.variable(:x)
    	Symbolics.@rule sin(~x)=>substitute(taylor(sin(x),x,0,0:5), Dict([x=>~x]))
	end
	SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(SymbolicUtils.RestartedChain([rcos, rsin])))
end

# ╔═╡ c8c98ec3-094e-4749-a611-2f1675f30fcb
sqrt3subs = Dict(
	GridOperatorAnalysis.sqrt3^2=>3, GridOperatorAnalysis.sqrt3^3=>3*GridOperatorAnalysis.sqrt3, GridOperatorAnalysis.sqrt3^4=>9, GridOperatorAnalysis.sqrt3^5=>9*GridOperatorAnalysis.sqrt3, GridOperatorAnalysis.sqrt3^6=>27,
	GridOperatorAnalysis.sqrt3^7=>27*GridOperatorAnalysis.sqrt3,
	GridOperatorAnalysis.sqrt3^8=>81,
);

# ╔═╡ 2e0a237f-fe84-4fd5-a26f-76d7976a16ef
if colpt_type(flow_t, :u⃗) == :edge
	(; jH) = inputs
	_fhmt = expand(fhmt)
	_fhmt = substitute(_fhmt, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
	_fhmt = expand(_fhmt)
	simplify(taylor_coeff(_fhmt, fflow.u⃗[jH], 1))
else
	(; jH, iTH, jTH) = inputs
	_fhmt = expand(fhmt[iTH])
	_fhmt = substitute(_fhmt, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
	_fhmt = expand(_fhmt)
	simplify(taylor_coeff(_fhmt, fflow.u⃗[jTH, jH], 1))
end

# ╔═╡ 4c1b9413-298c-42f5-aa72-c53c84f0f19d
fmu = let
	if colpt_type(flow_t, :u⃗) == :edge
		_fhmt = simplify(fhmt; rewriter=rtrig)
		_fhmt = expand(_fhmt)
		_fhmt = substitute(_fhmt, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
		_fhmt = taylor_coeff(_fhmt, k, 1)
		_fhmt = expand(_fhmt)
		simplify(taylor_coeff(_fhmt, fflow.u⃗[jH], 1))
	else
		_fhmt = simplify(fhmt[iTH]; rewriter=rtrig)
		_fhmt = expand(_fhmt)
		_fhmt = substitute(_fhmt, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
		_fhmt = taylor_coeff(_fhmt, k, 1)
		_fhmt = expand(_fhmt)
		simplify(taylor_coeff(_fhmt, fflow.u⃗[jTH,jH], 1))
	end
end

# ╔═╡ bf388183-ad70-4a97-b9ac-2c0f96ef87ef
simplify(fmu / Ū)

# ╔═╡ 38060711-bb39-4437-bcf9-a412a3eb4c08
let
	_fhst = expand(fhst)
	_fhst = substitute(_fhst, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
	_fhst = expand(_fhst)
	if colpt_type(flow_t, :u⃗) == :edge
		(; jH, bjH) = inputs
		(
			simplify(taylor_coeff(_fhst, fflow.u⃗[jH], 1)),
			simplify(taylor_coeff(_fhst, fflow.b[bjH], 1))
		)
	else
		(; jH, jTH, bjH) = inputs
		(
			simplify(taylor_coeff(_fhst, fflow.u⃗[jTH, jH], 1)),
			simplify(taylor_coeff(_fhst, fflow.b[bjH], 1))
		)
	end
end

# ╔═╡ 60f00b5c-e56f-45cb-b032-cd60d6384558
fsu, fsb = let
	(; bjH) = inputs
	_fhst = simplify(fhst; rewriter=rtrig)
	_fhst = expand(_fhst)
	_fhst = substitute(_fhst, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
	_fhst = taylor_coeff(_fhst, k, 1)
	_fhst = expand(_fhst)
	fsb   = taylor_coeff(_fhst, fflow.b[bjH], 1)
	fsb   = expand(simplify(fsb))
	fsu = if colpt_type(flow_t, :u⃗) == :edge
		(; jH) = inputs
		fs = taylor_coeff(_fhst, fflow.u⃗[jH], 1)
		expand(simplify(fs))
	else
		(; jTH, jH) = inputs
		fs = taylor_coeff(_fhst, fflow.u⃗[jTH,jH], 1)
		expand(simplify(fs))
	end
	(fsu, fsb)
end

# ╔═╡ 81781a43-27ba-4920-8323-473cf3feb976
expand(simplify(fsb / Ū))

# ╔═╡ Cell order:
# ╟─b9cca8c1-0d0f-48f9-b21b-8d6df2cb77aa
# ╠═f59a5438-8d5e-11f0-13fa-d9c703e5f87f
# ╟─59056484-9c6b-48ad-b741-2bc294d5cc6f
# ╟─0943c748-9fee-401e-851e-de2327d39706
# ╠═d2b72544-78fb-4a0b-b169-aa50ebaeb31d
# ╟─9abfbc35-bd37-4554-949c-27cd2bdfa1a7
# ╠═cfb78fd5-4435-4c0d-8cf2-1b5fffc2f91e
# ╠═fb5e7874-c067-4769-aefa-260fd3eca00b
# ╟─bd9551af-d36f-4b8c-99cf-02577b13b9b4
# ╠═8ae523f0-1eeb-4d62-b089-502c52dd1277
# ╠═780abc6b-d9ab-40d1-9a01-d12b1d3fc3ae
# ╟─3b152922-e343-42a1-b863-38814a80b0d7
# ╟─33e1ab03-9e1a-4014-84ed-1ad393f0445e
# ╟─f9adc84b-12c3-4a1c-8400-6ede0108c9e9
# ╠═a64854b8-b161-48d2-87da-ebed3e08671c
# ╠═18a2e213-d9a2-4094-b03a-f57df7c16e82
# ╠═f71432b9-cb86-4f5e-a6ba-4ad443140292
# ╠═ec914f25-22a9-486b-b893-a5484cf109f7
# ╟─583e619a-76da-4f63-8377-b5eac40c96af
# ╟─2e0a237f-fe84-4fd5-a26f-76d7976a16ef
# ╟─52def68d-5f55-48f8-8e9c-e2f387800bbc
# ╟─4c1b9413-298c-42f5-aa72-c53c84f0f19d
# ╠═bf388183-ad70-4a97-b9ac-2c0f96ef87ef
# ╟─990c73d7-0d2a-457c-81c7-88a1f80cb44a
# ╟─27ec48bb-d7bb-4b70-9932-08f0e0926504
# ╟─e9dd8811-c91e-45e4-8ced-afd7391f248a
# ╠═f80621ff-53de-4b21-b314-b8d4bbdbcc61
# ╠═534e0d3c-8453-4096-8aaf-354691a04826
# ╠═a30cee69-cce3-4044-900e-a3a6188b4653
# ╟─698f9004-f56b-41c0-b6a6-b55c505fb1a8
# ╟─38060711-bb39-4437-bcf9-a412a3eb4c08
# ╠═60f00b5c-e56f-45cb-b032-cd60d6384558
# ╠═81781a43-27ba-4920-8323-473cf3feb976
# ╟─7d55ea61-e0f1-4bd9-a36c-22857ad145e7
# ╟─8eea87eb-048e-4436-9d20-3bcc9e76c7b7
# ╟─f0f8c55d-82fa-4cf4-9678-555c1222e615
# ╟─9b8a0e33-2422-48e6-93d4-0313c557abe3
# ╟─e18fbb3e-02cb-4add-ab00-9e0cc7ae935f
# ╟─e322d734-1c1a-4ce1-a147-9f8420c64203
# ╟─3782fb16-26b7-4edd-8591-ef119b1cabef
# ╟─d424025e-3efd-4c4b-8d97-7369ff3618da
# ╟─423c5bfc-d6e2-46cc-a29f-05e3d9d3c84b
# ╟─ed74cf26-7ae6-43f6-9303-09c46b0fae0a
# ╠═c8c98ec3-094e-4749-a611-2f1675f30fcb
