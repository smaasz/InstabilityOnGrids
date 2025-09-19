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
	using Groebner
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
@variables f₀ g N² Ri le ϵ k l h a M² β θU

# ╔═╡ c6283dad-ae0a-4b8c-8c93-5ae02668a191
rewriter = let
	r = Symbolics.@acrule ~α * ~(~x) + ~β * ~(~x) => *(~α + ~β, ~(~x)...)
	rewriter = SymbolicUtils.Prewalk(SymbolicUtils.PassThrough(r))
end

# ╔═╡ fb5e7874-c067-4769-aefa-260fd3eca00b
bflow = eady_background_flow(Val(grid_t), a; f₀, N², Ri, θU, β);

# ╔═╡ 1b94a213-9311-4dad-bc3e-0ac2dcb71e14
Ū = (z + 4000 * (1//2 + β)) * -M²/f₀

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

# ╔═╡ 0f1cbe99-8b11-426d-8d28-77b6dffef277
md"""
Small wavenumber approximation: $(@bind doapprox CheckBox(default=false))
"""

# ╔═╡ 52def68d-5f55-48f8-8e9c-e2f387800bbc
md"""
#### Small wavenumber approximation
"""

# ╔═╡ 19e3a83c-d29a-408d-b965-74bcf6aa09e5


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
            	width: 200px;
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
end; top=300)

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
    fflow = flowft_t{Symbolics.Num}(; u⃗=u⃗̂, w=ŵ, ∫∇ᵀu⃗dz=∫∇ᵀu⃗̂dz, b=b̂, p=p̂, η=η̂)
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
	    elseif hmt_scheme == :fdv
	        TriB.u⃗ᵀ∇_fdv
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
begin
	hmts = []
	for iH=1:dims(colpt_type(flow_t, :u⃗))
		cp = colpts[colpt_type(flow_t, :u⃗)]
		pu⃗ = grid_t == :TriC ? bflow.u⃗[cp[1], cp[2], cp[3]] .+ ϵ*dflow.u⃗[cp[1], cp[2], cp[3]] : [bflow.u⃗[iTH, cp[1], cp[2], cp[3]] .+ ϵ*dflow.u⃗[iTH, cp[1], cp[2], cp[3]] for iTH=1:2]
		u⃗ᵀ∇ = gethmt(grid_t, hmt_scheme)
		push!(hmts, u⃗ᵀ∇(colptidx(0,0,iH,Val(colpt_type(flow_t, :u⃗))), cp, pu⃗, pu⃗))
	end
	hmt = hmts[inputs.iH]
end;

# ╔═╡ 50fc9497-99e4-44d7-b25f-ffccd14b145a
lhmts = if grid_t == :TriC
	[expand(taylor_coeff(hmts[iH], ϵ, 1)) for iH=1:dims(colpt_type(flow_t, :u⃗))]
else
	[[expand(taylor_coeff(hmts[iH][iTH], ϵ, 1)) for iTH=1:2] for iH=1:dims(colpt_type(flow_t, :u⃗))]
end;

# ╔═╡ 26b13d30-d8cf-4ba4-9f84-a591bb830630
begin
	fhmts = []
	for iH=1:dims(colpt_type(flow_t, :u⃗))
		if grid_t == :TriC
			fhmt = fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhmts[iH]; fflow, dflow, ϕ)
			push!(fhmts, fhmt)
		else
			fhmt = [fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhmts[iH][iTH]; fflow, dflow, ϕ) for iTH=1:2]
			push!(fhmts, fhmt)
		end
	end		
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
	hst = u⃗∇ᵀ(colptidx(0,0,biH,Val(colpt_type(flow_t, :b))), cpu⃗, cpb, pu⃗, pb)
	substitute(hst, GridOperatorAnalysis.sqrt3_subs)
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
	function p(x)
		 !isequal(x, θU)
	end
	rcos = let
    	x = Symbolics.variable(:x)
    	Symbolics.@rule cos(~x::p) => substitute(taylor(cos(x),x,0:6), Dict([x=>~x]))
	end
	rsin = let
    	x = Symbolics.variable(:x)
    	Symbolics.@rule sin(~x::p)=>substitute(taylor(sin(x),x,0,0:6), Dict([x=>~x]))
	end
	SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(SymbolicUtils.RestartedChain([rcos, rsin])))
end

# ╔═╡ 7b864ab1-b61c-481c-849a-9824c9abc984
rpyt = let
	r = Symbolics.@acrule sin(~x)^2 + cos(~x)^2 => one(~x)
	@show r(cos(l)^2+sin(l)^2)
	SymbolicUtils.Prewalk(SymbolicUtils.PassThrough(r))
end

# ╔═╡ c8c98ec3-094e-4749-a611-2f1675f30fcb
sqrt3subs = Dict(
	GridOperatorAnalysis.sqrt3^2=>3, GridOperatorAnalysis.sqrt3^3=>3*GridOperatorAnalysis.sqrt3, GridOperatorAnalysis.sqrt3^4=>9, GridOperatorAnalysis.sqrt3^5=>9*GridOperatorAnalysis.sqrt3, GridOperatorAnalysis.sqrt3^6=>27,
	GridOperatorAnalysis.sqrt3^7=>27*GridOperatorAnalysis.sqrt3,
	GridOperatorAnalysis.sqrt3^8=>81,
);

# ╔═╡ 192a4172-a2bb-416b-83df-194a090b093a
begin
	sfhmts = let
		colpt_t = colpt_type(flow_t, :u⃗)
		d = dims(colpt_t)
		if colpt_t == :edge
			zeros(Complex{Symbolics.Num}, d, d)
		else
			zeros(Complex{Symbolics.Num}, d, 2, d, 2)
		end
	end
	for iH=1:dims(colpt_type(flow_t, :u⃗))
		fhmt = let
			fhmt = expand.(fhmts[iH])
			if doapprox
				fhmt = expand.(simplify.(fhmt; rewriter=rtrig))
			else
				fhmt = substitute.(fhmt, Ref(Dict(l => h/a * 2/GridOperatorAnalysis.sqrt3 * l)))
			end
			fhmt = substitute.(fhmt, Ref(Dict(le=>a, √(f₀^2*N²/Ri) => M² ))) 
			fhmt = expand.(fhmt)
			fhmt = substitute.(fhmt, Ref(sqrt3subs))
		end
		for iTH=1:2
			for jH=1:dims(colpt_type(flow_t, :u⃗))
				for jTH=1:2
					u = fflow.u⃗[jTH, jH]
					sfhmt = taylor_coeff(fhmt[iTH], u, 1)
					sfhmts[iH, iTH, jH, jTH] = simplify(sfhmt; expand=true)
				end
			end
		end
	end
end

# ╔═╡ 88ad62f8-97f2-43b1-92bc-34ef51cf99f5
begin
	@variables k̃, l̃
	vals = let
		A = sfhmts[:,1,:,1] ./ Ū
		A = substitute.(A, Ref(Dict(cos(θU) => k̃/k, sin(θU) => l̃/l)))
		if doapprox
			A = substitute.(A, Ref(Dict(a=>2/GridOperatorAnalysis.sqrt3 * h)))
			A = simplify.(A; expand=true)
			A = substitute.(A, Ref(sqrt3subs))
		else
			A = substitute.(A, Ref(Dict(GridOperatorAnalysis.sqrt3=> 2 * h/a)))
		end
		A = simplify.(A; expand=true)
		A = k̃ * simplify.(taylor_coeff.(A, k̃, 1)) + l̃ * simplify.(taylor_coeff.(A, l̃, 1))
	end
end

# ╔═╡ 2e0a237f-fe84-4fd5-a26f-76d7976a16ef
if colpt_type(flow_t, :u⃗) == :edge
	(; jH) = inputs
	_fhmt = expand(fhmt)
	_fhmt = substitute(_fhmt, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
	_fhmt = expand(_fhmt)
	simplify(expand(taylor_coeff(_fhmt, fflow.u⃗[jH], 1)))
else
	(; jH, iTH, jTH) = inputs
	_fhmt = expand(fhmt[iTH])
	_fhmt = substitute(_fhmt, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
	_fhmt = expand(_fhmt)
	simplify(expand(taylor_coeff(_fhmt, fflow.u⃗[jTH, jH], 1)))
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
simplify(fmu / Ū / cos(θU))

# ╔═╡ 38060711-bb39-4437-bcf9-a412a3eb4c08
fsu_, fsb_ = let
	_fhst = expand(fhst)
	_fhst = substitute(_fhst, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
	_fhst = expand(_fhst)
	if colpt_type(flow_t, :u⃗) == :edge
		(; jH, bjH) = inputs
		(
			simplify(taylor_coeff(_fhst, fflow.u⃗[jH], 1); expand=true),
			simplify(taylor_coeff(_fhst, fflow.b[bjH], 1); expand=true)
		)
	else
		(; jH, jTH, bjH) = inputs
		(
			simplify(taylor_coeff(_fhst, fflow.u⃗[jTH, jH], 1); expand=true),
			simplify(taylor_coeff(_fhst, fflow.b[bjH], 1); expand=true)
		)
	end
end

# ╔═╡ 60f00b5c-e56f-45cb-b032-cd60d6384558
fsu, fsb = let
	fsus = []
	fsbs = []
	for var in [k, l]
		(; bjH) = inputs
		_fhst = simplify(fhst; rewriter=rtrig)
		_fhst = expand(_fhst)
		_fhst = substitute(_fhst, merge(sqrt3subs, Dict(le=>a, √(f₀^2*N²/Ri) => M²)))
		_fhst = taylor_coeff(_fhst, var, 1)
		_fhst = simplify(_fhst; expand=true)
		fsb   = taylor_coeff(_fhst, fflow.b[bjH], 1)
		fsb   = simplify(fsb; expand=true)
		push!(fsbs, fsb)
		fsu = if colpt_type(flow_t, :u⃗) == :edge
			(; jH) = inputs
			fs = taylor_coeff(_fhst, fflow.u⃗[jH], 1)
			simplify(fs; expand=true)
		else
			(; jTH, jH) = inputs
			fs = taylor_coeff(_fhst, fflow.u⃗[jTH,jH], 1)
			simplify(fs; expand=true)
		end
		push!(fsus, fsu)
	end
	#(simplify(cos(θU) * fsus[1] + sin(θU) * fsus[2]), cos(θU) * fsbs[1] + sin(θU) * fsbs[2])
	(fsus[1], fsbs[1])
end

# ╔═╡ 81781a43-27ba-4920-8323-473cf3feb976
simplify(simplify(fsb / Ū / cos(θU)))

# ╔═╡ f0e541cf-b2f5-40f5-8379-b1f8ef1a4fe8
let
	@variables a[1:300]
	r = Symbolics.@rule ~α * ~(~x) + ~β * ~(~x) => *(~α + ~β, ~(~x)...)
	rewriter = SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(r))
	t = +([a[i] * a[i+1] * a[i+2] * a[i+3] * a[i+4] for i=1:5:300]...)
	simplify(t; rewriter)
end

# ╔═╡ 97fca3d8-88cb-4cc2-851d-be78a5c761d2
let
	sys = GridOperatorAnalysis.TriALCCBoussinesqSys()
	push!(sys.momentum_transport_eq[1, 1], k)
	push!(sys.momentum_transport_eq[2, 1], l)
	(; state, momentum_transport_eq) = sys
	momentum_transport_eq
end

# ╔═╡ ab256f38-0986-468b-a22b-2eb07c88542e
sys = let
	import GridOperatorAnalysis: construct_lcc_boussinesq_sys
	construct_lcc_boussinesq_sys(TriAFlow, le, dflow; f₀, g, 𝕂ᵘ=0, 𝕂ᵇ=0)
end;

# ╔═╡ ace860ae-06c1-4650-b615-bc7edf7630e2
@show sys.momentum_transport_eq[1,1]

# ╔═╡ d94dff98-3128-4911-82e4-faccae824231
let
	@variables z::Complex w
	import LinearAlgebra: I
	function Base.sqrt(z::Complex{Symbolics.Num})
		@variables x
		vars = [Symbolics.get_variables(real(z)); Symbolics.get_variables(imag(z))]
		Dx = Symbolics.Differential(x)
		z̄ = Symbolics.value(substitute(z, Dict(vars .=> 0)))
		cs = [substitute(Symbolics.expand_derivatives(1/factorial(k)*(Dx^k)(sqrt(x))), Dict(x=>z̄)) for k=0:1]
		sqrtz = Symbolics.series(cs, z-z̄)
		Symbolics.simplify(sqrtz; expand=true)
	end
	function quadratic(c, b, a)
        discr = Symbolics.simplify(b^2 - 4*a*c; expand=true)
        Symbolics.simplify.([(-b + sqrt(discr))/(2*a), (-b - sqrt(discr))/(2*a)]; expand=true )
    end
	function symbolic_eigenvals(A::Matrix{<:Number})
	    @variables λ μ::Complex # eigenvalue  
	    # find eigenvalues first
	    p = Symbolics.expand(Symbolics.det(complex.(λ*I- A))) # polynomial to solve
		if Symbolics.degree(real(p), λ) ≤ 1 && Symbolics.degree(imag(p), λ) ≤ 1
			cs = Symbolics.taylor_coeff(p, λ, 0:1)
			p = Symbolics.series(cs, μ)
			GridOperatorAnalysis.symbolic_linear_solve(p ~ 0, μ)
		else
	    	quadratic(Symbolics.taylor_coeff(p, λ, 0:2)...) # solve polynomial
		end
	end
	#A = Complex{Symbolics.Num}[0; 1+im*(1-k);; -1+im*(1-k); 0]
	_vals = symbolic_eigenvals(taylor_coeff.(vals, k̃, 1))
	#_vals = symbolic_eigenvals(A)
	_vals = Symbolics.simplify.(_vals; expand=true)
end

# ╔═╡ Cell order:
# ╟─b9cca8c1-0d0f-48f9-b21b-8d6df2cb77aa
# ╠═f59a5438-8d5e-11f0-13fa-d9c703e5f87f
# ╟─59056484-9c6b-48ad-b741-2bc294d5cc6f
# ╟─0943c748-9fee-401e-851e-de2327d39706
# ╠═d2b72544-78fb-4a0b-b169-aa50ebaeb31d
# ╟─9abfbc35-bd37-4554-949c-27cd2bdfa1a7
# ╠═cfb78fd5-4435-4c0d-8cf2-1b5fffc2f91e
# ╟─c6283dad-ae0a-4b8c-8c93-5ae02668a191
# ╠═fb5e7874-c067-4769-aefa-260fd3eca00b
# ╠═1b94a213-9311-4dad-bc3e-0ac2dcb71e14
# ╠═8ae523f0-1eeb-4d62-b089-502c52dd1277
# ╠═780abc6b-d9ab-40d1-9a01-d12b1d3fc3ae
# ╟─3b152922-e343-42a1-b863-38814a80b0d7
# ╟─33e1ab03-9e1a-4014-84ed-1ad393f0445e
# ╟─f9adc84b-12c3-4a1c-8400-6ede0108c9e9
# ╠═a64854b8-b161-48d2-87da-ebed3e08671c
# ╠═50fc9497-99e4-44d7-b25f-ffccd14b145a
# ╠═18a2e213-d9a2-4094-b03a-f57df7c16e82
# ╠═f71432b9-cb86-4f5e-a6ba-4ad443140292
# ╠═26b13d30-d8cf-4ba4-9f84-a591bb830630
# ╠═ec914f25-22a9-486b-b893-a5484cf109f7
# ╟─583e619a-76da-4f63-8377-b5eac40c96af
# ╟─0f1cbe99-8b11-426d-8d28-77b6dffef277
# ╠═192a4172-a2bb-416b-83df-194a090b093a
# ╠═88ad62f8-97f2-43b1-92bc-34ef51cf99f5
# ╠═2e0a237f-fe84-4fd5-a26f-76d7976a16ef
# ╟─52def68d-5f55-48f8-8e9c-e2f387800bbc
# ╠═4c1b9413-298c-42f5-aa72-c53c84f0f19d
# ╠═bf388183-ad70-4a97-b9ac-2c0f96ef87ef
# ╠═19e3a83c-d29a-408d-b965-74bcf6aa09e5
# ╟─990c73d7-0d2a-457c-81c7-88a1f80cb44a
# ╟─27ec48bb-d7bb-4b70-9932-08f0e0926504
# ╟─e9dd8811-c91e-45e4-8ced-afd7391f248a
# ╠═f80621ff-53de-4b21-b314-b8d4bbdbcc61
# ╠═534e0d3c-8453-4096-8aaf-354691a04826
# ╠═a30cee69-cce3-4044-900e-a3a6188b4653
# ╟─698f9004-f56b-41c0-b6a6-b55c505fb1a8
# ╠═38060711-bb39-4437-bcf9-a412a3eb4c08
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
# ╠═ed74cf26-7ae6-43f6-9303-09c46b0fae0a
# ╠═7b864ab1-b61c-481c-849a-9824c9abc984
# ╠═c8c98ec3-094e-4749-a611-2f1675f30fcb
# ╠═f0e541cf-b2f5-40f5-8379-b1f8ef1a4fe8
# ╠═97fca3d8-88cb-4cc2-851d-be78a5c761d2
# ╠═ab256f38-0986-468b-a22b-2eb07c88542e
# ╠═ace860ae-06c1-4650-b615-bc7edf7630e2
# ╠═d94dff98-3128-4911-82e4-faccae824231
