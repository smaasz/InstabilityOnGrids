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
	using CairoMakie
	import GridOperatorAnalysis: eady_background_flow, bb, e, c, v, nS, dims, t, z, sqrt3
	import GridOperatorAnalysis: TriAFlow, TriAFlowFT, TriBFlow, TriBFlowFT, TriCFlow, TriCFlowFT, HexCFlow, HexCFlowFT
	import GridOperatorAnalysis: colpt_type, colptidx, compute_phases, sqrt3_subs
	import GridOperatorAnalysis: fourier_transform_expression
	import GridOperatorAnalysis.TriA
	import GridOperatorAnalysis.TriB
	import GridOperatorAnalysis.TriC
	import GridOperatorAnalysis.HexC
	import GridOperatorAnalysis
	import LinearAlgebra: eigen
	import OffsetArrays: no_offset_view
	import Symbolics: substitute, @variables, taylor_coeff, simplify, coeff, taylor, expand, det
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
__Grid__: $(@bind grid_t Select([:TriA, :TriB, :TriC, :HexC]; default=:TriA))
"""

# ╔═╡ cfb78fd5-4435-4c0d-8cf2-1b5fffc2f91e
@variables f₀ g N² Ri le ϵ k l h a M² β θU

# ╔═╡ 428cc192-03d0-497c-aee7-f7015c4c09d5
md"""
#### Background flow

The baroclinically unstable flow with a sloping buoyancy profile in across-front direction (first studied by Eady) is projected onto the respective discrete spaces.
The flow is parametrized by:
- the Richardson number $\mathrm{Ri}$, 
- the Brunt-Väisälä frequency $\mathrm{N^2}$, 
- the Coriolis frequency $f_0$ on the f-plane, 
- the flow direction $θ_U$, 
- and the flow shift $\beta$ (in $HM^2/f_0$).

The flow on the mesh is restriced to a small neighborhood around the lattices' origins.
"""

# ╔═╡ fb5e7874-c067-4769-aefa-260fd3eca00b
bflow = eady_background_flow(Val(grid_t), a; f₀, N², Ri, θU, β);

# ╔═╡ a697761a-3f01-46ac-ae31-4665bc7b6e47
md"""
__Vertical velocity profile in the flow direction__
"""

# ╔═╡ 2663542e-c719-42f0-bf00-790e8aa211fe
#substitute(bflow.u⃗[1,0,0,1] / cos(θU), Dict(√(f₀^2*N²/Ri) => M²))

# ╔═╡ 1b94a213-9311-4dad-bc3e-0ac2dcb71e14
Ū = (z + 4000 * (1//2 + β)) * -M²/f₀

# ╔═╡ 430cca4e-08b6-4d57-a4dc-422dcf8cc92e
md"""
#### Pertubation Flow

The background flow is perturbed by another flow scaled by the asymptotic parameter $\epsilon$.
"""

# ╔═╡ 77eb6488-c26a-416a-839b-1e617f0cdd50
md"""
#### Flow in Wavenumber Space
"""

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

# ╔═╡ 39031d65-863d-48cf-aeaf-520417eaf071
md"""
#### Linearization Procedure

First, evaluate the momentum transport operator at the perturbed background flow. To obtain the linearized expression in the perturbed flow, the nonlinear expression is (Taylor) expanded in $\epsilon$ and only the term that is of first-order in $\epsilon$ is retained. Since the background velocity is constant on horizontal planes, the linearized operator is also constant on horizontal planes.
"""

# ╔═╡ 286c65de-46f4-462d-903c-ad850bb1a3b1
md"""
#### Fourier Transform

Replacing the state variables' fields of the pertubation flow by plane waves with wavenumber $\vec{k} = \begin{bmatrix} k & l \end{bmatrix}^T$, the linearized discrete operator is transformed into a multiplication operator acting on the amplitudes of the respective plane wave fields. The multiplication operator is also called the __Fourier symbol__ of the linearized discrete operator.
"""

# ╔═╡ f71432b9-cb86-4f5e-a6ba-4ad443140292
ϕ = compute_phases(k, l, a);

# ╔═╡ 583e619a-76da-4f63-8377-b5eac40c96af
md"""
### Fourier Symbols
"""

# ╔═╡ 0f1cbe99-8b11-426d-8d28-77b6dffef277
md"""
Small wavenumber approximation: $(@bind doapprox CheckBox(default=true))
"""

# ╔═╡ 46afe77d-f904-40e9-92bb-cab83aeea51e
md"""
#### Eigenvalues of Fourier Symbols
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

# ╔═╡ 363d736c-c30b-4878-8f32-c066df479600
md"""
#### Linearization Procedure
"""

# ╔═╡ 386315ca-ef94-44f5-9df8-eb320be72b34
md"""
#### Fourier Transform
"""

# ╔═╡ 698f9004-f56b-41c0-b6a6-b55c505fb1a8
md"""
### Fourier Symbols
"""

# ╔═╡ 99c92534-3092-4995-8044-d8001aad8f5c
md"""
Small wavenumber approximation: $(@bind doapproxs CheckBox(default=true))
"""

# ╔═╡ ed799518-6bf1-4180-8206-0a3b6e51aed0
md"""
#### Advection of the buoyancy pertubation
"""

# ╔═╡ 05fb2f19-1b43-4e7a-84c5-1ea8a71b9284
md"""
#### Advection of the background buoyancy by the velocity pertubation
"""

# ╔═╡ 26fd8b72-c736-43db-9d32-f76cbe5dd901
md"""
#### Eigenvalues of Fourier Symbols
"""

# ╔═╡ bc991b63-8e0e-4e37-9369-15f2e11e195e
md"""
## Diffusion Operator
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

# ╔═╡ e18fbb3e-02cb-4add-ab00-9e0cc7ae935f
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ e322d734-1c1a-4ce1-a147-9f8420c64203
getflow(grid_t) =
	if grid_t == :TriA
		TriAFlow
	elseif grid_t == :TriB
		TriBFlow
	elseif grid_t == :TriC
		TriCFlow
	else
		HexCFlow
	end

# ╔═╡ 8ae523f0-1eeb-4d62-b089-502c52dd1277
begin
	flow_t = getflow(grid_t)
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
    dflow = flow_t{nS}(; u⃗ = du⃗, w = dw, ∫∇ᵀu⃗dz = ∫∇ᵀdu⃗dz, b = db, p = dp, η = dη)
end;

# ╔═╡ 7f20375f-6491-47ff-bcc0-fb77baa8bedf
hds = let
	cpb = colpts[colpt_type(flow_t, :b)]
	pb =  ϵ*dflow.b[cpb[1], cpb[2], cpb[3]] # bflow.b[cpb[1], cpb[2], cpb[3]] .+
	hds = []
	for biH=1:dims(colpt_type(flow_t, :b))
		hd = HexC.Δ(colptidx(0,0,biH,Val(colpt_type(flow_t, :b))), cpb, pb)
		hd = substitute(hd, sqrt3_subs)
		push!(hds, hd)
	end
	hds
end;

# ╔═╡ 1af2939c-6795-4472-ae44-8c9466a159ab
lhds = [
	expand(taylor_coeff(expand(hd), ϵ, 1)) for hd in hds
];

# ╔═╡ 3782fb16-26b7-4edd-8591-ef119b1cabef
getflowft(grid_t) =
	if grid_t == :TriA
		TriAFlowFT
	elseif grid_t == :TriB
		TriBFlowFT
	elseif grid_t == :TriC
		TriCFlowFT
	else
		HexCFlowFT
	end

# ╔═╡ 780abc6b-d9ab-40d1-9a01-d12b1d3fc3ae
begin
	flowft_t = getflowft(grid_t)
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
    fflow = flowft_t{Symbolics.Num}(; u⃗=u⃗̂, w=ŵ, ∫∇ᵀu⃗dz=∫∇ᵀu⃗̂dz, b=b̂, p=p̂, η=η̂)
end;

# ╔═╡ 4b2d9f6d-4438-436a-91cb-47fb66270841
fhds = [
	fourier_transform_expression(biH, colpt_type(flow_t, :b), lhds[biH]; fflow, dflow, ϕ) for biH=1:dims(colpt_type(flow_t, :b))
];

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
	else
		if hmt_scheme == :MPAS
	        HexC.u⃗ᵀ∇
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
		pu⃗ = colpt_type(flow_t, :u⃗) == :edge ? bflow.u⃗[cp[1], cp[2], cp[3]] .+ ϵ*dflow.u⃗[cp[1], cp[2], cp[3]] : [bflow.u⃗[iTH, cp[1], cp[2], cp[3]] .+ ϵ*dflow.u⃗[iTH, cp[1], cp[2], cp[3]] for iTH=1:2]
		u⃗ᵀ∇ = gethmt(grid_t, hmt_scheme)
		push!(hmts, u⃗ᵀ∇(colptidx(0,0,iH,Val(colpt_type(flow_t, :u⃗))), cp, pu⃗, pu⃗))
	end
end;

# ╔═╡ 50fc9497-99e4-44d7-b25f-ffccd14b145a
lhmts = if colpt_type(flow_t, :u⃗) == :edge
	[expand(taylor_coeff(hmts[iH], ϵ, 1)) for iH=1:dims(colpt_type(flow_t, :u⃗))]
else
	[[expand(taylor_coeff(hmts[iH][iTH], ϵ, 1)) for iTH=1:2] for iH=1:dims(colpt_type(flow_t, :u⃗))]
end;

# ╔═╡ 26b13d30-d8cf-4ba4-9f84-a591bb830630
begin
	fhmts = []
	for iH=1:dims(colpt_type(flow_t, :u⃗))
		if colpt_type(flow_t, :u⃗) == :edge
			fhmt = fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhmts[iH]; fflow, dflow, ϕ)
			push!(fhmts, fhmt)
		else
			fhmt = [fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhmts[iH][iTH]; fflow, dflow, ϕ) for iTH=1:2]
			push!(fhmts, fhmt)
		end
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
	else
		if hst_scheme == :low
	        HexC.u⃗∇ᵀ
	    else
	        throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
	    end
	end
end

# ╔═╡ f80621ff-53de-4b21-b314-b8d4bbdbcc61
hsts = let
	cpu⃗ = colpts[colpt_type(flow_t, :u⃗)]
	cpb = colpts[colpt_type(flow_t, :b)]
	pu⃗ = colpt_type(flow_t, :u⃗) == :edge ? bflow.u⃗[cpu⃗[1], cpu⃗[2], cpu⃗[3]] .+ ϵ*dflow.u⃗[cpu⃗[1], cpu⃗[2], cpu⃗[3]] : [bflow.u⃗[iTH, cpu⃗[1], cpu⃗[2], cpu⃗[3]] .+ ϵ*dflow.u⃗[iTH, cpu⃗[1], cpu⃗[2], cpu⃗[3]] for iTH=1:2]
	pb =  bflow.b[cpb[1], cpb[2], cpb[3]] .+ ϵ*dflow.b[cpb[1], cpb[2], cpb[3]]
	u⃗∇ᵀ = gethst(grid_t, hst_scheme)
	hsts = []
	for biH=1:dims(colpt_type(flow_t, :b))
		b = pb #- a^2/8*HexC.Δ(cpb, cpb, pb) #
		hst = u⃗∇ᵀ(colptidx(0,0,biH,Val(colpt_type(flow_t, :b))), cpu⃗, cpb, pu⃗, b)
		hst = substitute(hst, GridOperatorAnalysis.sqrt3_subs)
		push!(hsts, hst)
	end
	hsts
end;

# ╔═╡ 03ccdcf0-5dcb-4257-ace8-8f9ce22dee3e
lhsts = [
	expand(taylor_coeff(expand(hst), ϵ, 1)) for hst in hsts
];

# ╔═╡ 16bb4e26-9d82-4718-881f-66450ce90bf1
fhsts = [
	fourier_transform_expression(biH, colpt_type(flow_t, :b), lhsts[biH]; fflow, dflow, ϕ) for biH=1:dims(colpt_type(flow_t, :b))
];

# ╔═╡ 391de831-feea-436e-9be5-15686d9c9155
md"""
### Rewriter
"""

# ╔═╡ ed74cf26-7ae6-43f6-9303-09c46b0fae0a
rtrig = let
	function p(x)
		 !isequal(x, θU)
	end
	rcos = let
    	x = Symbolics.variable(:x)
    	Symbolics.@rule cos(~x::p) => substitute(taylor(cos(x),x,0:4), Dict([x=>~x]))
	end
	rsin = let
    	x = Symbolics.variable(:x)
    	Symbolics.@rule sin(~x::p)=>substitute(taylor(sin(x),x,0,0:4), Dict([x=>~x]))
	end
	SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(SymbolicUtils.RestartedChain([rcos, rsin])))
end

# ╔═╡ 192a4172-a2bb-416b-83df-194a090b093a
sfhmts = let
	colpt_t = colpt_type(flow_t, :u⃗)
	d       = dims(colpt_t)
	sfhmts = if colpt_t == :edge
		zeros(Complex{Symbolics.Num}, d, d)
	else
		zeros(Complex{Symbolics.Num}, d, 2, d, 2)
	end
	for iH=1:d
		fhmt = let
			fhmt = expand.(fhmts[iH])
			fhmt = if doapprox # rewrite trig-functions by Taylor polynomial
				simplify.(fhmt; rewriter=rtrig) .|> expand
			else
				substitute.(fhmt, Ref(Dict(l => h/a * 2/sqrt3 * l)))
			end
			fhmt = substitute.(fhmt, Ref(Dict(le=>a, √(f₀^2*N²/Ri) => M²))) .|> expand
			fhmt = substitute.(fhmt, Ref(sqrt3_subs)) .|> expand
		end
		if colpt_t == :edge
			for jH = 1:d
				u = fflow.u⃗[jH]
				sfhmt = taylor_coeff(fhmt, u, 1)
				sfhmts[iH, jH] = simplify(sfhmt; expand=true)
			end
		else
			for iTH=1:2
				for jH=1:d
					for jTH=1:2
						u = fflow.u⃗[jTH, jH]
						sfhmt = taylor_coeff(fhmt[iTH], u, 1)
						sfhmts[iH, iTH, jH, jTH] = simplify(sfhmt; expand=true)
					end
				end
			end
		end
	end
	sfhmts
end;

# ╔═╡ 4c11881b-8b92-4f9e-807d-0d2685893c82
begin
	@variables K
	Fn = let
		F = if colpt_type(flow_t, :u⃗) == :edge
			sfhmts ./ Ū
		elseif grid_t == :TriB && hmt_scheme == :avi
			reshape(sfhmts ./ Ū, 4, 4)
		else	
			sfhmts[:,1,:,1] ./ Ū
		end
		F = substitute.(F, Ref(Dict(k => K * cos(θU), l => K * sin(θU))))
		if doapprox
			F .= substitute.(F, Ref(Dict(a=>2/sqrt3 * h)))
			F .= simplify.(F; expand=true)
			F .= substitute.(F, Ref(sqrt3_subs))
		else
			F = expand(F)
			F = substitute.(F, Ref(Dict(sqrt3 => 2 * h/a)))
		end
	end
	Fn = simplify.(Fn ./ K; expand=true)
end

# ╔═╡ 78fa487a-5bdf-4aa3-b9c5-25d856697fae
begin
	sfbs, sfus = let
		colpt_t_b = colpt_type(flow_t, :b)
		colpt_t_u = colpt_type(flow_t, :u⃗)
		db = dims(colpt_t_b)
		du = dims(colpt_t_u)
		if colpt_type(flow_t, :u⃗) == :edge
			(zeros(Complex{Symbolics.Num},db,db), zeros(Complex{Symbolics.Num},db,du))
		else
			(zeros(Complex{Symbolics.Num},db,db), zeros(Complex{Symbolics.Num},db,du,2))
		end
	end
	for iH=1:dims(colpt_type(flow_t, :b))
		fhst = let
			fhst = expand.(fhsts[iH])
			if doapproxs
				fhst = expand.(simplify.(fhst; rewriter=rtrig))
			else
				fhst = substitute.(fhst, Ref(Dict(l => h/a * 2/GridOperatorAnalysis.sqrt3 * l)))
			end
			fhst = substitute.(fhst, Ref(Dict(le=>a, √(f₀^2*N²/Ri) => M²))) .|> expand
			fhst = substitute.(fhst, Ref(sqrt3_subs)) .|> expand
		end
		for jH=1:dims(colpt_type(flow_t, :b))
			b = fflow.b[jH]
			sfb = taylor_coeff(fhst, b, 1)
			sfbs[iH, jH] = simplify(sfb; expand=true)
		end
		for jH=1:dims(colpt_type(flow_t, :u⃗))
			if colpt_type(flow_t, :u⃗) == :edge
				u = fflow.u⃗[jH]
				sfu = taylor_coeff(fhst, u, 1)
				sfus[iH, jH] = simplify(sfu; expand=true)
			else
				for jTH=1:2
					u = fflow.u⃗[jTH, jH]
					sfu = taylor_coeff(fhst, u, 1)
					sfus[iH, jH, jTH] = simplify(sfu; expand=true)
				end
			end
		end
	end
end

# ╔═╡ 70137bf3-0f04-4863-ad91-2087ae3f1d3f
begin
	@variables K̃
	Fb = let
		F = sfbs./ Ū
		F = substitute.(F, Ref(Dict(k => K * cos(θU), l => K * sin(θU))))
		if doapproxs
			F .= substitute.(F, Ref(Dict(a=>2/sqrt3 * h)))
			F .= substitute.(F, Ref(Dict(h=> K̃ / K, sin(θU)^2=>1-cos(θU)^2, sin(θU)^4=>(1-cos(θU)^2)^2, sin(θU)^6=>(1-cos(θU)^2)^3))) #
			F .= simplify.(F; expand=true)
			F .= substitute.(F, Ref(sqrt3_subs))
		else
			F .= expand.(F)
			F .= substitute.(F, Ref(Dict(sqrt3 => 2 * h/a)))
		end
	end
	Fb .= simplify.(Fb ./ K; expand=true)
end

# ╔═╡ 7509a8b9-468b-4b2b-8a8c-175b3f045c4b
begin
	T = let
		F = sfus ./ M²
		F = substitute.(F, Ref(Dict(k => K * cos(θU), l => K * sin(θU))))
		if doapproxs
			F .= substitute.(F, Ref(Dict(a=>2/sqrt3 * h))) .|> expand
			F .= substitute.(F, Ref(Dict(
				h=> K̃ / K, 
				sin(θU)^2=>1-cos(θU)^2, 
				sin(θU)^3=>sin(θU)*(1-cos(θU)^2),
				sin(θU)^4=>(1-cos(θU)^2)^2,
				sin(θU)^5=>sin(θU)*(1-cos(θU)^2)^2, 
				sin(θU)^6=>(1-cos(θU)^2)^3
			)))
			F .= simplify.(F; expand=true)
			F .= substitute.(F, Ref(sqrt3_subs))
		else
			F .= expand.(F)
			F .= substitute.(F, Ref(Dict(sqrt3=> 2 * h/a)))
		end
	end
	T = if colpt_type(Val(grid_t), :u⃗) ≠ :edge
		simplify.(T ./ reshape([-sin(θU); cos(θU)], 1, 1, 2); expand=true)
		#(simplify.(T[:,:,1] ./ -sin(θU); expand=true), simplify.(T[:,:,2] ./ cos(θU); expand=true))
	else
		T
	end
	if colpt_type(Val(grid_t), :u⃗) ≠ :edge
		(T[:,:,1], T[:,:,2])
	else
		T
	end
end

# ╔═╡ 27dbf7f2-14db-40ac-9e59-d2cfb3864a8e
begin
	sfhds = let
		colpt_t_b = colpt_type(flow_t, :b)
		db = dims(colpt_t_b)
		zeros(Complex{Symbolics.Num},db,db)
	end
	for iH=1:dims(colpt_type(flow_t, :b))
		fhd = let
			fhd = expand.(fhds[iH])
			if doapproxs
				fhd = expand.(simplify.(fhd; rewriter=rtrig))
			else
				fhd = substitute.(fhd, Ref(Dict(l => h/a * 2/GridOperatorAnalysis.sqrt3 * l)))
			end
			fhd = substitute.(fhd, Ref(Dict(le=>a, √(f₀^2*N²/Ri) => M² ))) .|> expand
			fhd = substitute.(fhd, Ref(sqrt3_subs)) .|> expand
		end
		for jH=1:dims(colpt_type(flow_t, :b))
			b = fflow.b[jH]
			sfhd = taylor_coeff(fhd, b, 1)
			sfhds[iH, jH] = simplify(sfhd; expand=true)
		end
	end
end

# ╔═╡ 04d312fc-e353-4660-82ad-3dd048d79551
begin
	Fhd = let
		F = sfhds
		F = substitute.(F, Ref(Dict(k => K * cos(θU), l => K * sin(θU))))
		if doapproxs
			F .= substitute.(F, Ref(Dict(a=>2/sqrt3 * h)))
			F .= simplify.(F; expand=true)
			F .= substitute.(F, Ref(sqrt3_subs))
			F .= substitute.(F, Ref(Dict(h=> K̃ / K, cos(θU)^2=>1-sin(θU)^2, cos(θU)^4=>(1-sin(θU)^2)^2))) #
		else
			F = expand.(F)
			F = substitute.(F, Ref(Dict(sqrt3 => 2 * h/a)))
		end
	end
	@show Fhd = simplify.(Fhd ./ -K^2; expand=true)
end

# ╔═╡ 7b864ab1-b61c-481c-849a-9824c9abc984
rpyt = let
	r = Symbolics.@acrule sin(~x)^2 + cos(~x)^2 => one(~x)
	SymbolicUtils.Prewalk(SymbolicUtils.PassThrough(r))
end

# ╔═╡ d6ebd665-16c7-4850-8302-80356de08639
md"""
### Misc
"""

# ╔═╡ c9a610d4-ecd7-4798-8b5d-88c96ad6afae
rewriter = let
	r = Symbolics.@acrule ~α * ~(~x) + ~β * ~(~x) => *(~α + ~β, ~(~x)...)
	rewriter = SymbolicUtils.Prewalk(SymbolicUtils.PassThrough(r))
end

# ╔═╡ f0e541cf-b2f5-40f5-8379-b1f8ef1a4fe8
let
	@variables a[1:300]
	r = Symbolics.@rule ~α * ~(~x) + ~β * ~(~x) => *(~α + ~β, ~(~x)...)
	rewriter = SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(r))
	t = +([a[i] * a[i+1] * a[i+2] * a[i+3] * a[i+4] for i=1:5:300]...)
	simplify(t; rewriter)
end

# ╔═╡ 4c17d92e-cdf7-442f-9087-102ff8109a8a
function Base.cbrt(z::Complex{Symbolics.Num})
		@variables x
		vars = [Symbolics.get_variables(real(z)); Symbolics.get_variables(imag(z))]
		Dx = Symbolics.Differential(x)
		z̄ = Symbolics.value(substitute(z, Dict(vars .=> 0)))
		cs = [substitute(Symbolics.expand_derivatives(1/factorial(k)*(Dx^k)(cbrt(x))), Dict(x=>z̄)) for k=0:2]
		cbrtz = Symbolics.series(cs, z-z̄)
		Symbolics.simplify(cbrtz; expand=true)
	end

# ╔═╡ e667e7fd-ab40-49c9-b8ca-a4aff155e833
Base.cbrt(z::Complex{<:Real}) = z^(1/3)

# ╔═╡ dc22802f-307d-4bf1-a10d-3e1b0180783e
function cubic(d, c, b, a)
	d, c, b, a = complex.([d, c, b, a])
	Δ₀ = b^2 - 3 * a * c
	Δ₁ = 2*b^3 - 9 * a * b * c + 27 * a^2 * d
	C = cbrt((Δ₁ + sqrt(Δ₁^2 - 4 * Δ₀^3)) / 2)
	vars = [Symbolics.get_variables(real(C)); Symbolics.get_variables(imag(C))]
	substitute(C, Dict(vars .=> 0))
	Cs = [C * exp(im*2*π*k/3) for k=0:2]
	[-1/(3*a) * (b + C + Δ₀/C) for C in Cs]
end

# ╔═╡ e80b16a7-315f-4b3f-b5ce-4cd73b6a3fdb
compute_symbolic_eigenvals = let
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
	function compute_symbolic_eigenvals(A)
	    @variables λ μ::Complex # eigenvalue  
	    # find eigenvalues first
	    p = Symbolics.expand(Symbolics.det(complex.(λ*I- A); laplace=true)) # polynomial to solve
		if Symbolics.degree(real(p), λ) ≤ 1 && Symbolics.degree(imag(p), λ) ≤ 1
			cs = Symbolics.taylor_coeff(p, λ, 0:1)
			p = Symbolics.series(cs, μ)
			GridOperatorAnalysis.symbolic_linear_solve(p ~ 0, μ)
		elseif Symbolics.degree(real(p), λ) ≤ 2 && Symbolics.degree(imag(p), λ) ≤ 2
	    	quadratic(Symbolics.taylor_coeff(p, λ, 0:2)...) # solve polynomial
		else
			cubic(Symbolics.taylor_coeff(p, λ, 0:3)...)
		end
	end
end

# ╔═╡ 16535b97-0c7e-4c14-9e71-604b0e0c9a3d
if doapprox
	svs = compute_symbolic_eigenvals(Symbolics.wrap.(Fn))
	svs = Symbolics.simplify.(svs; expand=true)
else
	"Computation of symbolic eigenvalues not implemented in this case."
end

# ╔═╡ f1feff04-e468-4909-af4c-9bc5a032f407
K̃s, ωs = let
	@variables K̃
	K̃s = range(1e-20, π, 100) #floatmin(Float64)
	ωs = []
	if doapprox
		for sv in svs
			sv    = substitute(sv, Dict(
				θU    => 0, 
				h     => K̃ / K, 
				a     => K̃ / K * 2 / √3, 
				sqrt3 => √3))
			sv    = simplify(sv; expand=true)
			vfunc = Symbolics.build_function(sv, K̃; expression=Val(false))
			push!(ωs, vfunc.(K̃s))
		end
	else
		A = substitute.(Fn, Ref(Dict(
			θU => 0, 
			h => K̃ / K, 
			a => K̃ / K * 2 / √3, 
			sqrt3 => √3
		)))
		A = simplify.(A; expand=true)
		Afunc = Symbolics.build_function(A, K̃; expression=Val(false))[1]
		for k̃ in K̃s
			(; values) = eigen(Complex{Float64}.(Afunc(k̃)))
			push!(ωs, values)
		end
	end
	(K̃s, ωs)
end;

# ╔═╡ 54131367-66b2-4bc7-8820-0b6a5fc19ae3
let
	f = Figure(; fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / h⁻¹",
			 ylabel = L"ω / i\bar{U}K",
			 limits = (0, π, -1, 1.1), 
			 aspect = 1,
			 )
	if doapprox
		for ω in ωs
			lines!(ax, K̃s, imag.(ω), linewidth=3)
		end
	else
		for (k̃, ω) in zip(K̃s, ωs)
			for _ω in ω
				scatter!(ax, k̃, imag(_ω); color=:blue)
			end
		end
	end
	#axislegend()
	f
end

# ╔═╡ 6dab1873-d738-4e41-9d16-0f83fb168bd4
if doapproxs
	svbs = compute_symbolic_eigenvals(Symbolics.wrap.(Fb))
	svbs = Symbolics.simplify.(svbs; expand=true)
else
	"Computation of symbolic eigenvalues not implemented in this case."
end

# ╔═╡ 6a120e6b-ced8-414a-a56b-3b49ea21c8fa
let
	K̃s, ωs = let
		K̃s = range(1e-20, π, 100) #floatmin(Float64)
		ωs = []
		if doapproxs
			for sv in svbs
				sv = substitute(sv, Dict(
					θU    => 0, 
					h     => K̃ / K, 
					a     => K̃ / K * 2 / √3, 
					sqrt3 => √3))
				sv    = simplify(sv; expand=true)
				vfunc = Symbolics.build_function(sv, K̃; expression=Val(false))
				push!(ωs, vfunc.(K̃s))
			end
		else
			A = substitute.(Fb, Ref(Dict(
				θU => 0, 
				h => K̃ / K, 
				a => K̃ / K * 2 / √3, 
				sqrt3 => √3
			)))
			A = simplify.(A; expand=true)
			Afunc = Symbolics.build_function(A, K̃; expression=Val(false))[1]
			for k̃ in K̃s
				(; values) = eigen(Complex{Float64}.(Afunc(k̃)))
				push!(ωs, values)
			end
		end
		(K̃s, ωs)
	end
	f = Figure(; fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / h⁻¹",
			 ylabel = L"ω / i\bar{U}K",
			 limits = (0, π, -1, 1.1), 
			 aspect = 1,
			 )
	if doapproxs
		for ω in ωs
			lines!(ax, K̃s, imag.(ω), linewidth=3)
		end
	else
		for (k̃, ω) in zip(K̃s, ωs)
			for _ω in ω
				scatter!(ax, k̃, imag(_ω); color=:blue)
			end
		end
	end
	#axislegend()
	f
end

# ╔═╡ 5d83c523-626f-4b06-85ee-8b71423d398a
# ╠═╡ disabled = true
#=╠═╡
if doapproxs
	svus = compute_symbolic_eigenvals(Symbolics.wrap.(T[1]))
	svus = Symbolics.simplify.(svus; expand=true)
else
	"Computation of symbolic eigenvalues not implemented in this case."
end
  ╠═╡ =#

# ╔═╡ 9b3f3f72-ffb8-477d-b1a5-4a14d64d36e8
#=╠═╡
let
	K̃s, ωs = let
		K̃s = range(1e-20, π, 100) #floatmin(Float64)
		ωs = []
		if doapproxs
			for svu in svus
				svu    = substitute(svu, Dict(
					θU    => 0, 
					h     => K̃ / K, 
					a     => K̃ / K * 2 / √3, 
					sqrt3 => √3))
				svu    = simplify(svu; expand=true)
				vfunc = Symbolics.build_function(svu, K̃; expression=Val(false))
				push!(ωs, vfunc.(K̃s))
			end
		else
			A = substitute.(T[1], Ref(Dict(
				θU => 1e-30, 
				h => K̃ / K, 
				a => K̃ / K * 2 / √3, 
				sqrt3 => √3
			)))
			A = simplify.(A; expand=true)
			Afunc = Symbolics.build_function(A, K̃; expression=Val(false))[1]
			for k̃ in K̃s
				(; values) = eigen(Complex{Float64}.(Afunc(k̃)))
				push!(ωs, values)
			end
		end
		(K̃s, ωs)
	end
	f = Figure(; fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / h⁻¹",
			 ylabel = L"ω / -\mathrm{M²}\sin(θU)",
			 limits = (0, π, 0, 2), 
			 aspect = 1,
			 )
	if doapproxs
		for ω in ωs
			lines!(ax, K̃s, real.(ω), linewidth=3)
		end
	else
		for (k̃, ω) in zip(K̃s, ωs)
			for _ω in ω
				scatter!(ax, k̃, real(_ω); color=:blue)
			end
		end
	end
	#axislegend()
	f
end
  ╠═╡ =#

# ╔═╡ 30620b89-ced4-4353-bf3f-6d9c58794aa4
# ╠═╡ disabled = true
#=╠═╡
if doapproxs
	svvs = compute_symbolic_eigenvals(Symbolics.wrap.(T[2]))
	svvs = Symbolics.simplify.(svvs; expand=true)
else
	"Computation of symbolic eigenvalues not implemented in this case."
end
  ╠═╡ =#

# ╔═╡ 5361ebe0-a415-416c-b57f-051b0e6382f6
#=╠═╡
let
	K̃s, ωs = let
		K̃s = range(1e-20, π, 100) #floatmin(Float64)
		ωs = []
		if doapproxs
			for svv in svvs
				svv = substitute(svv, Dict(
					θU    => 0, 
					h     => K̃ / K, 
					a     => K̃ / K * 2 / √3, 
					sqrt3 => √3))
				svv   = simplify(svv; expand=true)
				vfunc = Symbolics.build_function(svv, K̃; expression=Val(false))
				push!(ωs, vfunc.(K̃s))
			end
		else
			A = substitute.(T[2], Ref(Dict(
				θU => 0, 
				h => K̃ / K, 
				a => K̃ / K * 2 / √3, 
				sqrt3 => √3
			)))
			A = simplify.(A; expand=true)
			Afunc = Symbolics.build_function(A, K̃; expression=Val(false))[1]
			for k̃ in K̃s
				(; values) = eigen(Complex{Float64}.(Afunc(k̃)))
				push!(ωs, values)
			end
		end
		(K̃s, ωs)
	end
	f = Figure(; fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / h⁻¹",
			 ylabel = L"ω / \mathrm{M²}\cos(θU)",
			 limits = (0, π, 0, 2), 
			 aspect = 1,
			 )
	if doapproxs
		for ω in ωs
			lines!(ax, K̃s, real.(ω), linewidth=3)
		end
	else
		for (k̃, ω) in zip(K̃s, ωs)
			for _ω in ω
				scatter!(ax, k̃, real(_ω); color=:blue)
			end
		end
	end
	#axislegend()
	f
end
  ╠═╡ =#

# ╔═╡ dcbd2629-eb5d-42b5-9234-d132106f5747
if doapproxs
	sevs = compute_symbolic_eigenvals(Symbolics.wrap.(Fhd*K̃^2))
	sevs = Symbolics.simplify.(sevs; expand=true)
else
	"Computation of symbolic eigenvalues not implemented in this case."
end

# ╔═╡ Cell order:
# ╟─b9cca8c1-0d0f-48f9-b21b-8d6df2cb77aa
# ╠═f59a5438-8d5e-11f0-13fa-d9c703e5f87f
# ╟─59056484-9c6b-48ad-b741-2bc294d5cc6f
# ╟─0943c748-9fee-401e-851e-de2327d39706
# ╠═d2b72544-78fb-4a0b-b169-aa50ebaeb31d
# ╟─9abfbc35-bd37-4554-949c-27cd2bdfa1a7
# ╠═cfb78fd5-4435-4c0d-8cf2-1b5fffc2f91e
# ╟─428cc192-03d0-497c-aee7-f7015c4c09d5
# ╠═fb5e7874-c067-4769-aefa-260fd3eca00b
# ╟─a697761a-3f01-46ac-ae31-4665bc7b6e47
# ╠═2663542e-c719-42f0-bf00-790e8aa211fe
# ╠═1b94a213-9311-4dad-bc3e-0ac2dcb71e14
# ╟─430cca4e-08b6-4d57-a4dc-422dcf8cc92e
# ╠═8ae523f0-1eeb-4d62-b089-502c52dd1277
# ╟─77eb6488-c26a-416a-839b-1e617f0cdd50
# ╠═780abc6b-d9ab-40d1-9a01-d12b1d3fc3ae
# ╟─3b152922-e343-42a1-b863-38814a80b0d7
# ╟─33e1ab03-9e1a-4014-84ed-1ad393f0445e
# ╟─f9adc84b-12c3-4a1c-8400-6ede0108c9e9
# ╟─39031d65-863d-48cf-aeaf-520417eaf071
# ╠═a64854b8-b161-48d2-87da-ebed3e08671c
# ╠═50fc9497-99e4-44d7-b25f-ffccd14b145a
# ╟─286c65de-46f4-462d-903c-ad850bb1a3b1
# ╠═f71432b9-cb86-4f5e-a6ba-4ad443140292
# ╠═26b13d30-d8cf-4ba4-9f84-a591bb830630
# ╟─583e619a-76da-4f63-8377-b5eac40c96af
# ╟─0f1cbe99-8b11-426d-8d28-77b6dffef277
# ╠═192a4172-a2bb-416b-83df-194a090b093a
# ╠═4c11881b-8b92-4f9e-807d-0d2685893c82
# ╟─46afe77d-f904-40e9-92bb-cab83aeea51e
# ╠═16535b97-0c7e-4c14-9e71-604b0e0c9a3d
# ╠═f1feff04-e468-4909-af4c-9bc5a032f407
# ╟─54131367-66b2-4bc7-8820-0b6a5fc19ae3
# ╟─990c73d7-0d2a-457c-81c7-88a1f80cb44a
# ╟─27ec48bb-d7bb-4b70-9932-08f0e0926504
# ╟─e9dd8811-c91e-45e4-8ced-afd7391f248a
# ╠═f80621ff-53de-4b21-b314-b8d4bbdbcc61
# ╟─363d736c-c30b-4878-8f32-c066df479600
# ╠═03ccdcf0-5dcb-4257-ace8-8f9ce22dee3e
# ╟─386315ca-ef94-44f5-9df8-eb320be72b34
# ╠═16bb4e26-9d82-4718-881f-66450ce90bf1
# ╟─698f9004-f56b-41c0-b6a6-b55c505fb1a8
# ╟─99c92534-3092-4995-8044-d8001aad8f5c
# ╠═78fa487a-5bdf-4aa3-b9c5-25d856697fae
# ╟─ed799518-6bf1-4180-8206-0a3b6e51aed0
# ╠═70137bf3-0f04-4863-ad91-2087ae3f1d3f
# ╟─05fb2f19-1b43-4e7a-84c5-1ea8a71b9284
# ╠═7509a8b9-468b-4b2b-8a8c-175b3f045c4b
# ╟─26fd8b72-c736-43db-9d32-f76cbe5dd901
# ╠═6dab1873-d738-4e41-9d16-0f83fb168bd4
# ╟─5d83c523-626f-4b06-85ee-8b71423d398a
# ╠═30620b89-ced4-4353-bf3f-6d9c58794aa4
# ╠═6a120e6b-ced8-414a-a56b-3b49ea21c8fa
# ╟─9b3f3f72-ffb8-477d-b1a5-4a14d64d36e8
# ╟─5361ebe0-a415-416c-b57f-051b0e6382f6
# ╟─bc991b63-8e0e-4e37-9369-15f2e11e195e
# ╠═7f20375f-6491-47ff-bcc0-fb77baa8bedf
# ╠═1af2939c-6795-4472-ae44-8c9466a159ab
# ╠═4b2d9f6d-4438-436a-91cb-47fb66270841
# ╠═27dbf7f2-14db-40ac-9e59-d2cfb3864a8e
# ╠═04d312fc-e353-4660-82ad-3dd048d79551
# ╠═dcbd2629-eb5d-42b5-9234-d132106f5747
# ╟─7d55ea61-e0f1-4bd9-a36c-22857ad145e7
# ╟─8eea87eb-048e-4436-9d20-3bcc9e76c7b7
# ╟─f0f8c55d-82fa-4cf4-9678-555c1222e615
# ╟─9b8a0e33-2422-48e6-93d4-0313c557abe3
# ╟─e18fbb3e-02cb-4add-ab00-9e0cc7ae935f
# ╟─e322d734-1c1a-4ce1-a147-9f8420c64203
# ╟─3782fb16-26b7-4edd-8591-ef119b1cabef
# ╟─d424025e-3efd-4c4b-8d97-7369ff3618da
# ╟─423c5bfc-d6e2-46cc-a29f-05e3d9d3c84b
# ╟─391de831-feea-436e-9be5-15686d9c9155
# ╠═ed74cf26-7ae6-43f6-9303-09c46b0fae0a
# ╠═7b864ab1-b61c-481c-849a-9824c9abc984
# ╟─d6ebd665-16c7-4850-8302-80356de08639
# ╠═c9a610d4-ecd7-4798-8b5d-88c96ad6afae
# ╠═f0e541cf-b2f5-40f5-8379-b1f8ef1a4fe8
# ╠═e80b16a7-315f-4b3f-b5ce-4cd73b6a3fdb
# ╠═4c17d92e-cdf7-442f-9087-102ff8109a8a
# ╠═dc22802f-307d-4bf1-a10d-3e1b0180783e
# ╠═e667e7fd-ab40-49c9-b8ca-a4aff155e833
