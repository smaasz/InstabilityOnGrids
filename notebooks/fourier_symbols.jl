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

# ╔═╡ d60c197f-28db-4834-b78f-aab7f23536f2
simplifyexpand(x) = simplify(x; expand=true)

# ╔═╡ 59056484-9c6b-48ad-b741-2bc294d5cc6f
TableOfContents(; depth=4)

# ╔═╡ 0943c748-9fee-401e-851e-de2327d39706
md"""
## Grid
"""

# ╔═╡ 9abfbc35-bd37-4554-949c-27cd2bdfa1a7
md"""
__Grid__: $(@bind grid_t Select([:TriA, :TriB, :TriC, :HexC]; default=:TriA))
"""

# ╔═╡ cfb78fd5-4435-4c0d-8cf2-1b5fffc2f91e
@variables f₀ g N² Ri le ϵ k l h a M² β θU 𝕂ᵘ 𝕂ᵇ

# ╔═╡ d2b72544-78fb-4a0b-b169-aa50ebaeb31d
colpts = (; vertex=v, cell=c, edge=e);

# ╔═╡ d32f2049-e2d4-4073-865f-585cc80d4d07
phase_subs = GridOperatorAnalysis.phasesubs(a, h, k, l);

# ╔═╡ 4a5bd3c3-bcfe-498e-8fa7-328cd12cea00
invphase_subs = GridOperatorAnalysis.invphasesubs(a, h, k, l);

# ╔═╡ 428cc192-03d0-497c-aee7-f7015c4c09d5
md"""
#### Background flow

The baroclinically unstable flow with a sloping buoyancy profile in across-front direction (first studied by Eady) is projected onto the respective discrete spaces.
The flow is parametrized by:
- the Richardson number $\mathrm{Ri}$, 
- the Brunt-Väisälä frequency $\mathrm{N^2}$, 
- the Coriolis frequency $f_0$ on the f-plane, 
- the flow direction $θ_U$, 
- and the flow shift $\beta$ (in $H\mathrm{M}^2/f_0$).

From the thermal wind balance, the derivative of the buoyancy in meridional direction may be computed

$$\mathrm{M}^2 = \sqrt{\frac{f₀^2* \mathrm{N}²}{\mathrm{Ri}}}.$$

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

# ╔═╡ da53d71b-e5b3-4ce0-8ce8-1a0449c4fd32
md"""
__Dissipation scheme__: $(@bind dissip_scheme Select([:biharmonic => "biharmonic", :harmonic => "harmonic"]))
"""

# ╔═╡ 77eb6488-c26a-416a-839b-1e617f0cdd50
md"""
#### Flow in Wavenumber Space
"""

# ╔═╡ 7741fc90-ad1f-4dab-baa5-bbf373592e13
@variables K̃

# ╔═╡ 0f1cbe99-8b11-426d-8d28-77b6dffef277
md"""
__Small wavenumber approximation__: $(@bind doapprox CheckBox(default=true))
"""

# ╔═╡ 77fe42fd-ddbb-495a-979e-0a852600a2a6
md"""
__Phase substitutions__: $(@bind dophasesubs CheckBox(default=false))
"""

# ╔═╡ e287982b-5fea-4d73-a0ce-d2627635ddb2
if grid_t == :TriC
L"""
\begin{bmatrix}
U\mathsf{\underline{G_x}} + V\mathsf{\underline{G_y}} + \left( U_z\mathsf{\underline{A^{(x)}}} + V_z\mathsf{\underline{A^{(y)}}}\right)\mathsf{\underline{W}} + f_0 \mathsf{\underline{M}}& \mathsf{\underline{g}}~\mathsf{\underline{P}} & g  \mathsf{g\otimes 1_V}\\
N^2 \mathsf{\underline{W}} + B_x \mathsf{\underline{Av^{(x)}}} + B_y \mathsf{\underline{Av^{(y)}}} & U \mathsf{\underline{\Gamma_x}} + V\mathsf{\underline{\Gamma_y}} & \underline{0}\\
\Delta_z \mathsf{D^{(x)} \otimes 1_V^T} & \mathsf{0 \otimes 1_V^T} & \mathsf{0}
\end{bmatrix}\begin{bmatrix} \underline{\boldsymbol{\hat{u}}} \\ \underline{\hat{b}} \\ \hat{\eta} \end{bmatrix}
"""
else
L"""
\begin{bmatrix}
U\mathsf{\underline{G_x^{(xx)}}} + V\mathsf{\underline{G_y^{(xx)}}} + U_z\mathsf{\underline{A^{(xx)}}} \mathsf{\underline{W^{(x)}}} & U\mathsf{\underline{G_x^{(xy)}}} + V\mathsf{\underline{G_y^{(xy)}}} + U_z\mathsf{\underline{A^{(xy)}}} \mathsf{\underline{W^{(y)}}} - f_0 \mathsf{\underline{I}} & \mathsf{\underline{g_x}}~\mathsf{\underline{P}} & g  \mathsf{g_x\otimes 1_V}\\
U\mathsf{\underline{G_x^{(yx)}}} + V\mathsf{\underline{G_y^{(yx)}}} + V_z\mathsf{\underline{A^{(yx)}}} \mathsf{\underline{W^{(x)}}} + f_0 \mathsf{\underline{I}} & U\mathsf{\underline{G_x^{(yy)}}} + V\mathsf{\underline{G_y^{(yy)}}} + V_z\mathsf{\underline{A^{(yy)}}} \mathsf{\underline{W^{(y)}}} & \mathsf{\underline{g_y}}~\mathsf{\underline{P}} & g\mathsf{g_y\otimes 1_V}\\
N^2 \mathsf{\underline{W^{(x)}}} + B_x \mathsf{\underline{Av^{(xx)}}} + B_y \mathsf{\underline{Av^{(yx)}}} & N^2 \mathsf{\underline{W^{(y)}}} + B_x \mathsf{\underline{Av^{(xy)}}} + B_y \mathsf{\underline{Av^{(yy)}}} & U \mathsf{\underline{\Gamma_x}} + V\mathsf{\underline{\Gamma_y}} & \underline{0}\\
\Delta_z \mathsf{D^{(x)} \otimes 1_V^T} & \Delta_z\mathsf{D^{(y)} \otimes 1_V^T} & \mathsf{0 \otimes 1_V^T} & \mathsf{0}
\end{bmatrix}\begin{bmatrix} \underline{\boldsymbol{\hat{u}}} \\ \underline{\boldsymbol{\hat{v}}} \\ \underline{\hat{b}} \\ \hat{\eta} \end{bmatrix}
"""
end

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

# ╔═╡ 729a626e-376c-4ad9-96e8-f2e641f4d2a4
md"""
#### ``\mathsf{G_x}``
"""

# ╔═╡ 30725e26-6949-4605-a275-3f72eca04ece
md"#### ``\mathsf{G_y}``"

# ╔═╡ 0e2ee007-f695-408b-9a8f-3619b72c137f
# ╠═╡ disabled = true
#=╠═╡
push!(Gs, (_gx, _gy))
  ╠═╡ =#

# ╔═╡ 33dcefea-8d67-4543-a241-3192ad3252d2
md"#### ``\mathsf{M}``"

# ╔═╡ d40d6c8c-f814-4c98-989c-a9b3d7913c0c
md"#### ``\mathsf{A^{(x)}}``"

# ╔═╡ 801f824d-9fc4-431c-bbbe-a72ed0d3b3ac
md"#### ``\mathsf{A^{(y)}}``"

# ╔═╡ 7d3ce7a0-eb0b-4801-8dec-84ba360dd477
md"#### ``\mathsf{g}``"

# ╔═╡ 4a18efe0-0d3e-481c-94a4-d1b3e4d56df3
md"#### ``\mathsf{D^{u}}``"

# ╔═╡ 46afe77d-f904-40e9-92bb-cab83aeea51e
md"""
#### Eigenvalues of Fourier Symbols
"""

# ╔═╡ 5708950e-80fb-4296-adfa-e2db01284d06
df = []

# ╔═╡ 4ffe2716-9a66-4ca5-bf4b-367e85399a5b
# ╠═╡ disabled = true
#=╠═╡
push!(df, hmt_scheme => (K̃s, ωs))
  ╠═╡ =#

# ╔═╡ 0890ade4-59b2-45cf-a71d-35f918c71883
md"""
θU: $(@bind sθU PlutoUI.Slider([1e-20, π/12, π/6]; show_value=true))
"""

# ╔═╡ 57b612b5-afb6-44db-ad3f-6c9882afae14
let
	f = Figure(; fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = "wavenumber / h⁻¹",
			 ylabel = L"ω / i\bar{U}K",
			 limits = (0, π, -2.1, 1.1), 
			 aspect = 1,
			 )
	colors = [:red, :green, :black, :blue]
	linestyles = [:solid, :solid, :dash, :solid]
	for (ischeme, (scheme, (K̃s, ωs))) in enumerate(df)
		linestyle = linestyles[ischeme]
		color = colors[ischeme]
		ωs = stack(ωs)
		ωs = sort(ωs, dims=1, by=imag)
		for i = 1:size(ωs, 1)
			lines!(ax, K̃s, imag.(ωs[i,:]); color, linestyle, label=String(scheme), linewidth=3)
		end
	end
	#axislegend(; merge=true)
	Legend(f[1,2], ax; merge=true)
	f
end

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

# ╔═╡ 77960a21-8754-427c-875b-e6578e0ed62f
md"#### ``\mathsf{\Gamma_x}``"

# ╔═╡ da3f6839-042f-4617-abbf-fb5a18fe772a
md"#### ``\mathsf{\Gamma_y}``"

# ╔═╡ 05fb2f19-1b43-4e7a-84c5-1ea8a71b9284
md"""
#### Advection of the background buoyancy by the velocity pertubation
"""

# ╔═╡ f4d63a77-6141-41df-bcdf-f7547d4e3cf5
md"#### ``\mathsf{Av^{(x)}}``"

# ╔═╡ fd7bfb37-1bfb-41ff-8dd1-0c566156693d
md"#### ``\mathsf{Av^{(y)}}``"

# ╔═╡ c6321e6c-e4ce-45fa-b884-69f14535db2d
md"#### ``\mathsf{D^{b}}``"

# ╔═╡ 26fd8b72-c736-43db-9d32-f76cbe5dd901
md"""
#### Eigenvalues of Fourier Symbols
"""

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

# ╔═╡ 932e3e09-c15f-4b0e-944a-1eecbcf378d7
dfs = []

# ╔═╡ 23a4913d-3686-42a0-a2ce-c34b298fa59f
md"""
θU: $(@bind tθU PlutoUI.Slider([0.0, π/12, π/6]; show_value=true))
"""

# ╔═╡ 0e90c7a9-fc7e-4890-b6e0-32ca12e8ae7d
let
	f = Figure(; fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = L"\text{wavenumber}~K / h^{-1}",
			 ylabel = L"ω / i\bar{U}K",
			 limits = (0, π, -0.1, 1.1), 
			 aspect = 1,
			 )
	colors = [:red, :green, :black, :blue]
	linestyles = [:solid, :solid, :dash, :solid]
	for (ischeme, (scheme, (K̃s, ωs))) in enumerate(dfs)
		linestyle = linestyles[ischeme]
		color = colors[ischeme]
		ωs = stack(ωs)
		ωs = sort(ωs, dims=1, by=imag)
		for i = 1:size(ωs, 1)
			lines!(ax, K̃s, imag.(ωs[i,:]); color, linestyle, label=scheme, linewidth=3)
		end
	end
	#axislegend(; merge=true)
	Legend(f[1,2], ax; merge=true)
	f
end

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

# ╔═╡ bc991b63-8e0e-4e37-9369-15f2e11e195e
md"""
## Diffusion Operator
"""

# ╔═╡ 3bf9a733-0328-460e-9988-114f95e24bb6
md"""
## Gradient Operator
"""

# ╔═╡ 986a7932-d8aa-4478-a631-2b028cbf1ed5
md"""
Small wavenumber approximation: $(@bind doapproxg CheckBox(default=true))
"""

# ╔═╡ ed499209-af9f-4539-bcb1-c47d93912208
md"""
## Divergence Operator
"""

# ╔═╡ 86daa502-9fb5-414a-9ea3-7294316abce1
md"""
## Vector Gradient Operator
"""

# ╔═╡ 7e1d7619-67ea-49bf-830f-2e63f224960d
md"""
Small wavenumber approximation: $(@bind doapproxv CheckBox(default=true))
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

# ╔═╡ 26c59f64-ad09-4ac6-a3ad-63f6372550a1
sys = GridOperatorAnalysis.construct_lcc_sys(bflow, dflow, a; f₀, g, 𝕂ᵘ, 𝕂ᵇ, hmt_scheme, hst_scheme, dissip_scheme);

# ╔═╡ 7f20375f-6491-47ff-bcc0-fb77baa8bedf
hds = let
	cpb = colpts[colpt_type(flow_t, :b)]
	pb =  ϵ*dflow.b[cpb[1], cpb[2], cpb[3]] # bflow.b[cpb[1], cpb[2], cpb[3]] .+
	hds = []
	for biH=1:dims(colpt_type(flow_t, :b))
		hd = TriB.Δ(colptidx(0,0,biH,Val(colpt_type(flow_t, :b))), cpb, pb)
		hd = substitute(hd, sqrt3_subs)
		push!(hds, hd)
	end
	hds
end;

# ╔═╡ 1af2939c-6795-4472-ae44-8c9466a159ab
lhds = [
	expand(taylor_coeff(expand(hd), ϵ, 1)) for hd in hds
];

# ╔═╡ 5ce14648-6fc3-4efe-b73d-755afe57eddc
hgs = let
	cpu⃗ = colpts[colpt_type(flow_t, :u⃗)]
	cpb = colpts[colpt_type(flow_t, :b)]
	pb  =  ϵ*dflow.b[cpb[1], cpb[2], cpb[3]]
	∇   = TriB.∇cv #gethst(grid_t, hst_scheme)
	hgs = []
	for uiH=1:dims(colpt_type(flow_t, :u⃗))
		b  = pb
		hg = ∇(colptidx(0,0,uiH,Val(colpt_type(flow_t, :u⃗))), cpb, b)
		hg = substitute(hg, GridOperatorAnalysis.sqrt3_subs)
		push!(hgs, hg)
	end
	hgs
end;

# ╔═╡ 8699ac4d-363d-4502-8353-5598e3e71591
lhgs = if colpt_type(flow_t, :u⃗) == :edge
	[expand(taylor_coeff(hgs[iH], ϵ, 1)) for iH=1:dims(colpt_type(flow_t, :u⃗))]
else
	[[expand(taylor_coeff(hgs[iH][iTH], ϵ, 1)) for iTH=1:2] for iH=1:dims(colpt_type(flow_t, :u⃗))]
end;

# ╔═╡ 5248fb8c-c5b7-4b60-a09e-30f7de30dc3d
hns = let
	cp = colpts[colpt_type(flow_t, :u⃗)]
	pu⃗  = colpt_type(flow_t, :u⃗) == :edge ? ϵ*dflow.u⃗[cp[1], cp[2], cp[3]] : [ϵ*dflow.u⃗[iTH, cp[1], cp[2], cp[3]] for iTH=1:2]
	∇ᵀ  = TriB.∇ᵀvc #gethst(grid_t, hst_scheme)
	hns = []
	for biH=1:dims(colpt_type(flow_t, :u⃗))
		u  = pu⃗
		hn = ∇ᵀ(colptidx(0,0,biH,Val(colpt_type(flow_t, :u⃗))), cp, u)
		hn = substitute(hn, GridOperatorAnalysis.sqrt3_subs)
		push!(hns, hn)
	end
	hns
end;

# ╔═╡ 9b9d9dc3-978a-42b0-a8e3-e934b27623ea
lhns = [
	expand(taylor_coeff(expand(hn), ϵ, 1)) for hn in hns
];

# ╔═╡ 21a6d07b-3ab4-46bc-b016-9d155e32fb0c
begin
	hvs = []
	for iH=1:dims(colpt_type(flow_t, :u⃗))
		cp = colpts[colpt_type(flow_t, :u⃗)]
		pu⃗ = colpt_type(flow_t, :u⃗) == :edge ? ϵ*dflow.u⃗[cp[1], cp[2], cp[3]] : [ϵ*dflow.u⃗[iTH, cp[1], cp[2], cp[3]] for iTH=1:2]
		∇ = TriB.∇cc
		push!(hvs, ∇(colptidx(0,0,iH,Val(colpt_type(flow_t, :u⃗))), cp, pu⃗))
	end
end;

# ╔═╡ 1932073b-e14c-40c3-bedc-706daf39e44c
lhvs = if colpt_type(flow_t, :u⃗) == :edge
	[expand(taylor_coeff(hvs[iH], ϵ, 1)) for iH=1:dims(colpt_type(flow_t, :u⃗))]
else
	[[expand(taylor_coeff(hvs[iH][iTH, jTH], ϵ, 1)) for iTH=1:2, jTH=1:2] for iH=1:dims(colpt_type(flow_t, :u⃗))]
end;

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

# ╔═╡ c4bdbd16-6527-4bca-808f-bb16ebb6e8c0
ftsys = GridOperatorAnalysis.fourier_transform_sys(Val(grid_t), sys; dflow, fflow, ϕ);

# ╔═╡ 7d618255-0b0c-408a-ae10-653383d6f4de
fsymbols = GridOperatorAnalysis.fouriersymbols(Val(grid_t), ftsys, fflow; k, l, f₀, g, N², M²=√(f₀^2*N²/Ri), β, θU, 𝕂ᵘ, 𝕂ᵇ, a, h, doapprox, dophasesubs);

# ╔═╡ f956d11d-37f5-4a7c-bfd9-471104e1d8a3
if colpt_type(flow_t, :u⃗) ≠ :edge
	let
		d = dims(colpt_type(flow_t, :u⃗))
		reshape(permutedims(fsymbols[:Gx], (2,1,4,3)), 2*d, 2*d)
	end
else
	fsymbols[:Gx]
end

# ╔═╡ c86c02aa-51a2-4af4-b6b0-13035262e909
if colpt_type(flow_t, :u⃗) ≠ :edge
	let
		d = dims(colpt_type(flow_t, :u⃗))
		fs = reshape(permutedims(fsymbols[:Gy], (2,1,4,3)), 2*d, 2*d)
		expand.(fs)
	end
else
	let
		fs = fsymbols[:Gy] #.|> simplifyexpand
		#simplify.(fs; rewriter=expanda)
	end
end

# ╔═╡ eeecfcd8-f994-4d0c-b84c-e1a4fbbf9d60
let
	d = dims(colpt_type(flow_t, :u⃗))
	if colpt_type(flow_t, :u⃗) ≠ :edge
		reshape(transpose(fsymbols[:A⁽ˣ⁾][:,:,1]), 2 * d, 1)
	else
		fs = substitute.(fsymbols[:A⁽ˣ⁾], Ref(Dict(a=>2/sqrt3 * h))) .|> simplifyexpand
		substitute.(fs, Ref(sqrt3_subs)) .|> simplifyexpand
	end
end

# ╔═╡ 9bcd1831-1391-466a-9cb8-8c547f82754f
let
	d = dims(colpt_type(flow_t, :u⃗))
	if colpt_type(flow_t, :u⃗) ≠ :edge
		expand.(reshape(transpose(fsymbols[:G][:,:,1]), 2 * d, 1))
	else
		#fs = substitute.(fsymbols[:G], Ref(Dict(a=>2/sqrt3 * h))) .|> simplifyexpand
		#substitute.(fs, Ref(sqrt3_subs)) .|> simplifyexpand
		expand.(fsymbols[:G])
	end
end

# ╔═╡ 05f47d93-061a-43c2-bdbe-f2edd5a3ae6f
fsymbols[:Γx]

# ╔═╡ 32ff2baf-1795-405e-be0f-ecd2398e78ce
substitute.(fsymbols[:Γy], Ref(Dict(a^2=>4//3*h^2))) .|> simplifyexpand

# ╔═╡ 8a580c44-1ed6-4cc1-91c2-8ff457d6d296
if colpt_type(flow_t, :u⃗) ≠ :edge
	substitute.(fsymbols[:Av⁽ˣ⁾][1,:,:], Ref(Dict(a=>2//sqrt3 * h))) .|> simplifyexpand
else
	fs = fsymbols[:Av⁽ˣ⁾]
	fs = substitute.(fs, Ref(Dict(a=>2/sqrt3 * h))) .|> simplifyexpand
	fs = substitute.(fs, Ref(sqrt3_subs)) * sqrt3 .|> simplifyexpand
	fs ./ sqrt3 .|> simplifyexpand
end

# ╔═╡ bdc0e050-28c9-44f5-acde-b1c63521d791
expand.(fsymbols[:Dᵇ])

# ╔═╡ 4b2d9f6d-4438-436a-91cb-47fb66270841
fhds = [
	fourier_transform_expression(biH, colpt_type(flow_t, :b), lhds[biH]; fflow, dflow, ϕ) for biH=1:dims(colpt_type(flow_t, :b))
];

# ╔═╡ 49c1fe6e-725f-4de1-922e-9d091d5b70e7
begin
	fhgs = []
	for iH=1:dims(colpt_type(flow_t, :u⃗))
		if colpt_type(flow_t, :u⃗) == :edge
			fhg = fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhgs[iH]; fflow, dflow, ϕ)
			push!(fhgs, fhg)
		else
			fhg = [fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhgs[iH][iTH]; fflow, dflow, ϕ) for iTH=1:2]
			push!(fhgs, fhg)
		end
	end		
end;

# ╔═╡ d04861c4-b170-49d4-84db-7a59bc5bafc4
fhns = [
	fourier_transform_expression(biH, colpt_type(flow_t, :b), lhns[biH]; fflow, dflow, ϕ) for biH=1:dims(colpt_type(flow_t, :b))
];

# ╔═╡ 8f8edfc5-c8a0-4fcc-a58b-d94e585082aa
begin
	fhvs = []
	for iH=1:dims(colpt_type(flow_t, :u⃗))
		if colpt_type(flow_t, :u⃗) == :edge
			fhv = fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhvs[iH]; fflow, dflow, ϕ)
			push!(fhvs, fhv)
		else
			fhv = [fourier_transform_expression(iH, colpt_type(flow_t, :u⃗), lhvs[iH][iTH, jTH]; fflow, dflow, ϕ) for iTH=1:2, jTH=1:2]
			push!(fhvs, fhv)
		end
	end		
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
	else
		if hmt_scheme == :MPAS
	        HexC.u⃗ᵀ∇
	    else
	        throw(ArgumentError("The scheme $hmt_scheme for the horizontal momentum transport is not implemented"))
	    end
	end
end

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
		b = pb# - a^2/24*TriC.ℳ(cpb, cpb, TriC.Δ(cpb, cpb, pb)) #
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
	Fhd = simplify.(Fhd ./ -K^2; expand=true)
end

# ╔═╡ 5401495e-95d7-4a67-b567-c812ff2daed5
sfhgs = let
	colpt_t   = colpt_type(flow_t, :u⃗)
	d         = dims(colpt_t)
	colpt_t_b = colpt_type(flow_t, :b)
	d_b       = dims(colpt_t_b)
	sfhgs = if colpt_t == :edge
		zeros(Complex{Symbolics.Num}, d, d_b)
	else
		zeros(Complex{Symbolics.Num}, 2, d, d_b)
	end
	for iH=1:d
		fhg = let
			fhg = expand.(fhgs[iH])
			fhg = if doapproxg # rewrite trig-functions by Taylor polynomial
				fhg = simplify.(fhg; rewriter=rtrig) .|> expand
				substitute.(fhg, Ref(Dict(sqrt3 => 2*h/a)))
			else
				fhg = substitute.(fhg, Ref(Dict(l => h/a * 2/sqrt3 * l)))
				substitute.(fhg, Ref(Dict(sqrt3 => 2*h/a)))
			end
			fhg = substitute.(fhg, Ref(Dict(le=>a, √(f₀^2*N²/Ri) => M²))) .|> expand
			fhg = substitute.(fhg, Ref(sqrt3_subs)) .|> expand
		end
		if colpt_t == :edge
			for jH = 1:d_b
				b = fflow.b[jH]
				sfhg = taylor_coeff(fhg, b, 1)
				sfhgs[iH, jH] = simplify(sfhg; expand=true)
			end
		else
			for iTH=1:2
				for jH=1:d_b
					b = fflow.b[jH]
					sfhg = taylor_coeff(fhg[iTH], b, 1)
					sfhgs[iTH, iH, jH] = simplify(sfhg; expand=true)
				end
			end
		end
	end
	simplify.(sfhgs; expand=true)
end;

# ╔═╡ cb54e7d8-81f9-4ccd-b04b-33e8c0c5642d
let
	sgu = sfhgs[:,1,1] #./ (im * [k; l])
	sgu = substitute.(sgu, Ref(Dict(a^2 => 4//3*h^2, a^4 => 16//9*h^4)))
	simplify.(sgu; expand=true)
end

# ╔═╡ cf8ee985-7fb5-489a-ae96-b475a4740ae4
sfhgs[1,:,1]

# ╔═╡ e7e6c299-c6be-4a7c-8e55-70c3061332e8
sfhns = let
	colpt_t   = colpt_type(flow_t, :u⃗)
	d         = dims(colpt_t)
	colpt_t_b = colpt_type(flow_t, :b)
	d_b       = dims(colpt_t_b)
	sfhns = if colpt_t == :edge
		zeros(Complex{Symbolics.Num}, d_b, d)
	else
		zeros(Complex{Symbolics.Num}, d_b, 2, d)
	end
	for iH=1:d_b
		fhn = let
			fhn = expand.(fhns[iH])
			fhn = if doapproxg # rewrite trig-functions by Taylor polynomial
				fhn = simplify.(fhn; rewriter=rtrig) .|> expand
			else
				fhn = substitute.(fhn, Ref(Dict(l => h/a * 2/sqrt3 * l)))
			end
			fhn = substitute.(fhn, Ref(Dict(le=>a, √(f₀^2*N²/Ri) => M²))) .|> expand
			fhn = substitute.(fhn, Ref(sqrt3_subs)) .|> expand
			substitute.(fhn, Ref(Dict(sqrt3 => 2*h/a)))
		end
		if colpt_t == :edge
			for jH = 1:d
				u = fflow.u⃗[jH]
				sfhn = taylor_coeff(fhn, u, 1)
				sfhns[iH, jH] = simplify(sfhn; expand=true)
			end
		else
			for jTH=1:2
				for jH=1:d
					u = fflow.u⃗[jTH, jH]
					sfhn = taylor_coeff(fhn, u, 1)
					sfhns[iH, jTH, jH] = simplify(sfhn; expand=true)
				end
			end
		end
	end
	simplify.(sfhns; expand=true)
end;

# ╔═╡ fb0cbfc6-a1b2-4ed3-b974-a184da507540
sfhns[1,1,:]

# ╔═╡ 2476af78-0079-4b3e-bdb4-5e58874d8380
isequal(sfhns[1,2,2], simplify(1//2*sfhgs[2,1,1]; expand=true))

# ╔═╡ f47ab6b0-e15d-49b8-b731-17dd973e4f63
sfhns[1,2,2], simplify(1//2*sfhgs[2,1,1]; expand=true)

# ╔═╡ 48926455-a6de-4351-a6d9-9b0074af0028
sfhvs = let
	colpt_t = colpt_type(flow_t, :u⃗)
	d       = dims(colpt_t)
	sfhvs = if colpt_t == :edge
		zeros(Complex{Symbolics.Num}, d, d)
	else
		zeros(Complex{Symbolics.Num}, d, 2, 2, d, 2)
	end
	for iH=1:d
		fhv = let
			fhv = expand.(fhvs[iH])
			fhv = if doapproxv # rewrite trig-functions by Taylor polynomial
				simplify.(fhv; rewriter=rtrig) .|> expand
			else
				substitute.(fhv, Ref(Dict(l => h/a * 2/sqrt3 * l)))
			end
			fhv = substitute.(fhv, Ref(Dict(le=>a, √(f₀^2*N²/Ri) => M²))) .|> expand
			fhv = substitute.(fhv, Ref(sqrt3_subs)) .|> expand
			substitute.(fhv, Ref(Dict(sqrt3 => 2*h/a)))
		end
		if colpt_t == :edge
			for jH = 1:d
				u = fflow.u⃗[jH]
				sfhv = taylor_coeff(fhmt, u, 1)
				sfhvs[iH, jH] = simplify(sfhv; expand=true)
			end
		else
			for iTH=1:2
				for ijTH=1:2
					for jH=1:d
						for jTH=1:2
							u = fflow.u⃗[jTH, jH]
							sfhv = taylor_coeff(fhv[iTH, ijTH], u, 1)
							sfhvs[iH, iTH, ijTH, jH, jTH] = simplify(sfhv; expand=true)
						end
					end
				end
			end
		end
	end
	sfhvs
end;

# ╔═╡ 84678921-87f9-4e3e-b571-25ebcf7ee6cc
let
	sv = sfhvs[:,2, 1,:,1]
	sv = substitute.(sv, Ref(Dict(a^2 => 4//3*h^2, a^4 => 16//9*h^4)))
	simplify.(sv; expand=true)
end

# ╔═╡ 7b864ab1-b61c-481c-849a-9824c9abc984
rpyt = let
	r = Symbolics.@acrule sin(~x)^2 + cos(~x)^2 => one(~x)
	SymbolicUtils.Prewalk(SymbolicUtils.PassThrough(r))
end

# ╔═╡ 3f4be713-595b-444a-9e76-9b3b04baed53
expanda = let
	function p(x::Int)
		x ≥ 2
	end
	r1 = Symbolics.@rule a^(~x::p) => 4//3 * a^(~x-2) * h^2
	r2 = Symbolics.@rule ~z * a^(~x::p) * h^(~y) => ~z * 4//3 * a^(~x-2) * h^(~y+2)
	r3 = SymbolicUtils.RestartedChain([r1, r2])
	SymbolicUtils.Postwalk(SymbolicUtils.PassThrough(SymbolicUtils.Fixpoint(r3)))
end

# ╔═╡ f457599b-08a8-4f60-bd4a-614a92091645
if colpt_type(flow_t, :u⃗) ≠ :edge
	let
		d = dims(colpt_type(flow_t, :u⃗))
		reshape(permutedims(fsymbols[:M], (2,1,4,3)), 2*d, 2*d)
	end
else
	let
		fs = simplify.(fsymbols[:M]; rewriter=expanda)
		fs = substitute.(fs, Ref(Dict(a=>2/sqrt3 * h))) * sqrt3 .|> simplifyexpand .|> expand
		fs ./ sqrt3 .|> simplifyexpand
	end
end

# ╔═╡ a1b6588c-bfbc-4578-8fe6-21a8b91c2312
let
	d = dims(colpt_type(flow_t, :u⃗))
	if colpt_type(flow_t, :u⃗) ≠ :edge
		reshape(transpose(fsymbols[:A⁽ʸ⁾][:,:,1]), 2 * d, 1)
	else
		simplify.(fsymbols[:A⁽ʸ⁾]; rewriter=expanda)
	end
end

# ╔═╡ 72c3124a-ef7c-4d15-aac4-9cae947c201a
let
	if colpt_type(flow_t, :u⃗) ≠ :edge
		d = dims(colpt_type(flow_t, :u⃗))
		fs = reshape(permutedims(fsymbols[:Dᵘ], (2,1,4,3)), 2*d, 2*d) * h^4 .|> simplifyexpand
		simplifyexpand.(fs ./ h^4)
	else
		fs = fsymbols[:Dᵘ]
		fs = simplify.(fs; rewriter=expanda) * h^4 .|> simplifyexpand .|> expand
		fs ./ h^4
	end
end

# ╔═╡ 564bd7ad-fc9f-4d69-a9c3-54d05f05fb5f
if colpt_type(flow_t, :u⃗) ≠ :edge
	substitute.(fsymbols[:Av⁽ʸ⁾][1,:,:], Ref(Dict(a=>2//sqrt3 * h))) .|> simplifyexpand
else
	simplify.(fsymbols[:Av⁽ʸ⁾]; rewriter=expanda)
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
if doapprox && colpt_type(flow_t, :u⃗) ≠ :edge
	svs = let
		Gx = fsymbols[:Gx][1,:,1,:] ./ (im*k) .|> simplifyexpand
		compute_symbolic_eigenvals(Gx)
	end
	svs = Symbolics.simplify.(svs; expand=true)
else
	"Computation of symbolic eigenvalues not implemented in this case."
end

# ╔═╡ f1feff04-e468-4909-af4c-9bc5a032f407
K̃s, ωs = let
	@variables K̃ K
	K̃s = range(1e-20, π * 2/√3, 100) #floatmin(Float64)
	ωs = []
	if doapprox && colpt_type(flow_t, :u⃗) ≠ :edge
		for sv in svs
			sv = substitute(sv, Dict(
				k => K̃ / h * cos(sθU),
				l => K̃ / h * sin(sθU),
				sqrt3 => √3)) |> simplifyexpand
			vfunc = Symbolics.build_function(sv, K̃; expression=Val(false))
			push!(ωs, vfunc.(K̃s))
		end
	else
		A = if colpt_type(flow_t, :u⃗) ≠ :edge
			cos(sθU) * reshape(fsymbols[:Gx], 4, 4) .+ sin(sθU)  * reshape(fsymbols[:Gy], 4, 4) .|> simplifyexpand
		else
			cos(sθU) * fsymbols[:Gx] .+ sin(sθU)  * fsymbols[:Gy] .|> simplifyexpand
		end
		A = substitute.(A, Ref(invphase_subs))
		A = substitute.(A, Ref(Dict(
			a => 2/√3 * h, sqrt3=>√3, k => K̃ / h * cos(sθU), l => K̃ / h * sin(sθU)
		))) .|> simplifyexpand
		A = substitute.(A, Ref(Dict(
			θU => sθU, 
			h => K̃ / K, 
			f₀ => -1e-4,
			g  => 1e9,
			N² => 1e-6,
			Ri => 1,
			M² => √((-1e-4)^2*1e-6/1),
			z  => 0,
			β  => 0
		))) .|> simplifyexpand
		A = A ./ (im * K) .|> simplifyexpand
		Afunc = Symbolics.build_function(A, K̃; expression=Val(false))[1]
		for k̃ in K̃s
			(; values) = eigen(Complex{Float64}.(Afunc(k̃)))
			push!(ωs, values)
		end
	end
	(K̃s, ωs)
end

# ╔═╡ 54131367-66b2-4bc7-8820-0b6a5fc19ae3
let
	f = Figure(; fontsize=36)
	ax = Axis(f[1,1];
			 xlabel = L"\text{wavenumber}~K / h^{-1}",
			 ylabel = L"ω / i\bar{U}K",
			 limits = (0, π * 2/√3, -2.1, 1.1), 
			 aspect = 1,
			 )
	if doapprox && colpt_type(flow_t, :u⃗) ≠ :edge
		for ω in ωs
			lines!(ax, K̃s, real.(ω), linewidth=3)
		end
	else
		_ωs = stack(ωs)
		_ωs = sort(_ωs, dims=1, by=real)
		for i = 1:size(_ωs, 1)
			lines!(ax, K̃s, real.(_ωs[i,:]); linewidth=3)
		end
		#for (k̃, ω) in zip(K̃s, ωs)
		#	for _ω in ω
		#		scatter!(ax, k̃, imag(_ω); color=:blue)
		#	end
		#end
	end
	#axislegend()
	f
end

# ╔═╡ 6dab1873-d738-4e41-9d16-0f83fb168bd4
if doapprox
	svbs = let
		@variables K
		fs = cos(θU) * fsymbols[:Γx] + sin(θU) * fsymbols[:Γy]
		fs = substitute.(fs, Ref(Dict(
			k => K̃ / h * cos(θU), l => K̃ / h * sin(θU)
		))) .|> expand
		fs = substitute.(fs, Ref(Dict(
			a => K̃ / K * 2//sqrt3, h => K̃ / K
		))) .|> expand
		fs = substitute.(fs, Ref(Dict(
			sin(θU)^2=>1-cos(θU)^2,
			sin(θU)^4=>(1-cos(θU)^2)^2,
			sin(θU)^6=>(1-cos(θU)^2)^3
		))) .|> simplifyexpand
		fs = fs ./ (im * K)
		svbs = compute_symbolic_eigenvals(fs)
		substitute.(svbs, Ref(Dict(
			sin(θU)^2=>1-cos(θU)^2,
			sin(θU)^4=>(1-cos(θU)^2)^2,
			sin(θU)^6=>(1-cos(θU)^2)^3, 
			sin(θU)^8=>(1-cos(θU)^2)^4
		)))
	end
	svbs = Symbolics.simplify.(svbs; expand=true)
else
	"Computation of symbolic eigenvalues not implemented in this case."
end

# ╔═╡ 6a120e6b-ced8-414a-a56b-3b49ea21c8fa
let
	K̃s, ωs = let
		K̃s = range(1e-20, π, 100) #floatmin(Float64)
		ωs = []
		if doapprox
			for sv in svbs
				sv = substitute.(sv, Ref(Dict(θU => sθU)))
				vfunc = Symbolics.build_function(sv, K̃; expression=Val(false))
				push!(ωs, vfunc.(K̃s))
			end
		else
			@variables K
			A = cos(tθU) * fsymbols[:Γx] .+ sin(tθU)  * fsymbols[:Γy] .|> simplifyexpand
			A = substitute.(A, Ref(invphase_subs))
			A = substitute.(A, Ref(Dict(
				a => 2/√3 * h, sqrt3=>√3, k => K̃ / h * cos(tθU), l => K̃ / h * sin(tθU)
			))) .|> simplifyexpand
			A = substitute.(A, Ref(Dict(
				θU => sθU, 
				h => K̃ / K, 
				f₀ => -1e-4,
				g  => 1e9,
				N² => 1e-6,
				Ri => 1,
				M² => √((-1e-4)^2*1e-6/1),
				z  => 0,
				β  => 0
			))) .|> simplifyexpand
			A = A ./ (im * K) .|> simplifyexpand
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
			 xlabel = L"\text{wavenumber}~K / h^{-1}",
			 ylabel = L"ω / i\bar{U}K",
			 limits = (0, π, -2.1, 1.1), 
			 aspect = 1,
			 )
	if doapprox
		for ω in ωs
			lines!(ax, K̃s, real.(ω), linewidth=3)
		end
	else
		ωs = stack(ωs)
		ωs = sort(ωs, dims=1, by=real)
		for i = 1:size(ωs, 1)
			lines!(ax, K̃s, real.(ωs[i,:]); linewidth=3)
		end
		#for (k̃, ω) in zip(K̃s, ωs)
		#	for _ω in ω
		#		scatter!(ax, k̃, imag(_ω); color=:blue)
		#	end
		#end
	end
	#axislegend()
	f
end

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
# ╠═d60c197f-28db-4834-b78f-aab7f23536f2
# ╟─59056484-9c6b-48ad-b741-2bc294d5cc6f
# ╟─0943c748-9fee-401e-851e-de2327d39706
# ╟─9abfbc35-bd37-4554-949c-27cd2bdfa1a7
# ╠═cfb78fd5-4435-4c0d-8cf2-1b5fffc2f91e
# ╠═d2b72544-78fb-4a0b-b169-aa50ebaeb31d
# ╠═d32f2049-e2d4-4073-865f-585cc80d4d07
# ╠═4a5bd3c3-bcfe-498e-8fa7-328cd12cea00
# ╟─428cc192-03d0-497c-aee7-f7015c4c09d5
# ╠═fb5e7874-c067-4769-aefa-260fd3eca00b
# ╟─a697761a-3f01-46ac-ae31-4665bc7b6e47
# ╠═2663542e-c719-42f0-bf00-790e8aa211fe
# ╠═1b94a213-9311-4dad-bc3e-0ac2dcb71e14
# ╟─430cca4e-08b6-4d57-a4dc-422dcf8cc92e
# ╠═8ae523f0-1eeb-4d62-b089-502c52dd1277
# ╟─da53d71b-e5b3-4ce0-8ce8-1a0449c4fd32
# ╠═26c59f64-ad09-4ac6-a3ad-63f6372550a1
# ╟─77eb6488-c26a-416a-839b-1e617f0cdd50
# ╠═7741fc90-ad1f-4dab-baa5-bbf373592e13
# ╠═780abc6b-d9ab-40d1-9a01-d12b1d3fc3ae
# ╠═c4bdbd16-6527-4bca-808f-bb16ebb6e8c0
# ╟─0f1cbe99-8b11-426d-8d28-77b6dffef277
# ╟─77fe42fd-ddbb-495a-979e-0a852600a2a6
# ╠═7d618255-0b0c-408a-ae10-653383d6f4de
# ╟─e287982b-5fea-4d73-a0ce-d2627635ddb2
# ╟─3b152922-e343-42a1-b863-38814a80b0d7
# ╟─33e1ab03-9e1a-4014-84ed-1ad393f0445e
# ╟─f9adc84b-12c3-4a1c-8400-6ede0108c9e9
# ╟─39031d65-863d-48cf-aeaf-520417eaf071
# ╟─286c65de-46f4-462d-903c-ad850bb1a3b1
# ╠═f71432b9-cb86-4f5e-a6ba-4ad443140292
# ╟─583e619a-76da-4f63-8377-b5eac40c96af
# ╟─729a626e-376c-4ad9-96e8-f2e641f4d2a4
# ╠═f956d11d-37f5-4a7c-bfd9-471104e1d8a3
# ╟─30725e26-6949-4605-a275-3f72eca04ece
# ╠═c86c02aa-51a2-4af4-b6b0-13035262e909
# ╠═0e2ee007-f695-408b-9a8f-3619b72c137f
# ╟─33dcefea-8d67-4543-a241-3192ad3252d2
# ╟─f457599b-08a8-4f60-bd4a-614a92091645
# ╟─d40d6c8c-f814-4c98-989c-a9b3d7913c0c
# ╟─eeecfcd8-f994-4d0c-b84c-e1a4fbbf9d60
# ╟─801f824d-9fc4-431c-bbbe-a72ed0d3b3ac
# ╟─a1b6588c-bfbc-4578-8fe6-21a8b91c2312
# ╟─7d3ce7a0-eb0b-4801-8dec-84ba360dd477
# ╟─9bcd1831-1391-466a-9cb8-8c547f82754f
# ╟─4a18efe0-0d3e-481c-94a4-d1b3e4d56df3
# ╠═72c3124a-ef7c-4d15-aac4-9cae947c201a
# ╟─46afe77d-f904-40e9-92bb-cab83aeea51e
# ╠═16535b97-0c7e-4c14-9e71-604b0e0c9a3d
# ╟─f1feff04-e468-4909-af4c-9bc5a032f407
# ╠═5708950e-80fb-4296-adfa-e2db01284d06
# ╠═4ffe2716-9a66-4ca5-bf4b-367e85399a5b
# ╟─0890ade4-59b2-45cf-a71d-35f918c71883
# ╟─54131367-66b2-4bc7-8820-0b6a5fc19ae3
# ╠═57b612b5-afb6-44db-ad3f-6c9882afae14
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
# ╟─77960a21-8754-427c-875b-e6578e0ed62f
# ╠═05f47d93-061a-43c2-bdbe-f2edd5a3ae6f
# ╟─da3f6839-042f-4617-abbf-fb5a18fe772a
# ╠═32ff2baf-1795-405e-be0f-ecd2398e78ce
# ╟─05fb2f19-1b43-4e7a-84c5-1ea8a71b9284
# ╟─f4d63a77-6141-41df-bcdf-f7547d4e3cf5
# ╟─8a580c44-1ed6-4cc1-91c2-8ff457d6d296
# ╟─fd7bfb37-1bfb-41ff-8dd1-0c566156693d
# ╟─564bd7ad-fc9f-4d69-a9c3-54d05f05fb5f
# ╟─c6321e6c-e4ce-45fa-b884-69f14535db2d
# ╠═bdc0e050-28c9-44f5-acde-b1c63521d791
# ╟─26fd8b72-c736-43db-9d32-f76cbe5dd901
# ╟─6dab1873-d738-4e41-9d16-0f83fb168bd4
# ╠═5d83c523-626f-4b06-85ee-8b71423d398a
# ╠═30620b89-ced4-4353-bf3f-6d9c58794aa4
# ╠═932e3e09-c15f-4b0e-944a-1eecbcf378d7
# ╟─6a120e6b-ced8-414a-a56b-3b49ea21c8fa
# ╟─23a4913d-3686-42a0-a2ce-c34b298fa59f
# ╟─0e90c7a9-fc7e-4890-b6e0-32ca12e8ae7d
# ╟─9b3f3f72-ffb8-477d-b1a5-4a14d64d36e8
# ╟─5361ebe0-a415-416c-b57f-051b0e6382f6
# ╟─bc991b63-8e0e-4e37-9369-15f2e11e195e
# ╠═7f20375f-6491-47ff-bcc0-fb77baa8bedf
# ╠═1af2939c-6795-4472-ae44-8c9466a159ab
# ╠═4b2d9f6d-4438-436a-91cb-47fb66270841
# ╠═27dbf7f2-14db-40ac-9e59-d2cfb3864a8e
# ╠═04d312fc-e353-4660-82ad-3dd048d79551
# ╠═dcbd2629-eb5d-42b5-9234-d132106f5747
# ╟─3bf9a733-0328-460e-9988-114f95e24bb6
# ╠═5ce14648-6fc3-4efe-b73d-755afe57eddc
# ╠═8699ac4d-363d-4502-8353-5598e3e71591
# ╠═49c1fe6e-725f-4de1-922e-9d091d5b70e7
# ╟─986a7932-d8aa-4478-a631-2b028cbf1ed5
# ╠═5401495e-95d7-4a67-b567-c812ff2daed5
# ╠═cb54e7d8-81f9-4ccd-b04b-33e8c0c5642d
# ╟─ed499209-af9f-4539-bcb1-c47d93912208
# ╠═5248fb8c-c5b7-4b60-a09e-30f7de30dc3d
# ╠═9b9d9dc3-978a-42b0-a8e3-e934b27623ea
# ╠═d04861c4-b170-49d4-84db-7a59bc5bafc4
# ╠═e7e6c299-c6be-4a7c-8e55-70c3061332e8
# ╠═fb0cbfc6-a1b2-4ed3-b974-a184da507540
# ╠═2476af78-0079-4b3e-bdb4-5e58874d8380
# ╠═f47ab6b0-e15d-49b8-b731-17dd973e4f63
# ╟─86daa502-9fb5-414a-9ea3-7294316abce1
# ╠═21a6d07b-3ab4-46bc-b016-9d155e32fb0c
# ╠═1932073b-e14c-40c3-bedc-706daf39e44c
# ╠═8f8edfc5-c8a0-4fcc-a58b-d94e585082aa
# ╟─7e1d7619-67ea-49bf-830f-2e63f224960d
# ╠═48926455-a6de-4351-a6d9-9b0074af0028
# ╠═84678921-87f9-4e3e-b571-25ebcf7ee6cc
# ╠═cf8ee985-7fb5-489a-ae96-b475a4740ae4
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
# ╠═3f4be713-595b-444a-9e76-9b3b04baed53
# ╟─d6ebd665-16c7-4850-8302-80356de08639
# ╠═c9a610d4-ecd7-4798-8b5d-88c96ad6afae
# ╠═f0e541cf-b2f5-40f5-8379-b1f8ef1a4fe8
# ╠═e80b16a7-315f-4b3f-b5ce-4cd73b6a3fdb
# ╠═4c17d92e-cdf7-442f-9087-102ff8109a8a
# ╠═dc22802f-307d-4bf1-a10d-3e1b0180783e
# ╠═e667e7fd-ab40-49c9-b8ca-a4aff155e833
