#------------------------------geometry-------------------------------
const sqrt3 = Symbolics.variable(:sqrt3)
const le    = Symbolics.variable(:le)
const lê    = Symbolics.variable(:lê)
const Ac    = Symbolics.variable(:Ac)
const Av    = Symbolics.variable(:Av)

const e = @syms e₁::Int e₂::Int e₃::Int
const c = @syms c₁::Int c₂::Int c₃::Int
const c̃ = @syms c̃₁::Int c̃₂::Int c̃₃::Int
const v = @syms v₁::Int v₂::Int
const ṽ = @syms ṽ₁::Int ṽ₂::Int

const n⃗ec, n⃗ec̃, n⃗eĉ = @variables n⃗ec[1:2] n⃗ec̃[1:2] n⃗eĉ[1:2]
const n⃗ev, n⃗eṽ, n⃗ev̂ = @variables n⃗ev[1:2] n⃗eṽ[1:2] n⃗ev̂[1:2]
const d = Symbolics.variable(:d)

function dot(a,b)
    substitute(a[1] * b[1] + a[2] * b[2], Dict([sqrt3 => √3]))
end

const n⃗e = ([ sqrt3//2;  1//2], [-sqrt3//2;  1//2], Num[ 0 ;-1//1])
const t⃗e = ([-1//2;  sqrt3//2], [-1//2; -sqrt3//2], Num[ 1//1;  0]) # e⃗₃ × n⃗e

const bb = le * Num[-1//2 -1; sqrt3//2 0]

_Ac(_le) = 1//2 * _le * √3/2 * _le
_Av(_le) = 2 * _Ac(_le)
_lê(le) = 2//3 * √3/2 * le

function EC(c)
    ifelse(c[3] == 1,
           [((c[1]+1, c[2]-1, 1), n⃗e[1]), ((c[1], c[2], 2), n⃗e[2]), ((c[1], c[2], 3), n⃗e[3])],
           [((c[1], c[2], 1), -n⃗e[1]), ((c[1]-1, c[2], 2), -n⃗e[2]), ((c[1], c[2], 3), -n⃗e[3])]
           )
end

function CE(e)
    ifelse(e[3] == 1,
           [((e[1]-1, e[2]+1, 1), n⃗e[1]), ((e[1], e[2], 2), -n⃗e[2])],
           ifelse(e[3] == 2,
                  [((e[1], e[2], 1), n⃗e[2]), ((e[1]+1, e[2], 2), -n⃗e[2])],
                  [((e[1], e[2], 1), n⃗e[3]), ((e[1], e[2], 2), -n⃗e[3])]
                  )
           )
end

function VE(e)
    ifelse(e[3] == 1,
           [((e[1], e[2]), -t⃗e[1]), ((e[1]-1, e[2]),  t⃗e[1])],
           ifelse(e[3] == 2,
                  [((e[1], e[2]), -t⃗e[2]), ((e[1]+1, e[2]-1), t⃗e[2])],
                  [((e[1], e[2]), t⃗e[3]), ((e[1], e[2]-1), -t⃗e[3])]
                  )
           )
end

function EV(v)
    [
        ((v[1], v[2], 1), -t⃗e[1]), ((v[1]+1, v[2], 1), t⃗e[1]),
        ((v[1], v[2], 2), -t⃗e[2]), ((v[1]-1, v[2]+1, 2), t⃗e[2]),
        ((v[1], v[2], 3),  t⃗e[3]), ((v[1], v[2]+1, 3), -t⃗e[3])
    ]
end

function VC(c)
    ifelse(c[3] == 1,
           [(c[1], c[2]), (c[1]+1, c[2]-1), (c[1], c[2]-1)],
           [(c[1], c[2]), (c[1]-1, c[2]), (c[1], c[2]-1)]
           )
end

function CV(v)
    return [
        (v[1], v[2], 1), (v[1], v[2]+1, 1), (v[1]-1, v[2]+1, 1),
        (v[1], v[2], 2), (v[1], v[2]+1, 2), (v[1]+1, v[2], 2)
    ]
end

function VEup(e, n⃗ev)
    ifelse(e[3] == 1,
           ifelse(dot(t⃗e[1], n⃗ev) < 0,
                  [((e[1]+1,e[2]), -1), ((e[1],e[2]), 1)],
                  [((e[1]-2,e[2]), -1), ((e[1]-1,e[2]), 1)]
                  ),
           ifelse(e[3] == 2,
                  ifelse(dot(t⃗e[2], n⃗ev) < 0,
                         [((e[1]-1,e[2]+1), -1), ((e[1],e[2]), 1)],
                         [((e[1]+2,e[2]-2), -1), ((e[1]+1,e[2]-1), 1)]
                         ),
                  ifelse(dot(t⃗e[3], n⃗ev) > 0,
                         [((e[1], e[2]+1), -1), ((e[1], e[2]), 1)],
                         [((e[1],e[2]-2), -1), ((e[1],e[2]-1), 1)])
                  )
           )
end

#------------------------------operators------------------------------

# gradient operators


function ∇cv(_c, _v, p)
    ve = ∑((_v, n⃗ev), VE(e), p // 2)
    [ 1//Ac * ∑((e, n⃗ec), EC(_c), ve * le * n⃗ec[iTH] ) for iTH = 1:2]
end

# implement ∇ec(e,c,F,n⃗ec) = [ ∑((c, d), CE(e, n⃗ec), d * F[jTH]) for jTH=1:2 ]?
# or d⃗ᵀ∇(e,c,F) = [ ∑((c,d), CE(e), d * F[jTH]) for jTH=1:2 ] and CE(e) are ordered according to d⃗ₑ = l⃗ₑ × k⃗, note n⃗ec n⃗ecᵀ∇ = d⃗ₑ d⃗ₑᵀ∇ is actually independent of c
# function ∇cc(_cout, _cin, F)
#     dF = [ ∑((_cin, d), CE(e, n⃗ec), d * F[jTH]) for jTH=1:2 ]
#     [∑((e, n⃗ec), EC(_cout), n⃗ec[iTH] * dF[jTH]) for iTH=1:2, jTH=1:2]
# end

function ∇cc(_cout, _cin, F⃗)
    dF = [ ∑((_cin, n⃗ec̃), CE(e), dot(n⃗ec, n⃗ec) * F⃗[jTH]) for jTH=1:2 ]
    [∑((e, n⃗ec), EC(_cout), n⃗ec[iTH] * dF[jTH]) for iTH=1:2, jTH=1:2]
end

# divergence operators

function ∇ᵀvc(_v, _c, F⃗)
    fe = ∑((_c, n⃗ec), CE(e), lê//2 * dot(n⃗ev, F⃗))
    1//Av * ∑((e, n⃗ev), EV(_v), fe)
end

# curl operators

function curl_vc(_v, _c, u⃗)
    s = ∑((_c, n⃗ec), CE(e), 1//2 * (-n⃗ev[2] * u⃗[1] + n⃗ev[1] * u⃗[2]))
    1//Av * ∑((e, n⃗ev), EV(_v), le * s)
end

# average operators

function av_vc(_v, _c, u)
    ∑(_c, CV(_v), u * Ac // (3 * Av))
end

function av_cv(_c, _v, p)
    ∑(_v, VC(_c), p * Av // (6 * Ac))
end

# horizontal momentum advection

function ∇ᵀcc(_cout, _cin, u⃗ad, u⃗tr)
    [av_cv(_cout, v, ∇ᵀvc(v, _cin, u⃗ad .* u⃗tr[iTH])) for iTH=1:2]
end

# function u⃗ᵀ∇(_cout, _cin, u⃗tr, u⃗ad)
#     ĉ = @syms ĉ₁::Int ĉ₂::Int ĉ₃::Int
#     F = ∑(_cin, CE(e), le//2 * d * n⃗ec' * u⃗ad)
#     Ru = [u⃗tr[iTH] + le//2 * d * n⃗ec' * ∇(c, _cin, u⃗tr[iTH]) for iTH=1:2]
#     Fu = ∑((_cin,d), CE(e), d * 0.5*(F + abs(F)) * ( Ru[iTH] - substitute(u⃗tr[iTH], Dict( _cin .=> _cout)) )) 
#     [1//Ac * ∑((e, n⃗ec), EC(_cout), Fu) for iTH=1:2]
# end

function u⃗ᵀ∇(_cout, _cin, u⃗tr, u⃗ad) # c ≡ _cout, ĉ ≡ _cin
    F = ∑((_cin, n⃗eĉ), CE(e), le//2 * dot(n⃗ec̃, u⃗ad))
    Ru = [evalat(c̃, _cin, u⃗tr[iTH]) + le//2 * dot(n⃗ec̃, ∇(c̃, _cin, u⃗tr[iTH])) for iTH=1:2]
    Fu = ∑((c̃, n⃗ec̃), CE(e), dot(n⃗ec, n⃗ec̃) * 0.5*(F + abs(F)) * ( Ru[iTH] - evalat(_cout, _cin, u⃗tr[iTH]) )) 
    [1//Ac * ∑((e, n⃗ec), EC(_cout), Fu) for iTH=1:2]
end

function u⃗ᵀ∇_vector_invariant(_cout, _cin, u⃗tr, u⃗ad)
    u⃗⊥ = [-u⃗tr[2]; u⃗tr[1]]
    ∇KE = ∇cv(_cout, v, sum([av_vc(v, _cin, u⃗tr[jTH]) * av_vc(v, _cin, u⃗ad[jTH]) for jTH=1:2]))
    ω = av_cv(_cout, v, curl_vc(v, _cin, u⃗tr))
    [ω * u⃗⊥[iTH] + ∇KE[iTH] for iTH=1:2]
end

# function ∇ᵀ(_vout, _c, _vin, u⃗, b; γ=3//4)
#     d̃ = Symbolics.variable(:d̃; T=Real)
#     qe = ∑(_c, CE(e), le//(2*sqrt3) * n⃗ev' * u⃗)
#     ∇beᵘ = ∑((_vin, d̃), VEup(e,n⃗ev), d̃ * b//le)
#     ∇beᶜ = ∑((_vin, d̃), VEce(e,n⃗ev), d̃ * b//le)
#     ∇beᵈ = ∑((_vin, d̃), VEdo(e,n⃗ev), d̃ * b//le)
#     ∇be⁺ = 2//3 * ∇beᶜ + 1//3 * ∇beᵘ
#     ∇be⁻ = 2//3 * ∇beᶜ + 1//3 * ∇beᵈ
#     be⁺ = b + le//2 * ∇be⁺ 
#     be⁻ = b - le//2 * ∇be⁻
#     be = ifelse(d == 1, be⁺, be⁻)
#     fe = ∑((_vin, d), VEce(e,n⃗ev), (qe + d * (1-γ) * abs(qe)) * be)
#     1//Av * ∑((e, n⃗ev), EV(_vout), fe)
# end

function ∇ᵀ(_vout, _c, _vin, u⃗, b; γ=3//4) # v ≡ _vout
    F = ∑((_c, n⃗ec), CE(e), lê//2 * dot(n⃗ev, u⃗))
    ∇beᶜ = 1//le * ∑((_vin, n⃗ev̂), VE(e), dot(n⃗eṽ, n⃗ev̂) * b)
    #∇beᵘ = n⃗eṽ' * ∇cv(cup(ṽ, n⃗eṽ), _vin, b)
    ∇beᵘ = 1//le * ∑((_vin, d), VEup(e,n⃗ev), d * b)
    ∇be  = 2//3 * ∇beᶜ + 1//3 * ∇beᵘ
    Rbe  = evalat(ṽ, _vin, b) + le//2 * ∇be
    1//Av * ∑((e, n⃗ev), EV(_vout), ∑((ṽ, n⃗eṽ), VE(e), dot(n⃗ev, n⃗eṽ) * (F + (1-γ) * abs(F))/2 * Rbe))
end

const nS = 3

# function bous_sys(cout, vout, cin, vin, u⃗, w, η, p, b, t, z)
#     @variables f₀
#     @variables t z
#     @variables (u(t,z))[1:2,-nS:nS,-nS:nS,-nS:nS] (b(t,z))[-nS:nS,-nS:nS] (η(t))[-nS:nS, -nS:nS]
#     @variables (w(t,z))[-nS:nS, -nS:nS] (p(t,z))[-nS:nS, -nS:nS]
#     vin  = @syms vin₁::Int vin₂::Int
#     vout = @syms vout₁::Int vout₂::Int
#     cin  = @syms cin₁::Int cin₂::Int cin₃::Int
#     cout = @syms cout₁::Int cout₂::Int cout₃::Int
#     ∂ₜ = Differential(t)
#     ∂₃ = Differential(z)

#     u⃗ = [u[iTH, cin[1], cin[2], cin[3]] for iTH=1:2]
#     b = b[vin[1], vin[2]]
#     η = η[vin[1], vin[2]]
#     w = w[vin[1], vin[2]]
#     p = p[vin[1], vin[2]]
    

#     eqs = [
#         ∇ᵀvc(vout, cin, u⃗) + ∂₃(evalat(vout, vin, w)) ~ 0,
#         ∂ₜ.(evalat(cout, cin, u⃗)) .+ ∇ᵀcc(cout, cin, u⃗, u⃗) .+ ∂₃.(evalat(cout, cin, u⃗) .* av_cv(cout, vin, w)) .+ f₀ * [0 -1; 1 0] * evalat(cout, cin, u⃗) .~ -∇cv(cout, vin, p) .- ∇cv(cout, vin, η),
#         ∂₃(evalat(vout, vin, p)) - evalat(vout, vin, b) ~ 0, 
#         ∂ₜ(evalat(vout, vin, b)) + ∇ᵀ(vout, cin, vin, u⃗, b) + ∂₃(evalat(vout, vin, w * b)) ~ 0,
#         ∂ₜ(evalat(vout, vin, η)) + ∇ᵀvc(vout, cin, ∫dz.(u⃗)) ~ 0
#     ]
#     eqs = evalat((0,0), vout, eqs)
#     eqs = evalat((0,0,1), cout, eqs)
# end

@kwdef struct BoussinesqSystem
    momentum_transport_eqs::Vector{Equation}
    hydrostatic_balance_eq::Equation
    continuity_eq::Equation
    buoyancy_transport_eq::Equation
    surface_elevation_eq::Equation
end

function BoussinesqSystem(cout, vout, cin, vin, t, z, u⃗, w, η, p, b, ∫∇ᵀu⃗dz, f₀)
    ∂ₜ = Differential(t)
    ∂₃ = Differential(z)

    momentum_transport_eqs = collect(∂ₜ.(evalat(cout, cin, u⃗)) .+ ∇ᵀcc(cout, cin, u⃗, u⃗) .+ ∂₃.(evalat(cout, cin, u⃗) .* av_cv(cout, vin, w)) .+ f₀ * [0 -1; 1 0] * evalat(cout, cin, u⃗) .~ -∇cv(cout, vin, p) .- ∇cv(cout, vin, η))
    hydrostatic_balance_eq = ∂₃(evalat(vout, vin, p)) - evalat(vout, vin, b) ~ 0
    continuity_eq          = ∇ᵀvc(vout, cin, u⃗) + ∂₃(evalat(vout, vin, w)) ~ 0
    buoyancy_transport_eq  = ∂ₜ(evalat(vout, vin, b)) + ∇ᵀ(vout, cin, vin, u⃗, b) + ∂₃(evalat(vout, vin, w * b)) ~ 0
    surface_elevation_eq   = ∂ₜ(evalat(vout, vin, η)) + evalat(vout, vin, ∫∇ᵀu⃗dz) ~ 0
    
    BoussinesqSystem(momentum_transport_eq, hydrostatic_balance_eq, continuity_eq, buoyancy_transport_eq, surface_elevation_eq)
end

linearize_eq(eq, ϵ) = Symbolics.taylor_coeff(simplify(expand_derivatives(eq.lhs)), ϵ, 1) ~ Symbolics.taylor_coeff(simplify(expand_derivatives(eq.rhs)), ϵ, 1)

function linearize_bous_sys(t, z; u⃗̄, b̄, η̄, w̄, p̄, du⃗, db, dη, dw, dp, f₀) #∫∇ᵀu⃗̄dz, ∫∇ᵀdu⃗dz
    vin  = @syms vin₁::Int vin₂::Int
    vout = @syms vout₁::Int vout₂::Int
    cin  = @syms cin₁::Int cin₂::Int cin₃::Int
    cout = @syms cout₁::Int cout₂::Int cout₃::Int
    @variables ϵ
 
    u⃗      = [u⃗̄[iTH, cin[1], cin[2], cin[3]] + ϵ * du⃗[iTH, cin[1], cin[2], cin[3]] for iTH=1:2]
    b      = b̄[vin[1], vin[2]] + ϵ * db[vin[1], vin[2]]
    η      = η̄[vin[1], vin[2]] + ϵ * dη[vin[1], vin[2]]
    w      = w̄[vin[1], vin[2]] + ϵ * dw[vin[1], vin[2]]
    p      = p̄[vin[1], vin[2]] + ϵ * dp[vin[1], vin[2]]
    ∫∇ᵀu⃗dz = ∫∇ᵀu⃗̄dz[vin[1], vin[2]] + ϵ * ∫∇ᵀdu⃗dz[vin[1], vin[2]]

    bous_sys = BoussinesqSystem(cout, vout, cin, vin, t, z, u⃗, w, η, p, b, ∫∇ᵀu⃗dz, f₀)
    # make assumption that linearized system has horizontally constant coefficients; how to assert this?

    momentum_transport_eqs = Equation[]
    for iTH=1:2
        for iHc=1:2
            eq = evalat((0,0,iHc), cout, bous_sys.momentum_transport_eqs[iTH])
            eq = linearize_eq(eq, ϵ)
            push!(momentum_transport_eqs, eq)
        end
    end

    hydrostatic_balance_eq = let
        eq = evalat((0,0), vout, bous_sys.hydrostatic_balance_eq)
        eq = linearize_eq(eq, ϵ)
    end

    continuity_eq = let
        eq = evalat((0,0), vout, bous_sys.continuity_eq)
        eq = linearize_eq(eq, ϵ)
    end

    buoyancy_transport_eq = let
        eq = evalat((0,0), vout, bous_sys.buoyancy_transport_eq)
        eq = linearize_eq(eq, ϵ)
    end
    
    surface_elevation_eq = let
        eq = evalat((0,0), vout, bous_sys.hydrostatic_balance_eq)
        eq = linearize_eq(eq, ϵ)
    end

    BoussinesqSystem(momentum_transport_eq, hydrostatic_balance_eq, continuity_eq, buoyancy_transport_eq, surface_elevation_eq)
end

function fourier_transform_eqs(eqs, k, l; u⃗̂, b̂, η̂, ŵ, p̂, ∫∇ᵀu⃗̂dz, du⃗, db, dη, dw, dp, ∫∇ᵀdu⃗dz)
    ϕ = exp.(im * [k * x + l * y for x = -nS:nS, y = -nS:nS])
    ϕcv = [exp(im * (k * 1//3 + l * -2//3)), exp(im * (k * -1//3 + l * -1//3))]
    ϕvc = [exp(im * (k * -1//3 + l * 2//3)), exp(im * (k * 1//3 + l * 1//3))]
    ϕcc = [1                                exp(im * (k * 0 + l * -1//sqrt3));
             exp(im * (k * 0 + l * 1//sqrt3)) 1]

    bous_sys_components = Dict{Symbol, BoussinesqSystem}()
    for component in [real, imag]
        subs_c = [
            merge(
                Dict(Symbolics.scalarize(du⃗) .=> component.([ϕ[x,y] * ϕcc[iHc,jHc] * u⃗̂[jTH,jHc] for jTH=1:2,x=1:size(ϕ,1),y=1:size(ϕ,2),jHc=1:2])),
                Dict(Symbolics.scalarize(db) .=> component.(ϕ * ϕcv[iHc] * b̂)),
                Dict(Symbolics.scalarize(dη) .=> component.(ϕ * ϕcv[iHc] * η̂)),
                Dict(Symbolics.scalarize(dp) .=> component.(ϕ * ϕcv[iHc] * p̂)),
                Dict(Symbolics.scalarize(dw) .=> component.(ϕ * ϕcv[iHc] * ŵ)),
                Dict(Symbolics.scalarize(∫∇ᵀdu⃗dz) .=> component.(ϕ * ϕcv[iHc] * ∫∇ᵀu⃗̂dz)),
            )
            for iHc = 1:2]
        momentum_transport_eqs = [
            [substitute(eq, subs_c[1]) for eq in bous_sys.momentum_transport_eqs[1:2]];
            [substitute(eq, subs_c[2]) for eq in bous_sys.momentum_transport_eqs[3:4]]
        ]
        
        subs_v = merge(
            Dict(Symbolics.scalarize(du⃗) .=> component.([ϕ[x,y] * ϕvc[jHc] * u⃗̂[jTH,jHc] for jTH=1:2,x=1:size(ϕ,1),y=1:size(ϕ,2),jHc=1:2])),
            Dict(Symbolics.scalarize(db) .=> component.(ϕ  * b̂)),
            Dict(Symbolics.scalarize(dη) .=> component.(ϕ  * η̂)),
            Dict(Symbolics.scalarize(dp) .=> component.(ϕ  * p̂)),
            Dict(Symbolics.scalarize(dw) .=> component.(ϕ  * ŵ)),
            Dict(Symbolics.scalarize(∫∇ᵀdu⃗dz) .=> component.(ϕ * ∫∇ᵀu⃗̂dz)),
        )
        bous_eqs = Dict{Symbol, Equation}()
        for eqn in [:hydrostatic_balance_eq, :continuity_eq, :buoyancy_transport_eq, :surface_elevation_eq]
            eq = get_property(bous_sys, eqn)
            bous_eqs[eqn] = substitute(eq, subs_v)
        end
        bous_sys_components[Symbol(component)] = BoussinesqSystem(momentum_transport_eqs, bous_eqs...)
    end
    NamedTuple(bous_sys_components)
end

function vertical_discretization()
    new_bous_sys_components = Dict{Symbol, ODESystem}()
    for (component, bous_sys) in pairs(bous_sys_components)

        ∫̲∇ᵀu⃗̂dz = let # integrate from -H to z
            ∇ᵀu⃗̂ = -symbolic_linear_solve(bous_sys.continuity_eq, ∂₃(ŵ))
            ∇ᵀu̲⃗̂ = [substitute(∇ᵀu⃗̂z, Dict([ û[iTH] => û̲[iV, iTH], v̂[iTH] => v̲̂[iV, iTH] for iTH=1:2])) for iV=1:Nz]
            dz * cumsum(∇ᵀu̲⃗̂)
        end
        ŵ̲ = -∫̲∇ᵀu⃗̂dz
        p̲̂ = let
            _b̂    = symbolic_linear_solve(bous_sys.hydrostatic_balance_eq, ∂₃(p̂))
            _b̲̂    = [substitute(_b̂, Dict([b̂ => b̲̂[iV]])) for iV=1:Nz]
            avb̲̂   = [[0.5 * _b̲̂[i] + 0.5 * _b̲̂[i+1] for i=1:Nz-1]; [0.0]]
            ∫̲b̂dz  =  -dz * reverse(cumsum(reverse(avb̲̂))) #integrate from 0 to z
            ∫̲b̂dz .- mean(∫̲b̂dz)
        end       

        for i=1:Nz
            subs = Dict([
                ŵ      => ŵ̲[iV],
                û[1]   => û̲[iV,1],
                û[2]   => û̲[iV,2],
                v̂[1]   => v̲̂[iV,1],
                v̂[2]   => v̲̂[iV,2],
                b̂      => b̲̂[iV],
                p̂      => p̲̂[iV],
                ∫∇ᵀu⃗̂dz => ∫̲∇ᵀu⃗̂dz[end]
            ])
            for eq in bous_sys.momentum_transport_eqs
                push!(new_eqs, substitute(eq, subs))
            end
            push!(new_eqs, substitute(bous_sys.buoyancy_eq, subs))
        end

        #create ODESystem
    end



        
        new_bous_sys_components[component] = 1
    end
    NamedTuple(new_bous_sys_components)
end
    

    
    for i=1:Nz
        subs = Dict([
            ŵ      => ŵ̲[iV],
            û[1]   => û̲[iV,1],
            û[2]   => û̲[iV,2],
            v̂[1]   => v̲̂[iV,1],
            v̂[2]   => v̲̂[iV,2],
            b̂      => b̲̂[iV],
            p̂      => p̲̂[iV],
            ∫∇ᵀu⃗̂dz => ∫̲∇ᵀu⃗̂dz[end]
        ])
        for eq in eqs
            push!(new_eqs, substitute(eq, subs))
        end
    end
end

function linearize_bous_sys(t, z; u⃗̄, b̄, η̄, w̄, p̄, du⃗, db, dη, dw, dp, f₀) #∫∇ᵀu⃗̄dz, ∫∇ᵀdu⃗dz
    vin  = @syms vin₁::Int vin₂::Int
    vout = @syms vout₁::Int vout₂::Int
    cin  = @syms cin₁::Int cin₂::Int cin₃::Int
    cout = @syms cout₁::Int cout₂::Int cout₃::Int
    @variables ϵ
    ∂ₜ = Differential(t)
    ∂₃ = Differential(z)

    u⃗ = [u⃗̄[iTH, cin[1], cin[2], cin[3]] + ϵ * du⃗[iTH, cin[1], cin[2], cin[3]] for iTH=1:2]
    b = b̄[vin[1], vin[2]] + ϵ * db[vin[1], vin[2]]
    η = η̄[vin[1], vin[2]] + ϵ * dη[vin[1], vin[2]]
    w = w̄[vin[1], vin[2]] + ϵ * dw[vin[1], vin[2]]
    p = p̄[vin[1], vin[2]] + ϵ * dp[vin[1], vin[2]]
    #∫∇ᵀu⃗dz = ∫∇ᵀu⃗̄dz[vin[1], vin[2]] + ϵ * ∫∇ᵀdu⃗dz[vin[1], vin[2]]
    
    eqs = [
        collect(∂ₜ.(evalat(cout, cin, u⃗)) .+ ∇ᵀcc(cout, cin, u⃗, u⃗) .+ ∂₃.(evalat(cout, cin, u⃗) .* av_cv(cout, vin, w)) .+ f₀ * [0 -1; 1 0] * evalat(cout, cin, u⃗) .~ -∇cv(cout, vin, p) .- ∇cv(cout, vin, η))...,
        ∇ᵀvc(vout, cin, u⃗) + ∂₃(evalat(vout, vin, w)) ~ 0,
        ∂₃(evalat(vout, vin, p)) - evalat(vout, vin, b) ~ 0, 
        ∂ₜ(evalat(vout, vin, b)) + ∇ᵀ(vout, cin, vin, u⃗, b) + ∂₃(evalat(vout, vin, w * b)) ~ 0,
        #∂ₜ(evalat(vout, vin, η)) + evalat(vout, vin, ∫∇ᵀu⃗dz) ~ 0
    ]
    eqs = [reduce(vcat, [evalat((0,0,iHc), cout, eqs[iTH]) for iHc=1:2,iTH=1:2]); [evalat((0,0), vout, eq) for eq in eqs[3:end]]]

    eqs = [Symbolics.taylor_coeff(simplify(expand_derivatives(eq.lhs)), ϵ, 1) ~ Symbolics.taylor_coeff(simplify(expand_derivatives(eq.rhs)), ϵ, 1) for eq in eqs]

    @show expand_derivatives.(eqs)
    expand_derivatives.(eqs)
end

function fourier_transform_eqs(eqs, k, l; u⃗̂, b̂, η̂, ŵ, p̂, du⃗, db, dη, dw, dp) # ∫∇ᵀu⃗̂dz, ∫∇ᵀdu⃗dz
    ϕ = exp.(im * [k * x + l * y for x = -nS:nS, y = -nS:nS])
    ϕcv = [exp(im * (k * 1//3 + l * -2//3)), exp(im * (k * -1//3 + l * -1//3))]
    ϕvc = [exp(im * (k * -1//3 + l * 2//3)), exp(im * (k * 1//3 + l * 1//3))]
    ϕcc = [1                                exp(im * (k * 0 + l * -1//sqrt3));
             exp(im * (k * 0 + l * 1//sqrt3)) 1]

    subs_c = [
        merge(
            Dict(Symbolics.scalarize(du⃗) .=> imag.([ϕ[x,y] * ϕcc[iHc,jHc] * u⃗̂[jTH,jHc] for jTH=1:2,x=1:size(ϕ,1),y=1:size(ϕ,2),jHc=1:2])),
            Dict(Symbolics.scalarize(db) .=> imag.(ϕ * ϕcv[iHc] * b̂)),
            Dict(Symbolics.scalarize(dη) .=> imag.(ϕ * ϕcv[iHc] * η̂)),
            Dict(Symbolics.scalarize(dp) .=> imag.(ϕ * ϕcv[iHc] * p̂)),
            Dict(Symbolics.scalarize(dw) .=> imag.(ϕ * ϕcv[iHc] * ŵ)),
            #Dict(Symbolics.scalarize(∫∇ᵀdu⃗dz) .=> imag.(ϕ * ϕcv[iHc] * ∫∇ᵀu⃗̂dz)),
        )
        for iHc = 1:2]
    subs_v = merge(
        Dict(Symbolics.scalarize(du⃗) .=> imag.([ϕ[x,y] * ϕvc[jHc] * u⃗̂[jTH,jHc] for jTH=1:2,x=1:size(ϕ,1),y=1:size(ϕ,2),jHc=1:2])),
        Dict(Symbolics.scalarize(db) .=> imag.(ϕ  * b̂)),
        Dict(Symbolics.scalarize(dη) .=> imag.(ϕ  * η̂)),
        Dict(Symbolics.scalarize(dp) .=> imag.(ϕ  * p̂)),
        Dict(Symbolics.scalarize(dw) .=> imag.(ϕ  * ŵ)),
        #Dict(Symbolics.scalarize(∫∇ᵀdu⃗dz) .=> imag.(ϕ * ∫∇ᵀu⃗̂dz)),
    )

    eqs = [[substitute(eq, subs_c[1]) for eq in eqs[1:2]]; [substitute(eq, subs_c[2]) for eq in eqs[3:4]]; [substitute(eq, subs_v) for eq in eqs[5:end]]]
end

# let ∇ᵀu⃗̂ := 
# let ∫∇ᵀu⃗̂dz := dz * cumsum(∇ᵀu⃗̂)
# loop over i=1:Nz and replace
# 
# ŵ(t,z)      => -∫̲∇ᵀu⃗̂dz[i]
# û(t,z)[iTH] => û̲(t)[i, iTH]
# v̂(t,z)[iTH] => v̲̂(t)[i, iTH]
# b̂(t,z)      => b̲̂(t)[i]
# η̂(t)        => η̂(t)
# ∫∇ᵀu⃗̂dz      => ∫̲∇ᵀu⃗̂dz[end]
function vertical_discretization()
    ∫̲∇ᵀu⃗̂dz = let # integrate from -H to z
        ∇ᵀu⃗̂ = -symbolic_linear_solve(continuity_eq, ∂₃(ŵ))
        ∇ᵀu̲⃗̂ = [substitute(∇ᵀu⃗̂z, Dict([ û[iTH] => û̲[iV, iTH], v̂[iTH] => v̲̂[iV, iTH] for iTH=1:2])) for iV=1:Nz]
        dz * cumsum(∇ᵀu̲⃗̂)
    end
    ŵ̲ = -∫̲∇ᵀu⃗̂dz
    p̲̂ = let
        _b̂    = symbolic_linear_solve(hydrstatic_balance_eq, ∂₃(p̂))
        _b̲̂    = [substitute(_b̂, Dict([b̂ => b̲̂[iV]])) for iV=1:Nz]
        avb̲̂   = [[0.5 * _b̲̂[i] + 0.5 * _b̲̂[i+1] for i=1:Nz-1]; [0.0]]
        ∫̲b̂dz  =  -dz * reverse(cumsum(reverse(avb̲̂))) #integrate from 0 to z
        ∫̲b̂dz .- mean(∫̲b̂dz)
    end
    
    for i=1:Nz
        subs = Dict([
            ŵ      => ŵ̲[iV],
            û[1]   => û̲[iV,1],
            û[2]   => û̲[iV,2],
            v̂[1]   => v̲̂[iV,1],
            v̂[2]   => v̲̂[iV,2],
            b̂      => b̲̂[iV],
            p̂      => p̲̂[iV],
            ∫∇ᵀu⃗̂dz => ∫̲∇ᵀu⃗̂dz[end]
        ])
        for eq in eqs
            push!(new_eqs, substitute(eq, subs))
        end
    end
end

function eady_problem(ū, v̄; f₀, M², N², H=4000)
    #@variables f₀ M2 N2
    @variables t z

    u⃗̄ = [ū; v̄] .* ones(1:1, -nS:nS, -nS:nS, 1:2) * (z-H/2) * M²
    db1 = M² * sqrt3//2
    b̄ = OffsetArray([x * db1 for x=-nS:nS, y=-nS:nS], -nS:nS, -nS:nS) * z * N²
    p̄ = OffsetArray([x * db1 for x=-nS:nS, y=-nS:nS], -nS:nS, -nS:nS) * 1//2 * (z-H)^2
    η̄ = zeros(-nS:nS, -nS:nS)
    w̄ = zeros(-nS:nS, -nS:nS)
    # assert geostrophic balance!
    
    @variables (du⃗(t,z))[1:2,-nS:nS,-nS:nS,1:2] (db(t,z))[-nS:nS,-nS:nS] (dη(t))[-nS:nS, -nS:nS]
    #@variables (∫∇ᵀdu⃗dz(t,z))[-nS:nS,-nS:nS]
    @variables (dw(t,z))[-nS:nS, -nS:nS] (dp(t,z))[-nS:nS, -nS:nS]

    eqs = linearize_bous_sys(t, z; u⃗̄, b̄, η̄, w̄, p̄, du⃗, db, dη, dw, dp, f₀)

    @variables k l
    @variables (û(t,z))[1:2] (v̂(t,z))[1:2] b̂(t,z) η̂(t)
    #@variables ∫∇ᵀu⃗̂dz(t,z)
    @variables ŵ(t,z) p̂(t,z)
    u⃗̂ = [û v̂]

    bb = le * Num[-1//2 -1; sqrt3//2 0]
    k̃ = bb[1,1] * k + bb[2,1] * l
    l̃ = bb[1,2] * k + bb[2,2] * l
    eqs = fourier_transform_eqs(eqs, k̃, l̃; u⃗̂, b̂, η̂, ŵ, p̂, du⃗, db, dη, dw, dp)
    eqs = substitute.(eqs, Ref(Dict([Av => _Av(le), Ac => _Ac(le), lê => _lê(le), sqrt3 => √3])))
    eqs = [expand_derivatives(eq.lhs) ~ expand_derivatives(eq.rhs) for eq in eqs]
    eqs = simplify.(eqs; expand=true)
end

function linearized_bous_sys()
    @variables f₀ M2
    @variables t z
    @variables (du(t,z))[1:2,-nS:nS,-nS:nS,1:2] (db(t,z))[-nS:nS,-nS:nS] (dη(t))[-nS:nS, -nS:nS]
    @variables (dw(t,z))[-nS:nS, -nS:nS] (dp(t,z))[-nS:nS, -nS:nS]
    vin  = @syms vin₁::Int vin₂::Int
    vout = @syms vout₁::Int vout₂::Int
    cin  = @syms cin₁::Int cin₂::Int cin₃::Int
    cout = @syms cout₁::Int cout₂::Int cout₃::Int
    ∂ₜ = Differential(t)
    ∂₃ = Differential(z)

    @variables ū v̄ ϵ
    db1 = M2 * sqrt3//2
    b̄ = OffsetArray([x * db1 for x=-nS:nS, y=-nS:nS], -nS:nS, -nS:nS)
    u⃗̄ = [1; 0] .* ones(1:1, -nS:nS, -nS:nS, 1:2)
    
    u⃗ = [u⃗̄[iTH, cin[1], cin[2], cin[3]] + ϵ * du[iTH, cin[1], cin[2], cin[3]] for iTH=1:2]
    b = b̄[vin[1], vin[2]] + ϵ * db[vin[1], vin[2]]
    η = ϵ * dη[vin[1], vin[2]]
    w = ϵ * dw[vin[1], vin[2]]
    p = ϵ * dp[vin[1], vin[2]]
    

    eqs = [
        collect(∂ₜ.(evalat(cout, cin, u⃗)) .+ ∇ᵀcc(cout, cin, u⃗, u⃗) .+ ∂₃.(evalat(cout, cin, u⃗) .* av_cv(cout, vin, w)) .+ f₀ * [0 -1; 1 0] * evalat(cout, cin, u⃗) .~ -∇cv(cout, vin, p) .- ∇cv(cout, vin, η))...,
        ∇ᵀvc(vout, cin, u⃗) + ∂₃(evalat(vout, vin, w)) ~ 0,
        ∂₃(evalat(vout, vin, p)) - evalat(vout, vin, b) ~ 0, 
        ∂ₜ(evalat(vout, vin, b)) + ∇ᵀ(vout, cin, vin, u⃗, b) + ∂₃(evalat(vout, vin, w * b)) ~ 0,
        #∂ₜ(evalat(vout, vin, η)) + ∇ᵀvc(vout, cin, ∫dz.(u⃗)) ~ 0
    ]
    eqs = [reduce(vcat, [evalat((0,0,iHc), cout, eqs[iTH]) for iHc=1:2,iTH=1:2]); [evalat((0,0), vout, eq) for eq in eqs[3:end]]]

    eqs = [Symbolics.taylor_coeff(simplify(expand_derivatives(eq.lhs)), ϵ, 1) ~ Symbolics.taylor_coeff(simplify(expand_derivatives(eq.rhs)), ϵ, 1) for eq in eqs]

    @variables k l
    @variables (û(t,z))[1:2]::Complex (v̂(t,z))[1:2]::Complex b̂(t,z)::Complex η̂(t)::Complex
    @variables ŵ(t,z)::Complex p̂(t,z)::Complex
    u⃗̂ = [û v̂]
    ϕ = exp.(im * [k * x + l * y for x = -nS:nS, y = -nS:nS])
    ϕcv = [exp(im * (k * 1//3 + l * -2//3)), exp(im * (k * -1//3 + l * -1//3))]
    ϕvc = [exp(im * (k * -1//3 + l * 2//3)), exp(im * (k * 1//3 + l * 1//3))]
    ϕcc = [1                                exp(im * (k * 0 + l * -1//sqrt3));
             exp(im * (k * 0 + l * 1//sqrt3)) 1]

    subs_c = [
        merge(
            Dict(Symbolics.scalarize(du) .=> [ϕ[x,y] * ϕcc[iHc,jHc] * u⃗̂[jTH,jHc] for jTH=1:2,x=1:size(ϕ,1),y=1:size(ϕ,2),jHc=1:2]),
            Dict(Symbolics.scalarize(db) .=> ϕ * ϕcv[iHc] * b̂),
            Dict(Symbolics.scalarize(dη) .=> ϕ * ϕcv[iHc] * η̂),
            Dict(Symbolics.scalarize(dp) .=> ϕ * ϕcv[iHc] * p̂),
            Dict(Symbolics.scalarize(dw) .=> ϕ * ϕcv[iHc] * ŵ),
        )
        for iHc = 1:2]
    subs_v = merge(
            Dict(Symbolics.scalarize(du) .=> [ϕ[x,y] * ϕvc[jHc] * u⃗̂[jTH,jHc] for jTH=1:2,x=1:size(ϕ,1),y=1:size(ϕ,2),jHc=1:2]),
            Dict(Symbolics.scalarize(db) .=> ϕ  * b̂),
            Dict(Symbolics.scalarize(dη) .=> ϕ  * η̂),
            Dict(Symbolics.scalarize(dp) .=> ϕ  * p̂),
            Dict(Symbolics.scalarize(dw) .=> ϕ  * ŵ),
        )

    eqs = [[substitute(eq, subs_c[1]) for eq in eqs[1:2]]; [substitute(eq, subs_c[2]) for eq in eqs[3:4]]; [substitute(eq, subs_v) for eq in eqs[5:end]]]
    
    
    # eqs = [simplify(substitute(eqs[1], subs_v); expand=true), substitute(eqs[2], subs_c[1]), substitute(eqs[3], subs_c[1]), substitute(eqs[4], subs_v), substitute(eqs[5], subs_v)]
    eqs
end
