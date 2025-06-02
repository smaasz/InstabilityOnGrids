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

function dot(a,b)
    substitute(a[1] * b[1] + a[2] * b[2], Dict([sqrt3 => √3]))
end

n⃗e = ([ sqrt3//2;  1//2], [-sqrt3//2;  1//2], Num[ 0 ;-1//1])
t⃗e = ([-1//2;  sqrt3//2], [-1//2; -sqrt3//2], Num[ 1//1;  0]) # e⃗₃ × n⃗e

_lê(le) = 2//3 * sqrt3//2 * le

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
    ve = ∑(_v, VE(e), p // 2)
    [ 1//Ac * ∑((e, n⃗ec), EC(_c), ve * le * n⃗ec[iTH] ) for iTH = 1:2]
end

# implement ∇ec(e,c,F,n⃗ec) = [ ∑((c, d), CE(e, n⃗ec), d * F[jTH]) for jTH=1:2 ]?
# or d⃗ᵀ∇(e,c,F) = [ ∑((c,d), CE(e), d * F[jTH]) for jTH=1:2 ] and CE(e) are ordered according to d⃗ₑ = l⃗ₑ × k⃗, note n⃗ec n⃗ecᵀ∇ = d⃗ₑ d⃗ₑᵀ∇ is actually independent of c
# function ∇cc(_cout, _cin, F)
#     dF = [ ∑((_cin, d), CE(e, n⃗ec), d * F[jTH]) for jTH=1:2 ]
#     [∑((e, n⃗ec), EC(_cout), n⃗ec[iTH] * dF[jTH]) for iTH=1:2, jTH=1:2]
# end

function ∇cc(_cout, _cin, F⃗)
    dF = [ ∑((_cin, n⃗ec̃), CE(e), (n⃗ec' * n⃗ec) * F[jTH]) for jTH=1:2 ]
    [∑((e, n⃗ec), EC(_cout), n⃗ec[iTH] * dF[jTH]) for iTH=1:2, jTH=1:2]
end

# divergence operators

function ∇ᵀvc(_v, _c, F⃗)
    fe = ∑((_c, n⃗ec), CE(e), lê//2 * n⃗ev' * F⃗)
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
    F = ∑((_cin, n⃗eĉ), CE(e), le//2 * n⃗ec̃' * u⃗ad)
    Ru = [evalat(c̃, _cin, u⃗tr[iTH]) + le//2 * n⃗ec̃' * ∇(c̃, _cin, u⃗tr[iTH]) for iTH=1:2]
    Fu = ∑((c̃, n⃗ec̃), CE(e), n⃗ec' * n⃗ec̃ * 0.5*(F + abs(F)) * ( Ru[iTH] - evalat(_cout, _cin, u⃗tr[iTH]) )) 
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
    F = ∑((_c, n⃗ec), CE(e), lê//2 * n⃗ev' * u⃗)
    ∇beᶜ = 1//le * ∑((_vin, n⃗ev̂), VE(e), n⃗eṽ' * n⃗ev̂ * b)
    #∇beᵘ = n⃗eṽ' * ∇cv(cup(ṽ, n⃗eṽ), _vin, b)
    ∇beᶜ = 1//le * ∑((_vin, d), VEce(e,n⃗ev), d * b)
    ∇be  = 2//3 * ∇beᶜ + 1//3 * ∇beᵘ
    Rbe  = evalat(ṽ, _vin, b) + le//2 * ∇be
    1//Av * ∑((e, n⃗ev), EV(_vout), ∑((ṽ, n⃗eṽ), VE(e), n⃗ev' * n⃗eṽ * (F + (1-γ) * abs(F))/2 * Rbe))
end

const nS = 3

function bous_sys()
    @variables t z
    @variables (u(t,z))[1:2,-nS:nS,-nS:nS,-nS:nS] (b(t,z))[-nS:nS,-nS:nS] (η(t))[-nS:nS, -nS:nS]
    @variables (w(t,z))[-nS:nS, -nS:nS] (p(t,z))[-nS:nS, -nS:nS]
    vin  = @syms vin₁::Int vin₂::Int
    vout = @syms vout₁::Int vout₂::Int
    cin  = @syms cin₁::Int cin₂::Int cin₃::Int
    cout = @syms cout₁::Int cout₂::Int cout₃::Int
    ∂ₜ = Differential(t)
    ∂₃ = Differential(z)

    u⃗ = [u[iTH, c[1], c[2], c[3]] for iTH=1:2]
    b = b[v[1], v[2]]
    η = η[v[1], v[2]]
    

    eqs = [
        ∇ᵀvc(vout, cin, u⃗) + ∂₃(evalat(vout, vin, w)) ~ 0,
        ∂ₜ.(evalat(cout, cin, u⃗)) + ∇ᵀcc(cout, cin, u⃗, u⃗) + ∂₃(evalat(cout, cin, u⃗) .* av_cv(c̃, v, w)) + f₀ * [0 -1; 1 0] * evalat(cout, cin, u⃗) + ∇cv(cout, vin, p) + ∇cv(cout, vin, η) .~ 0,
        ∂₃(evalat(vout, vin, p)) - evalat(vout, vin, b) ~ 0, 
        ∂ₜ(evalat(vout, vin, b)) + ∇ᵀ(vout, cin, vin, u⃗, b) + ∂₃(evalat(vout, vin, w * b)) ~ 0,
        ∂ₜ(evalat(vout, vin, η)) + ∇ᵀvc(vout, cin, ∫dz.(u⃗)) ~ 0,
    ]
end
