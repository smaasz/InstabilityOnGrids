
function ∇ᵀce(_c, _e, u)
    1//Ac * ∑((_e, n⃗ec), EC(_c), le * n⃗ec' * n⃗e * u)
end

function curl_ve(_v, _e, u)
    1//Av * ∑((_e, n⃗ev), EV(_v), le * n⃗ev' * t⃗e * u)
end

function ∇ec(_e, _c, p)
    1//lê * ∑((_c, n⃗ec), CE(_e), n⃗ec' * n⃗e * p) 
end

function ∇⊥ev(_e, _v, ω)
    1//le * ∑((_v, n⃗ev), VE(_e), n⃗ev' * t⃗e * ω)
end

function Pce(_c, _e, u)
    1//Av * ∑(_e, E(_c), le * lê//2 * n⃗e * u)
end

function Pᵀec(_e, _c, F⃗)
    1//2 * ∑(_c, CE(_e), n⃗e' * F⃗)
end

function P̃ve(_v, _e, u)
    1//Ac * ∑(_e, EV(_v), le//2 * lê * n⃗e * u)
end

function P̃⁺ev(_e, _v, G⃗)
    1//2 * ∑(_v, VE(_e), t⃗e' * G⃗)
end

function ℳ(_eout, _ein, u)
    Pᵀec(_eout, c, Pce(c, _ein, u))
end

function ℳ(_eout, _ein, _c , u, b)
    Pᵀec(_eout, _c, b .* Pce(_c, _ein, u))
end

function ℳ̃(_eout, _ein, u)
    P̃⁺ev(_eout, v, P̃ve(v, _ein, u))
end

function ℳ̃(_eout, _ein, _v, u, ω)
    P̃⁺ev(_eout, _v, ω .* P̃ve(_v, _ein, u))
end

function ∇ᵀee(_eout, _ein, utr, uad)
    KE = Pce(c, _ein, utr)' * Pce(c, _ein, uad)
    ∇KE = ℳ(_eout, e, ∇ec(e, c, KE))
    ω = curl_ve(v, _ein, utr)
    ℳ̃(_eout, _ein, v, uad, ω) + ∇KE
end

function ∇ᵀcec(_cout, _e, _cin, u, b)
    ∇ᵀce(_cout, e, ℳ(e, _e, _cin, u, b))
end
