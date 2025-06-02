
function ∇ᵀce(_c, _e, u)
    1//Ac * ∑((_e, n⃗ec), EC(_c), le * n⃗ec' * n⃗e * u)
end

function ∇ec(_e, _c, p)
    1//lê * ∑((_c, n⃗ec), CE(_e), n⃗ec' * n⃗e * p)
end

function ⊥(_eout, _ein, u)
    1//lê * ∑((_ein, w), w * le * u)
end

function ℳ̃(_eout, _ein, _v, u, ω)
    q̃ = av_ec(_ein, _v, ω)
    1//lê * ∑((_ein, w), ECP(_eout), w * le * u * (q̃ + evalat(q̃, _ein, _eout))//2)
end

function curl_ve(_v, _e, u)
    1//Av * ∑((_e, n⃗ev), EV(_v), le * n⃗ev' * t⃗e * u)
end

function ∇ᵀve(_v, _e, u)
    1//Av * ∑(c, CV(_v), Rcv * Ac * ∇ᵀce(c, _e, u))
end

function KEce(_c, _e, u)
    1//Av * ∑(_e, EC(_c), le * lê//2 * u * u)
end

function u⃗ᵀ∇_vector_invariant(_eout, _ein, u)
    ℳ̃(_eout, _ein, v, u, curl_ve(v, _ein, u)) + ∇ec(_eout, c, KE(c, _ein, u))
end

function ∇ᵀcec(_cout, _e, _cin, u, b)
    ∇ᵀce(_cout, _e, u * av_ec(_e, _cin, b))
end
