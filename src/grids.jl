module FourierSymbols

using TensorKit
using BlockTensorKit
using BlockTensorKit: ⊕
#using TensorOperations
using LinearAlgebra: Diagonal, I, diagm
using MultiplesOfPi

# grids on
## rectangular mesh
const rHc = ℂ^1
const rHe = ℂ^2
const rHv = ℂ^1

## triangular mesh
const tHc = ℂ^2
const tHe = ℂ^3
const tHv = ℂ^1

const tTH = ℂ^2
const e₁ = Tensor(ComplexF64[1,0], tTH)
const e₂ = Tensor(ComplexF64[0,1], tTH)

#const nu = [cos(π/6) -cos(π/6) 0.0; 0.5 0.5 -1.0]
const nu = [√3/2 -√3/2 0.0; 0.5 0.5 -1.0]
const mu = [0 -1; 1 0] * nu

Base.setindex!(d::Dict{CartesianIndex{2}, S}, v, k::Tuple{Int, Int}) where S = setindex!(d, v, CartesianIndex(k))


function bte(a, i)
    @assert 0 < i <= 3
    b = zeros(2,2)
    if i == 1
        b = [mu[:,1]*a  -mu[:,3]*a]
    elseif i == 2
        b = [mu[:,2]*a  mu[:,3]*a]
    elseif i == 3
        b = [mu[:,3]*a  mu[:,1]*a]
    end
    b
end

function btc(a, i)
    @assert 0 < i <= 2
    b = zeros(2,2)
    if i == 1
        b = [mu[:,1]*a -mu[:,3]*a]
    elseif i == 2
        b = [mu[:,1]*a -mu[:,3]*a]
    end
    b
end

btv(a) = [mu[:,1]*a -mu[:,3]*a]


bte(a,i,j) = bte(a,i)[:,j]
btc(a,i,j) = btc(a,i)[:,j]
btv(a,j)   = btv(a)[:,j]

function rtce(a,i,j)
    @assert 0 < i <= 2 && 0 < j <= 3
    r = zeros(2)
    if (i,j) == (1,1)
        r = -nu[:,1] * h(a)/3
    elseif (i,j) == (1,2)
        r = -nu[:,2] * h(a)/3
    elseif (i,j) == (1,3)
        r = -nu[:,3] * h(a)/3
    elseif (i,j) == (2,1)
        r = nu[:,1] * h(a)/3
    elseif (i,j) == (2,2)
        r = nu[:,2] * h(a)/3
    elseif (i,j) == (2,3)
        r = nu[:,3] * h(a)/3
    end
    r
end

rtec(a,j,i) = -rtce(a,i,j)

function rtcv(a,i)
    @assert 0 < i <= 2
    r = zeros(2)
    if i == 1
        r =  nu[:,1] * 2/3 * h(a)
    elseif i == 2
        r = -nu[:,1] * 2/3 * h(a)
    end
    r
end

rtvc(a,i) = -rtcv(a,i)

function rtve(a,i)
    @assert 0 < i <= 3
    r = zeros(2)
    if i == 1
        r = mu[:,1] * a/2
    elseif i == 2
        r = mu[:,2] * a/2
    elseif i == 3
        r = mu[:,3] * a/2
    end
    r
end

rtev(a,i) = -rtve(a,i)

function rtcc(a,i,j)
    @assert 0 < i <= 2 && 0 < j <= 2
    r = zeros(2)
    if (i,j) == (1,2)
        r = -nu[:,1] * 2/3 * h(a)
    elseif (i,j) == (2,1)
        r = nu[:,1] * 2/3 * h(a)
    end
    r
end
    
rot(θ) = TensorMap([cos(θ) -sin(θ); sin(θ) cos(θ)], tTH ← tTH)
h(a) = √3/2 * a

fouriertransform(stencil, k::Tensor; b₁, b₂, r=[0.;0.]) = fouriertransform(stencil, k.data; b₁, b₂, r)

function fouriertransform(stencil, k; b₁, b₂, r=[0.;0.])
    @assert isodd(size(stencil, 1)) && isodd(size(stencil,2))
    N₁ = size(stencil,1) ÷ 2
    N₂ = size(stencil,2) ÷ 2

    # fourier transform on grid + shift in space by r
    sum( stencil .* expm1.(-im * [k' * (b₁ * n₁ + b₂ * n₂ + r) for n₁ in -N₁:N₁, n₂ in -N₂:N₂]) ) + sum(stencil)
end

# Fourier Symbols on triangular B-grid
function vector_ec_up(k;a,θ=0)
    k = rot(θ) * k
    #t = zeros(ComplexF64,3,2,2)
    t = zeros(ComplexF64, tHe ← tTH ⊗ tHc)
    for iHe = 1:3
        for jTH = 1:2
            t[iHe,jTH,1] = fouriertransform([mu[jTH,iHe]], -k; b₁=btc(a,1,1), b₂=btc(a,1,2), r=rtce(a,1,iHe))
        end
    end
    t
end

function vector_ec_down(k;a,θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHe ← tTH ⊗ tHc)
    for iHe = 1:3
        for jTH = 1:2
            t[iHe,jTH,2] = fouriertransform([mu[jTH,iHe]], -k; b₁=btc(a,2,1), b₂=btc(a,2,2), r=rtce(a,2,iHe))
        end
    end
    t
end

vector_ec_av(k;a,θ=0.0) = 0.5 * vector_ec_up(k;a,θ) + 0.5 * vector_ec_down(k;a,θ)

# divergence operators

function div_ce(k;a,θ=0.0)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHc ← tHe)
    for iHe = 1:3
        t[1,iHe] = fouriertransform([1], -k; b₁=bte(a,iHe,1), b₂=bte(a,iHe,2), r=rtec(a,iHe,1)) / (h(a) / 2)
    end
    t[2,:] .= -conj.(t[1,:])
    t
end

function div_ve(k;a,θ=0)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHv ← tHe)
    for iHe = 1:3
        t[1,iHe] = fouriertransform([-2/3; 2/3; 0], -k; b₁=bte(a,iHe,1), b₂=bte(a,iHe,2), r=rtev(a,iHe)) / a
    end
    t
end
 
# function div_vc(k;a,θ=0.)
#     k = rot(θ) * k
#     div_ve_ = div_ve(k;a)
#     vector_ec_av_ = vector_ec_av(k;a)
#     @tensor t[iHv;jTH jHc] := div_ve_[iHv,jHe] * vector_ec_av_[jHe,jTH,jHc]
#     t
# end

function div_vc(k;a,θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHv ← tTH ⊗ tHc)
    for jTH = 1:2
        t[1,jTH,1] = fouriertransform([0. 0. nu[jTH,3]; 0. nu[jTH,1] nu[jTH,2]; 0. 0. 0.], -k; b₁=btc(a,1,1), b₂=btc(a,1,2), r=rtcv(a,1)) / (2 * h(a))
    end
    t[1,:,2] = -conj.(t[1,:,1])
    t
end

function div_cv(k;a,θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHc ← tTH ⊗ tHv)
    for jTH = 1:2
        t[1,jTH,1] = fouriertransform([0. 0. 0.; -nu[jTH,2] -nu[jTH,1] 0.; -nu[jTH,3] 0. 0.], -k; b₁=btv(a,1), b₂=btv(a,2), r=rtvc(a,1)) / h(a)
    end
    t[2,:,1] = -conj.(t[1,:,1])
    t
end

function div_cc(k;a,θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHc ← tTH ⊗ tHc)
    for jTH = 1:2
        t[1,jTH,2] = fouriertransform([0. 0. nu[jTH,3]; 0. nu[jTH,1] nu[jTH,2]; 0. 0. 0.], -k; b₁=btc(a,2,1), b₂=btc(a,2,2), r=rtcc(a,2,1)) / h(a)
    end
    t[2,:,1] = -conj.(t[1,:,2])
    t
end

# gradient operators

function grad_cv(k;a,θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tTH ⊗ tHc ← tHv)
    for jTH = 1:2
        t[jTH,1,1] = fouriertransform([0 0 0; -nu[jTH,2] -nu[jTH,1] 0; -nu[jTH,3] 0 0] , -k; b₁=btv(a,1), b₂=btv(a,2), r=rtvc(a,1)) / h(a)
    end
    t[:,2,1] = -conj.(t[:,1,1])
    t
end

function grad_cc(k;a,θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tTH ⊗ tHc ← tHc)
    for jTH = 1:2
        t[jTH,1,2] = fouriertransform([0. 0. nu[jTH,3]; 0. nu[jTH,1] nu[jTH,2]; 0. 0. 0.], -k; b₁=btc(a,2,1), b₂=btc(a,2,2), r=rtcc(a,2,1)) ./ h(a)
    end
    t[:,2,1] = -conj.(t[:,1,2])
    t
end

# curl operators

function curl_vc(k;a,θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHv ← tTH ⊗ tHc)
    for jTH = 1:2
        t[1,jTH,1] = fouriertransform([0. 0. mu[jTH,3]; 0. mu[jTH,1] mu[jTH,2]; 0. 0. 0.], -k; b₁=btc(a,1,1), b₂=btc(a,1,2), r=rtcv(a,1)) ./ (2 * h(a))
    end
    t[1,:,2] = -conj.(t[1,:,1])
    t
end

# Laplacian

function laplacian_cc(k; a, θ=0.0)
    k = rot(θ) * k
    t = zeros(Complex{BigFloat}, tHc ← tHc)
    t[1,1] = -12. / a^2
    t[2,2] = -12. / a^2
    t[1,2] = 4. * fouriertransform(BigFloat[0. 0. 1.0; 0. 1.0 1.0; 0. 0. 0.] ./ a^2, -k; b₁=btc(a,2,1), b₂=btc(a,2,2), r=rtcc(a,2,1))
    t[2,1] = conj(t[1,2])
    t
end

# linearized horizontal momentum advection

function hma_advective_form_cc(k, U::Tensor{T, ComplexSpace, 1, Vector{T}} where T <: Number; a, θ=0., θU=0.)
    k = rot(θ) * k
    d = rot(θU) * e₁
    Nz = dim(U)
    V = codomain(U)
    grad_cc_ = grad_cc(k;a)
    pw_mult = pointwise_mult(V)
    @tensor begin
        ū_mult[iV;jV jTH] := pw_mult[iV,jV,kV] * U[kV] * d'[jTH]
        t[iV,iTH,iHc;jV,jTH,jHc] := ū_mult[iV,kV,kTH] * grad_cc_[kTH,iHc,jHc] * id(tTH)[iTH,jTH] * id(V)[kV,jV]
    end
    t
end

function hma_vec_inv_form_cc(k, U::Tensor{T, ComplexSpace, 1, Vector{T}} where T <: Number; a, θ=0., θU=0.)
    k = rot(θ) * k
    d = rot(θU) * e₁
    Nz = dim(U)
    V = codomain(U)
    av_cv_ = average_cv(k; a)
    curl_vc_ = curl_vc(k; a)
    av_vc_ = average_vc(k; a)
    grad_cv_ = grad_cv(k; a)
    pw_mult = pointwise_mult(V)
    @tensor begin
        ū_mult1[iV,iTH;jV] := pw_mult[iV,jV,kV] * U[kV] * d[iTH]
        ū_mult2[iV;jV jTH] := pw_mult[iV,jV,kV] * U[kV] * d'[jTH]
        t1[iV,iTH,iHc;jV,jTH,jHc] := rot(Pi//2)[iTH,kTH] * ū_mult1[iV,kTH,kV] * av_cv_[iHc, kHv] * curl_vc_[kHv,jTH,jHc] * id(V)[kV,jV]
        t2[iV,iTH,iHc;jV,jTH,jHc] := ū_mult2[iV,kV,jTH] * grad_cv_[iTH,iHc,kHv] * av_vc_[kHv,jHc] * id(V)[kV,jV]
    end
    t1 + t2
end

#function hma_

# function hma_advective_form_cc(k, U::Number; a, θ=0., θU=0.)
#     ū = rot(θU)[:,1] * U
#     t = zeros(ComplexF64, 2, 2, 2, 2)
#     grad_cc_ = grad_cc(k;a,θ)
#     @tensor t[iTH,iHc,jTH,jHc] = ū[kTH] * grad_cc_[kTH,iHC,jHc] * I(2)[iTH,jTH]
#     t
# end

# horizontal average operators
function average_vc(k; a, θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHv ← tHc)
    t[1,1] = fouriertransform([0. 0. 1/6; 0. 1/6 1/6; 0. 0. 0.], -k; b₁=btc(a,1,1), b₂=btc(a,1,2), r=rtcv(a,1))
    t[1,2] = conj(t[1,1])
    t
end

function average_cv(k; a, θ=0.)
    k = rot(θ) * k
    t = zeros(ComplexF64, tHc ← tHv)
    t[1,1] = fouriertransform([0. 0. 0.; 0. 1/3 1/3; 0. 1/3 0.],-k; b₁=btv(a,1), b₂=btv(a,2), r=rtvc(a,1))
    t[2,1] = conj(t[1,1])
    t
end

# vertical average
average_vertical_pd(Nz) = TensorMap(diagm(0 => 0.5 * ones(ComplexF64,Nz), -1=> 0.5 * ones(ComplexF64, Nz-1)), ℂ^Nz ← ℂ^Nz)
#average_vertical_pd(Nz) = TensorMap(1/2 * [zeros(1,Nz-1); diagm(ones(Nz-1))] + 1/2 * [diagm(ones(Nz-1)); zeros(1,Nz-1)], ℂ^Nz ← ℂ^(Nz-1))
average_vertical_dp(Nz) = TensorMap(1/2 * [zeros(Nz-1,1) diagm(ones(Nz-1))] + 1/2 * [diagm(ones(Nz-1)) zeros(Nz-1,1)], ℂ^(Nz-1) ← ℂ^Nz)
derivative_vertical_dp(Nz, dz) = TensorMap(([zeros(Nz-1,1) diagm(ones(Nz-1))] - [diagm(ones(Nz-1)) zeros(Nz-1,1)]) / dz, ℂ^(Nz-1) ← ℂ^Nz)
derivative_vertical_pd(Nz, dz) = TensorMap(([diagm(ones(Nz-1)); zeros(1,Nz-1)] - [zeros(1,Nz-1); diagm(ones(Nz-1))]) / dz, ℂ^Nz ← ℂ^(Nz-1))

function vma(k, dū_dz::Tensor{T, ComplexSpace, 1, Vector{T}} where T <: Number; a, θ=0.0, θU=0.0)
    k = rot(θ) * k
    d = rot(θU) * e₁
    V = codomain(dū_dz)
    Nz = dim(dū_dz)
    average_vertical_pd_ = average_vertical_pd(Nz)
    average_cv_ = average_cv(k;a)
    pw_mult = pointwise_mult(V)
    @tensor begin
        dū_dz_mult[iV,iTH;jV] := pw_mult[iV,jV,kV] * dū_dz[kV] * d[iTH]
        t[iV,iTH,iHc;jṼ,jHv] := dū_dz_mult[iV,iTH,kV] * average_cv_[iHc,jHv] * average_vertical_pd_[kV,jṼ]
    end
    t
end

# function vma(k, dū_dz::Number; a, θ=0.0, θU=0.0)
#     d = rot(θU) * e₁
#     average_vertical_pd_ = average_vertical_pd(Nz)
#     average_cv_ = average_cv(k;a,θ)
#     pw_mult = pointwise_mult(Nz)
#     @tensor begin
#         dū_dz_mult[iV,jV, jTH] := pw_mult[iV,jV,kV] * ones(Nz)[kV] * dū_dz * d[jTH]
#         t[iV,iHc,jṼ,jTH,jHv] := dū_dz_mult[iV,kV,jTH] * average_cv_[IHc,jHv] * average_vertical_pd[kV,jṼ]
#     end
#     t
# end

function hsa(k, M2, U::Tensor{T, ComplexSpace, 1, Vector{T}} where T <: Number; a, θ=0.0, θU=0.0)
    k = rot(θ) * k
    d = rot(θU) * e₁
    @tensor ū[iV,iTH] := U[iV] * d[iTH]
    ∇b̄ = convert(Array, M2 * rot(θU) * e₂)
    V = codomain(U)
    Nz = dim(U)

    t_ = zeros(ComplexF64, tHe ← tHv)
    t_[3,1] = fouriertransform([0. 0. 0. 0. 0.; 0. 0. 0. 0. 0.; -1/12 7/12 7/12 -1/12 0.; 0. 0. 0. 0. 0.; 0. 0. 0. 0. 0.], -k; b₁=btv(a,1), b₂=btv(a,2), r=rtve(a,3))
    t_[1,1] = fouriertransform([0. 0. 0. 0. 0.; 0. 0. 0. 0. 0.; -1/12 7/12 7/12 -1/12 0.; 0. 0. 0. 0. 0.; 0. 0. 0. 0. 0.], -rot(-2//3 * Pi)*k; b₁=btv(a,1), b₂=btv(a,2), r=rtve(a,3))
    t_[2,1] = fouriertransform([0. 0. 0. 0. 0.; 0. 0. 0. 0. 0.; -1/12 7/12 7/12 -1/12 0.; 0. 0. 0. 0. 0.; 0. 0. 0. 0. 0.], -rot(2//3 * Pi)*k; b₁=btv(a,1), b₂=btv(a,2), r=rtve(a,3))
    
    div_ve_ = div_ve(k; a)
    pw_mult = pointwise_mult(V, tHe)
    mu_ = Tensor(mu, tTH' ⊗ tHe)
    @tensor begin
        Q_mult[iV,iHe;jV,jHe] := pw_mult[iV,iHe,jV,jHe,kV,kHe] * ū[kV,kTH] * mu_[kTH,kHe]
        tbb[iV,iHv;jV,jHv] := div_ve_[iHv, iHe] * Q_mult[iV,iHe,kV,kHe] * t_[kHe,jHv] * id(V)[kV,jV]
    end

    tbu_ = zeros(ComplexF64, tHv ← tTH ⊗ tHc)
    v1 = mu[:,2] * mu[:,2]' + mu[:,3] * mu[:,3]'
    v2 = mu[:,3] * mu[:,3]' + mu[:,1] * mu[:,1]'
    v3 = mu[:,1] * mu[:,1]' + mu[:,2] * mu[:,2]'
    for jTH = 1:2
        tbu_[1, jTH, 1] = fouriertransform([0. 0. 1/6 * ∇b̄' * v3[:,jTH]; 0. 1/6 * ∇b̄' * v1[:,jTH] 1/6 * ∇b̄' * v2[:,jTH]; 0. 0. 0.], -k; b₁=btc(a,1,1), b₂=btc(a,1,2), r=rtcv(a,1))
    end
    tbu_[1, :, 2] = conj.(tbu_[1, :, 1])
    @tensor tbu[iV, iHv; jV, jTH, jHc] := tbu_[iHv, jTH, jHc] * id(V)[iV, jV]
    
    # vector_ec_av_ = vector_ec_av(k;a)
    # pw_mult = pointwise_mult(3)
    # @tensor begin
    #     mult[iHe, jHe] := pw_mult[iHe, jHe, kHe] * ∇b̄[kTH] * mu[kTH,kHe]
    #     tbu[iV,iHv,jV,jTH,jHc] := div_ve_[iHv,iHe] * mult[iHe,jHe] * vector_ec_av_[jHe, jTH, jHc] * I(Nz)[iV,jV]
    # end
    (tbu, tbb)
end



function pointwise_mult(s...)
    ps = ⊗(s...)
    t = zeros(ComplexF64, ps ← ps ⊗ ps)
    for i in Iterators.product([1:n for n in dims(ps)]...)
        t[i...,i...,i...] = 1.0
    end
    t
end

function build_system(k; a, Ri, N=1.0e-3, f₀=-1.0e-4, H=4.0e3, Nz=32, θ=0.0, θU=0.0)
    M2 = √(N^2 * f₀^2 / Ri)
    dz = H / Nz
    V = ℂ^Nz
    Ṽ = ℂ^(Nz-1)
    U = Tensor(-M2 * dz / f₀ * (collect(1:Nz) .- (Nz + 1)/2 * ones(Nz)), V)
    dū_dz = -M2 / f₀ * ones(V)
    db̄_dz = N^2 * ones(V)
    db̄_dy = M2

    @tensor coriolis_cc[iV, iTH, iHc; jV, jTH, jHc] := f₀ * rot(Pi//2)[iTH, jTH] * id(V)[iV, jV] * id(tHc)[iHc, jHc]

    iu, iw, ip, iϕ, ib = 1:5
    blockmaps = Dict{CartesianIndex{2}, TensorMap{ComplexF64, ComplexSpace,1,1, Vector{ComplexF64}}}()
    # du/dt    
    Auu = -hma_advective_form_cc(k, U; a, θ, θU) - coriolis_cc
    Auw = -vma(k, dū_dz; a, θ, θU)
    @tensor Aup[iV, iTH, iHc; jV, jHv] := -grad_cv(k; a, θ)[iTH, iHc, jHv] * id(V)[iV,jV]
    @tensor Auϕ[iV, iTH, iHc; jHv] := -grad_cv(k; a, θ)[iTH, iHc, jHv] * ones(V)[iV]
    # 0.0 * dw/dt
    @tensor Awp[iṼ, iHv; jV, jHv] := -derivative_vertical_dp(Nz, dz)[iṼ, jV] * id(tHv)[iHv, jHv]
    @tensor Awb[iṼ, iHv; jV, jHv] := average_vertical_dp(Nz)[iṼ, jV] * id(tHv)[iHv, jHv]
    # 0.0 * dp/dt
    @tensor Apu[iV, iHv; jV, jTH, jHc] := div_vc(k; a, θ)[iHv, jTH, jHc] * id(V)[iV, jV]
    @tensor Apw[iV, iHv; jṼ, jHv] := derivative_vertical_pd(Nz, dz)[iV, jṼ] * id(tHv)[iHv, jHv]

    # 0.0 * dϕ/dt
    @tensor Aϕp[iHv; jV jHv] := id(tHv)[iHv, jHv] * ones(V)'[jV]

    # db/dt
    Abu, Abb = hsa(k, M2, U; a, θ, θU)
    Abu = - Abu
    Abb = -Abb
    @tensor Abw[iV, iHv; jṼ, jHv] := -N^2 * average_vertical_pd(Nz)[iV, jṼ] * id(tHv)[iHv, jHv]

    uspace = V ⊗ tTH ⊗ tHc
    wspace = Ṽ ⊗ tHv
    pspace = V ⊗ tHv
    ϕspace = tHv
    bspace = V ⊗ tHv
    
    eu = unitary(fuse(uspace)←uspace)
    ew = unitary(fuse(wspace)←wspace)
    ep = unitary(fuse(pspace)←pspace)
    eϕ = unitary(fuse(ϕspace)←ϕspace)
    eb = unitary(fuse(bspace)←bspace)

    blockmaps = Dict([
        CartesianIndex(iu, iu) => eu * Auu * eu',
        CartesianIndex(iu, iw) => eu * Auw * ew',
        CartesianIndex(iu, ip) => eu * Aup * ep',
        CartesianIndex(iu, iϕ) => eu * Auϕ * eϕ',
        CartesianIndex(iw, ip) => ew * Awp * ep',
        CartesianIndex(iw, ib) => ew * Awb * eb',
        CartesianIndex(ip, iu) => ep * Apu * eu',
        CartesianIndex(ip, iw) => ep * Apw * ew',
        CartesianIndex(iϕ, ip) => eϕ * Aϕp * ep',
        CartesianIndex(ib, iu) => eb * Abu * eu',
        CartesianIndex(ib, ib) => eb * Abb * eb',
        CartesianIndex(ib, iw) => eb * Abw * ew',
    ])
    
    statespace = fuse(uspace) ⊕ fuse(wspace) ⊕ fuse(pspace) ⊕ fuse(ϕspace) ⊕ fuse(bspace)
    A = SparseBlockTensorMap{TensorMap{ComplexF64, ComplexSpace,1,1, Vector{ComplexF64}}}(blockmaps, statespace ← statespace)

    massblockmaps = Dict([
        CartesianIndex(iu,iu) => id(ComplexF64, fuse(uspace)),
        CartesianIndex(ib,ib) => id(ComplexF64, fuse(bspace)),
    ])
    B = SparseBlockTensorMap{TensorMap{ComplexF64, ComplexSpace,1,1, Vector{ComplexF64}}}(massblockmaps, statespace ← statespace)
    
    (A,B)
end


function build_system_free_surface(k; a, Ri, N=1.0e-3, f₀=-1.0e-4, g=9.81, H=4.0e3, Nz=32, θ=0.0, θU=0.0, u0)
    k = rot(θ) * k
    θ = 0.0
    M2 = √(N^2 * f₀^2 / Ri)
    dz = H / Nz
    V = ℂ^Nz
    Ṽ = ℂ^Nz
    U = Tensor(-M2 * dz / f₀ * (collect(1:Nz) .- (Nz + 1)/2 * ones(Nz)) .+ u0, V)
    dū_dz = -M2 / f₀ * ones(V)
    db̄_dz = N^2 * ones(V)
    db̄_dy = M2

    eNz = Tensor(id(V)[:,Nz], V)
    wu = -dz * TensorMap(ComplexF64[i <= j ? 1. : 0 for j = 1:Nz, i = 1:Nz], Ṽ ← V) ⊗ div_vc(k;a,θ)
    pb = let
        Δ = diagm(0 => -ones(ComplexF64, Nz), 1 => ones(ComplexF64, Nz-1))
        Δ[end,:] .= 1.0
        av = 0.5 * diagm(0 => ones(ComplexF64, Nz), 1 => ones(ComplexF64, Nz-1))
        av[end,:] .= 0.0
        TensorMap(dz * inv(Δ) * av, V ← V) ⊗ id(tHv)
    end

    iu, ib, iη = 1:3
    blockmaps = Dict{CartesianIndex{2}, TensorMap{ComplexF64, ComplexSpace,1,1, Vector{ComplexF64}}}()
    # du/dt
    Auu = let
        @tensor coriolis_cc[iV, iTH, iHc; jV, jTH, jHc] := f₀ * rot(Pi//2)[iTH, jTH] * id(V)[iV, jV] * id(tHc)[iHc, jHc]
        Auw = -vma(k, dū_dz; a, θ, θU)
        -hma_advective_form_cc(k, U; a, θ, θU) - coriolis_cc + Auw * wu
    end
    Aub = let
        @tensor Aup[iV, iTH, iHc; jV, jHv] := -grad_cv(k;a,θ)[iTH, iHc, jHv] * id(V)[iV,jV]
        Aup * pb
    end
    Auη = -g * ones(V) ⊗ grad_cv(k;a,θ)

    # db/dt
    Abu, Abb = hsa(k, M2, U; a, θ, θU)
    Abu = let
        @tensor Abw[iV, iHv; jṼ, jHv] := -N^2 * average_vertical_pd(Nz)[iV, jṼ] * id(tHv)[iHv, jHv]
        -Abu + Abw * wu
    end
    Abb = -Abb

    # dη/dt
    Aηu = (eNz' ⊗ id(tHv)) * wu
    Aηη = (eNz' ⊗ id(tHv)) * Abb * (eNz ⊗ id(tHv))

    uspace = V ⊗ tTH ⊗ tHc
    bspace = V ⊗ tHv
    ηspace = tHv    
    
    eu = unitary(fuse(uspace)←uspace)
    eb = unitary(fuse(bspace)←bspace)
    eη = unitary(fuse(ηspace)←ηspace)

    blockmaps = Dict([
        CartesianIndex(iu, iu) => eu * Auu * eu',
        CartesianIndex(iu, ib) => eu * Aub * eb',
        CartesianIndex(iu, iη) => eu * Auη * eη',
        CartesianIndex(ib, iu) => eb * Abu * eu',
        CartesianIndex(ib, ib) => eb * Abb * eb',
        CartesianIndex(iη, iu) => eη * Aηu * eu',
        #CartesianIndex(iη, iη) => eη * Aηη * eη',
    ])
    statespace = fuse(uspace) ⊕ fuse(bspace) ⊕ fuse(ηspace)
    SparseBlockTensorMap{TensorMap{ComplexF64, ComplexSpace,1,1, Vector{ComplexF64}}}(blockmaps, statespace ← statespace)
end

export div_vc, div_cc, grad_cv, grad_cc, curl_vc, hma_advective_form_cc, average_cv, vma, pointwise_mult

end
