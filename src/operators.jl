module GridOperators

using Symbolics
import Base: getindex

function getindex(A::AbstractArray, ii::Int, jj::Symbolics.Symbolic{<:Symbolics.Integer}, kk::Symbolics.Symbolic{<:Symbolics.Integer})
    Symbolics.Term{Symbolics.symeltype(A)}(Symbolics.getindex, [A, ii, jj, kk])
end

function getindex(A::AbstractArray, ii::Int, jj::Symbolics.Symbolic{<:Symbolics.Integer}, kk::Symbolics.Symbolic{<:Symbolics.Integer}, ll::Symbolics.Symbolic{<:Symbolics.Integer})
    Symbolics.Term{Symbolics.symeltype(A)}(Symbolics.getindex, [A, ii, jj, kk, ll])
end

function ∑(e, E::Vector{<:Tuple}, exp)
    @assert length(e) == length(E[1])
    sum([substitute(exp, merge([Dict(e[i] .=> e_[i]) for i in 1:length(e)]...)) for e_ in E])
end

function ∑(e, E::SymbolicUtils.BasicSymbolic{Vector{T}}, exp) where T <: Tuple
    if istree(E)
        @assert operation(E) == ifelse
        return ifelse(E.arguments[1], ∑(e, Symbolics.unwrap(E.arguments[2]), exp), ∑(e, Symbolics.unwrap(E.arguments[3]), exp)) # unwrap necessary because of a bug
    else
        throw(DomainError(E))
    end
end

function evalat(a::Vector, _xin, _xout)
    substitute.(a, Ref(Dict( _xin .=> _xout)))
end

function evalat(a, _xin, _xout)
    substitute(a, Dict(_xin .=> _xout))
end

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

module TriB
include("TriB/operators.jl")
end

module TriC
include("TriC/operators.jl")
end

module HexC
include("HexC/operators.jl")
end

end
