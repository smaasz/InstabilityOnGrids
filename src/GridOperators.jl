module GridOperators

using Symbolics
import OffsetArrays: OffsetArray
import Base: getindex

# function getindex(A::AbstractArray, ii::Int, jj::Symbolics.Symbolic{<:Symbolics.Integer}, kk::Symbolics.Symbolic{<:Symbolics.Integer})
#     Symbolics.Term{Symbolics.symeltype(A)}(Symbolics.getindex, [A, ii, jj, kk])
# end

function getindex(A::OffsetArray, ii::Int, jj::Symbolics.Symbolic{<:Symbolics.Integer}, kk::Symbolics.Symbolic{<:Symbolics.Integer}, ll::Symbolics.Symbolic{<:Symbolics.Integer})
    Symbolics.Term{Symbolics.symeltype(A)}(Symbolics.getindex, [A, ii, jj, kk, ll])
end

function ∑(e, E::Vector{<:Tuple}, exp)
    @assert length(e) == length(E[1])
    sum([substitute(exp, merge([Dict(e[1] .=> e_[1]), Dict(e[2] => e_[2])]...)) for e_ in E])
end

function ∑(e, E::SymbolicUtils.BasicSymbolic{Vector{T}}, exp) where T <: Tuple
    if istree(E)
        @assert operation(E) == ifelse
        return ifelse(E.arguments[1], ∑(e, Symbolics.unwrap(E.arguments[2]), exp), ∑(e, Symbolics.unwrap(E.arguments[3]), exp)) # unwrap necessary because of a bug
    else
        throw(DomainError(E))
    end
end

function evalat(_xout, _xin, a::Vector)
    substitute.(a, Ref(Dict( _xin .=> _xout)))
end

function evalat(_xout, _xin, a)
    substitute(a, Dict(_xin .=> _xout))
end

function evalat(_xout, _xin, a::Equation)
    substitute(a.lhs, Dict(_xin .=> _xout)) ~ substitute(a.rhs, Dict(_xin .=> _xout))
end

export evalat, ∑ 

module TriB
import ..GridOperators: evalat, ∑
import OffsetArrays: OffsetArray
using Symbolics
include("TriB/operators.jl")
export bous_sys, eady_problem
end

module TriC
using ..GridOperators
using Symbolics
include("TriC/operators.jl")
end

module HexC
using ..GridOperators
using Symbolics
include("HexC/operators.jl")
end

end
