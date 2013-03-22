
# Immutable types as vectors

module MyVec

export IVector

immutable IVector{T<:Real}
    x::T
    y::T
    z::T
end

importall Base
function getindex{T}(a::IVector{T}, i::Integer)
    if i == 1
        return a.x
    elseif i == 2
        return a.y
    elseif i == 3
        return a.z
    else
        throw("Invalid index for IVector: $i")
    end
end

getindex{T<:Real}(A::Array{IVector{T}}, i::Integer, j::Integer) = A[i][j]

convert{T<:Real}(::Type{IVector{T}}, a::Array{T,1})= size(a,1) != 3 ? throw("Conversion error") : IVector(a[1], a[2], a[3])

zero{T<:Real}(::Type{IVector{T}}) = IVector{T}(zero(T), zero(T), zero(T))

norm{T<:Real}(V::IVector{T}) = sqrt(V.x^2 + V.y^2 + V.z^2)

+{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x+b.x, a.y+b.y, a.z+b.z)
-{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x-b.x, a.y-b.y, a.z-b.z)
.*{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x*b.x, a.y*b.y, a.z*b.z)
./{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x/b.x, a.y/b.y, a.z/b.z)
.^{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x^b.x, a.y^b.y, a.z^b.z)

.*{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x*b, a.y*b, a.z*b)
.*{T<:Real}(a::Real, b::IVector{T}) = b*a
*{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x*b, a.y*b, a.z*b)
*{T<:Real}(a::Real, b::IVector{T}) = b*a
./{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x/b, a.y/b, a.z/b)
./{T<:Real}(a::Real, b::IVector{T}) = IVector{T}(a/(b.x), a/(b.y), a/(b.z))
+{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x+b, a.y+b, a.z+b)
+{T<:Real}(a::Real, b::IVector{T}) = b+a
-{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x-b, a.y-b, a.z-b)
-{T<:Real}(a::Real, b::IVector{T}) = IVector{T}(a-b.x, a-b.y, a-b.z)

*{T<:Real}(A::Array{IVector{T}}, x::Real) = reshape([A[i]*x for i=1:length(A)], size(A))
*{T<:Real}(x::Real, A::Array{IVector{T}}) = A*x
+{T<:Real}(A::Array{IVector{T}}, x::Real) = reshape([A[i]+x for i=1:length(A)], size(A))
+{T<:Real}(x::Real, A::Array{IVector{T}}) = A*x


end # module
