
# Immutable types as vectors

module MyVec

export IVector

immutable IVector{T<:Real}
    x::T
    y::T
    z::T
end

import Base.*
import Base.getindex
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

+{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x+b.x, a.y+b.y, a.z+b.z)
-{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x-b.x, a.y-b.y, a.z-b.z)
.*{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x*b.x, a.y*b.y, a.z*b.z)
./{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x/b.x, a.y/b.y, a.z/b.z)
.^{T<:Real}(a::IVector{T}, b::IVector{T}) = IVector{T}(a.x^b.x, a.y^b.y, a.z^b.z)

*{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x*b, a.y*b, a.z*b)
*{T<:Real}(a::Real, b::IVector{T}) = b*a
/{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x/b, a.y/b, a.z/b)
/{T<:Real}(a::Real, b::IVector{T}) = IVector{T}(a/(b.x), a/(b.y), a/(b.z))
+{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x+b, a.y+b, a.z+b)
+{T<:Real}(a::Real, b::IVector{T}) = b+a
-{T<:Real}(a::IVector{T}, b::Real) = IVector{T}(a.x-b, a.y-b, a.z-b)
-{T<:Real}(a::Real, b::IVector{T}) = IVector{T}(a-b.x, a-b.y, a-b.z)

end # module
