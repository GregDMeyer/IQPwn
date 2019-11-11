"""
Data structures and functions for operating on vectors over GF2.
We implement these ourselves because the BitArray type from Julia does not operate over this field.
E.g., summing vectors is not mod 2. However this type borrows ideas from the BitArray type.
"""

# several of the functions here are derived from the BitArray definition file,
# from Julia's source: https://github.com/JuliaLang/julia/blob/v1.2.0/base/bitarray.jl
# which is licensed under the MIT license.

"""
    GF2Array{N} <: AbstractArray{Bool, N}
"""
mutable struct GF2Array{N} <: AbstractArray{Bool, N}
    # unlike in BitArray, we keep an underlying array of UInt64. The main difference is
    # that we're always aligned so that each column of our GF2Vector starts at the beginning
    # of a UInt64

    chunks::Array{UInt64, N}
    len::Int
    dims::NTuple{N,Int}

    function GF2Array{N}(::UndefInitializer, dims::Vararg{Int,N}) where N
        n = 1
        i = 1
        for d in dims
            d >= 0 || throw(ArgumentError("dimension size must be ≥ 0, got $d for dimension $i"))
            n *= d
            i += 1
        end
        nc = num_bit_chunks(dims[1])  # number of chunks for one column
        chunks = fill(zero(UInt64), (nc, dims[2:end]...))
        b = new(chunks, n)
        N != 1 && (b.dims = dims)
        return b
    end
end

GF2Array(::UndefInitializer, dims::Integer...) = GF2Array(undef, map(Int,dims))
GF2Array{N}(::UndefInitializer, dims::Integer...) where {N} = GF2Array{N}(undef, map(Int,dims))
GF2Array(::UndefInitializer, dims::NTuple{N,Integer}) where {N} = GF2Array{N}(undef, map(Int, dims)...)
GF2Array{N}(::UndefInitializer, dims::NTuple{N,Integer}) where {N} = GF2Array{N}(undef, map(Int, dims)...)

const GF2Vector = GF2Array{1}
const GF2Matrix = GF2Array{2}

@inline _div64(l) = l >> 6
num_bit_chunks(n::Int) = _div64(n+63)

import Base: size, length, getindex, setindex!, IndexStyle
length(x::GF2Array) = x.len
size(x::GF2Vector) = (x.len,)
size(x::GF2Array) = x.dims

# all we need for our purposes is:
#  - efficient inner product (with special case of mod 2)
#  - efficient vector sum

"""
    dot(a::GF2Vector, b::GF2Vector)

Inner product modulo 2
"""
function dot(a::GF2Vector, b::GF2Vector)
    length(a) != length(b) && throw(ArgumentError("vectors have different lengths $(length(a)) and $(length(b))"))
    t = zero(UInt64)
    @inbounds for i in 1:length(a.chunks)
        t ⊻= a.chunks[i] & b.chunks[i]
    end
    return isodd(count_ones(t))
end

"""
    dotcol(a, b, i)

Dot of column i of matrix b to a
"""
function dotcol(a::GF2Vector, b::GF2Matrix, i::Int)
    length(a) != size(b)[1] && throw(ArgumentError("dimension mismatch: vector $(length(A)) and matrix $(size(B))"))
    @boundscheck i > size(b)[2] && throw(BoundsError("column index $i greater than number of columns $(size(b)[2])"))
    unsafe_dotcol(a, b, i)
end

"""
    unsafe_dotcol(a, b, i)
"""
function unsafe_dotcol(a::GF2Vector, b::GF2Matrix, i::Int)
    chunkoffset = (i-1)*size(b.chunks)[1]
    t = zero(UInt64)
    for i in 1:length(a.chunks)
        t ⊻= a.chunks[i] & b.chunks[i + chunkoffset]
    end
    return isodd(count_ones(t))
end

"""
    add!(a, b)

In-place addition of vector b to a
"""
function add!(a::GF2Vector, b::GF2Vector)
    length(a) != length(b) && throw(ArgumentError("vectors have different lengths $(length(A)) and $(length(B))"))
    unsafe_add!(a, b)
end

"""
    unsafe_add!(a, b)

In-place addition of vector b to a, without checking that they are the same size
"""
function unsafe_add!(a::GF2Vector, b::GF2Vector)
    @inbounds for i in 1:length(a.chunks)
        a.chunks[i] ⊻= b.chunks[i]
    end
    return a
end

"""
    add!(a, b)

In-place addition of vector b to a
"""
function add!(a::GF2Vector, b::SubArray{Bool, 1, <:GF2Array, I, L}) where I where L
    length(a) != length(b) && throw(ArgumentError("vectors have different lengths $(length(A)) and $(length(B))"))
    unsafe_add!(a, b)
end

"""
    unsafe_add!(a, b)

In-place addition of vector b to a, without checking that they are the same size
"""
function unsafe_add!(a::GF2Vector, b::SubArray{Bool, 1, <:GF2Array, I, L}) where I where L
    chunkoffset = getchunkidx(b.offset1, size(a.chunks)[1], length(b))
    @inbounds for i in 1:length(a.chunks)
        a.chunks[i] ⊻= b.parent.chunks[i + chunkoffset]
    end
    return a
end

"""
compute which chunk of the parent to read from
"""
@inline function getchunkidx(i::Int, nc::Int, d1::Int)
     return ((i-1)÷d1)*nc + 1
end

"""
    addcol!(a, b, i)

In-place addition of column i of matrix b to a
"""
function addcol!(a::GF2Vector, b::GF2Matrix, i::Int)
    length(a) != size(b)[1] && throw(ArgumentError("dimension mismatch: vector $(length(A)) and matrix $(size(B))"))
    @boundscheck i > size(b)[2] && throw(BoundsError("column index $i greater than number of columns $(size(b)[2])"))
    unsafe_addcol!(a, b, i)
end

"""
    unsafe_addcol!(a, b, i)
"""
function unsafe_addcol!(a::GF2Vector, b::GF2Matrix, i::Int)
    chunkoffset = (i-1)*size(b.chunks)[1]
    @inbounds for i in 1:length(a.chunks)
        a.chunks[i] ⊻= b.chunks[i + chunkoffset]
    end
    return a
end

"""
    addcol!(A, i, j)

In-place addition of column j of A to column i of the same matrix
"""
function addcol!(A::GF2Matrix, i::Int, j::Int)
    @boundscheck i > size(b)[2] && throw(BoundsError("column index i=$i greater than number of columns $(size(b)[2])"))
    @boundscheck j > size(b)[2] && throw(BoundsError("column index j=$j greater than number of columns $(size(b)[2])"))
    unsafe_addcol!(A, b, i)
end

"""
    unsafe_addcol!(A, i, j)
"""
function unsafe_addcol!(A::GF2Matrix, i::Int, j::Int)
    nchk = size(A.chunks)[1]
    ichunkoffset = (i-1)*nchk
    jchunkoffset = (j-1)*nchk
    @inbounds for chunkidx in 1:nchk
        A.chunks[chunkidx + ichunkoffset] ⊻= A.chunks[chunkidx + jchunkoffset]
    end
    return A
end

"""
Constructor from BitArray
"""
GF2Array(A::AbstractArray{<:Any,N}) where {N} = GF2Array{N}(A)
function GF2Array{N}(a::BitArray{N}) where N
    r = GF2Array(undef, size(a)...)
    l = length(r)
    d1 = size(r)[1]
    nc = size(r.chunks)[1]
    l == 0 && return r

    idx = 1
    @inbounds begin
        for i = 1:length(r.chunks)
            stop = i%nc == 0 ? d1%64 : 64
            c = UInt64(0)
            for j = 0:stop-1
                c |= UInt64(a[idx]) << j
                idx += 1
            end
            r.chunks[i] = c
        end
    end

    return r
end

"""
Constructor from arbitrary array
"""
function GF2Array(a::AbstractArray{T,N}) where N where T
    b = BitArray(a)
    return GF2Array{N}(b)
end

# function Base.show(io::IO, a::GF2Vector)
#     vstr = ""
#     for i in 1:length(a.chunks)-1
#         vstr *= string(a.chunks[i], base=2, pad=64)
#     end
#     vstr *= string(a.chunks[end], base=2, pad=((length(a)-1)%64)+1)
#     print(io, "GF2Vector([$(join(vstr, ' '))])")
# end

isassigned(x::GF2Array, i::Int) = 1 <= i <= length(x)
IndexStyle(::Type{<:GF2Array}) = IndexLinear()

@inline function get_chunks_id(i::Integer, d1::Integer, nc::Integer)
    local_idx = (i-1)%d1
    ((i-1)÷d1)*nc + (local_idx÷64) + 1, local_idx%64
end

@inline function unsafe_bitgetindex(chunks::Array{UInt64}, i::Int, d1::Int)
    i1, i2 = get_chunks_id(i, d1, size(chunks)[1])
    u = UInt64(1) << i2
    @inbounds r = (chunks[i1] & u) != 0
    return r
end

@inline function getindex(x::GF2Array, i::Int)
    @boundscheck checkbounds(x, i)
    unsafe_bitgetindex(x.chunks, i, size(x)[1])
end

@inline function unsafe_bitsetindex!(chunks::Array{UInt64}, x::Bool, i::Int, d1::Int)
    i1, i2 = get_chunks_id(i, d1, size(chunks)[1])
    _unsafe_bitsetindex!(chunks, x, i1, i2)
end

@inline function _unsafe_bitsetindex!(chunks::Array{UInt64}, x::Bool, i1::Int, i2::Int)
    u = UInt64(1) << i2
    @inbounds begin
        c = chunks[i1]
        chunks[i1] = ifelse(x, c | u, c & ~u)
    end
end

@inline function setindex!(a::GF2Array, x, i::Int)
    @boundscheck checkbounds(a, i)
    unsafe_bitsetindex!(a.chunks, convert(Bool, x), i, size(a)[1])
    return a
end
