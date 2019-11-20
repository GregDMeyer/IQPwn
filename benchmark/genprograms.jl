"""
This file contains routines for generating X-programs as described in [Shepherd and Bremner 2009]
"""

include("../solver/gf2.jl")

using Primes

"""
    jacobisymbol(n, k)

Compute the Jacobi symbol of integers n and k.
Implementation adapted from https://en.wikipedia.org/wiki/Jacobi_symbol
"""
function jacobisymbol(n, k)
    (k < 1 || iseven(k)) && throw(ArgumentError("k must be an odd positive number"))

    n = mod(n, k)
    t = 1
    while n > 0
        while iseven(n)
            n ÷= 2
            r = k%8
            if r == 3 || r == 5
                t = -t
            end
        end
        n, k = k, n
        if n%4 == k%4 == 3
            t = -t
        end
        n %= k
    end
    if k == 1
        return t
    else
        return 0
    end
end
@assert jacobisymbol(3, 13) == 1
@assert jacobisymbol(6, 9) == 0
@assert jacobisymbol(-10, 13) == 1
@assert jacobisymbol(-3, 9) == 0

"""
    eliminatecolumns(A, cov)

Full Gauss-Jordan elimination of BitArray A and its covector over GF2
"""
function eliminatecolumns(A::GF2Matrix, cov)
    nrows = size(A)[1]
    ncols = size(A)[2]

    rank = 0
    pivrow = 0
    while rank < ncols
        # find a column with a 1 in the right spot
        rank += 1
        pivrow += 1
        while pivrow <= nrows && !any(A[pivrow, rank:end])
            pivrow += 1
        end
        if pivrow > nrows
            break
        end

        pivcol = rank
        while pivcol <= ncols && !(A[pivrow, pivcol])
            pivcol += 1
        end

        # swap it into the appropriate spot
        if pivcol != rank
            # XOR swap
            unsafe_addcol!(A, pivcol, rank)
            unsafe_addcol!(A, rank, pivcol)
            unsafe_addcol!(A, pivcol, rank)

            tmp = cov[pivcol]
            cov[pivcol] = cov[rank]
            cov[rank] = tmp
        end

        # now eliminate
        for colidx in rank+1:ncols
            if A[pivrow, colidx]
                unsafe_addcol!(A, colidx, rank)
                cov[rank] ⊻= cov[colidx]
            end
        end
    end

    # do the reduction step
    for pivcol in 2:rank
        # eliminate backwards
        pivrow = findfirst(A[:, pivcol])
        for colidx in pivcol-1:-1:1
            if A[pivrow, colidx]
                unsafe_addcol!(A, colidx, pivcol)
                cov[pivcol] ⊻= cov[colidx]
            end
        end
    end
end

"""
    genprogram(q)

Generate an X-program using a quadratic residue code with prime q.
This function performs the same operations as the C-code provided at
quantumchallenge.wordpress.com to generate X-programs.
"""
function genprogram(q; sortrows=true)
    @assert isprime(q)
    @assert q%8 == 7

    k = (q+1)÷2
    nrows = 2q
    ncols = k+1

    P = GF2Array(undef, (nrows, ncols))
    key = GF2Array(undef, (ncols,))

    # first column
    P[1:q, 1]     .= true
    P[q+1:end, 1] .= false
    key[1]         = true

    # the first q rows of the rest of the columns are the quadratic residue code
    # compute the Jacobi symbol of the indices
    for j in 1:k
        for i in 1:q
            P[i, j+1] = jacobisymbol(i-j+1, q) == 1
        end
    end

    # the rest of the rows are random bits
    P[q+1:end, 2:end] .= bitrand(q, k)

    # the key has all falses for the rest of the columns
    key[2:end] .= false

    # now shuffle the rows
    for i in 1:nrows
        newi = rand(1:nrows)
        tmp = P[newi, :]
        P[newi, :] = P[i, :]
        P[i, :] = tmp
    end

    # Gaussian eliminate on columns
    eliminatecolumns(P, key)

    if sortrows
        P = sortslices(P, dims=1)
    end

    return P, key
end
