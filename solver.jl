
include("gf2.jl")

using Random
using Base64

"""
    readprogram(filename)

Read an X-program from the file filename;
return a GF2Array containing *its transpose* (for efficiency).
"""
function readprogram(filename)
    lines = eachline(filename)
    nr = parse(Int, iterate(lines)[1][6:end])
    nc = parse(Int, iterate(lines)[1][6:end])

    # store it transposed for column-major efficiency
    rtn = GF2Array(undef, (nc, nr))

    for (i,line) in enumerate(lines)
        if (i > nr)
            break  # there are a few extra lines at the end
        end

        for (j,val) in enumerate(split(line, ' '))
            if (j > nc)
                break
            end

            rtn[j,i] = parse(Bool, val)
        end
    end

    rtn
end

"""
    arraytofile(bits, filename)

Write a BitArray bits to filename, in the format
that can be read by the challenge code to check
samples.
"""
function arraytofile(bits, filename)
    open(filename, write=true) do f
        write(f, "nr = $(size(bits)[2])\n")
        write(f, "nc = $(size(bits)[1])\n")
        for row in eachcol(bits) # we are transposed
            for c in row
                write(f, "$(Int(c)) ")
            end
            write(f, "\n")
        end
        write(f, "=====\n\n")
    end
end

"""
    papersoln(P, samples)

Given a program P as a GF2Array and a desired
number of samples nsamples, generate samples by
the classical method from the original paper, with
a success probability of 0.75 per sample.
"""
function papersoln(P, nsamples)
    n = size(P)[1]
    ncol = size(P)[2]

    rtn = falses((n, nsamples))
    for i in 1:nsamples
        d = GF2Array(bitrand(n))
        e = GF2Array(bitrand(n))
        tmp = GF2Array(falses(n))
        for colidx in 1:ncol
            if (unsafe_dotcol(d, P, colidx) + unsafe_dotcol(e, P, colidx) < 2)
                unsafe_addcol!(tmp, P, colidx)
            end
        end
        rtn[:, i] .= tmp
    end
    rtn
end

"""
    backsolve(system)

Given a system of equations over F2, return the set
of vectors that solve the system (perhaps more than 1
if the system is underconstrained). system should be a
BitArray in upper triangular form, with 1s along the diagonal
for every nonzero row (the matrix may contain some zero rows
for underconstrained systems).
"""
function backsolve(system)
    n = size(system)[2]
    solvemat = copy(system)

    # eliminate back up the chain
    for rowi in n:-1:1
        # if there is not already a row here, we need to guess
        if !(solvemat[rowi, rowi])

            # put a 1 in this slot
            solvemat[rowi, rowi] = true

            # append guesses to the end of the matrix
            s = copy(solvemat[n+1:end, :])
            s[:, rowi] .= true
            solvemat = [solvemat; s]
        end

        # eliminate back up the system
        for rowj in rowi-1:-1:1
            if solvemat[rowi, rowj]
                solvemat[:, rowj] .⊻= solvemat[:, rowi]
            end
        end
    end

    solvemat[n+1:end, :]
end

"""
    gensystem(P, maxiters)

Given an X-program P, generate a system of equations for s
in upper-triangular format. This routine has a 50% chance of
yielding a correct system of equations.
"""
function gensystem(P, maxiters)
    n = size(P)[1]
    ncol = size(P)[2]

    system = GF2Array(falses((n+1, n)))

    d = GF2Array(bitrand(n))

    rank = 0
    for i in 1:maxiters

        # generate a sample
        sample = GF2Array(falses(n))
        e = GF2Array(bitrand(n))
        for colidx in 1:ncol
            if (unsafe_dotcol(d, P, colidx) + unsafe_dotcol(e, P, colidx) < 2)
                unsafe_addcol!(sample, P, colidx)
            end
        end

        # add a 1 to the end
        sample = GF2Array([sample; [true]])

        # eliminate on it until it can fill a row that is empty
        for colidx in 1:n
            if sample[colidx]
                if system[colidx, colidx] # there is already data here, orthogonalize
                    addcol!(sample, system, colidx)
                else # we found a new linearly independent row!
                    system[:, colidx] .= sample
                    rank += 1
                    break
                end
            end
        end

        # check if we're done
        if rank == n
            break
        end
    end
    system
end

"""
    checkkey(P, s)

Check that a possibly correct key s is correct for the program P,
by extracting the sub-matrix corresponding to s and checking that
as a generator matrix its codewords are -1 or 0 mod 4
"""
function checkkey(P, s)
    Ns = 40
    n = size(P)[1]
    ncol = size(P)[2]
    for i in 1:Ns
        d = GF2Array(bitrand(n))
        tot = 0
        for colidx in 1:ncol
            if (unsafe_dotcol(s, P, colidx))
                tot += unsafe_dotcol(d, P, colidx)
            end
        end
        wtmod4 = tot%4
        if (wtmod4 != 0 && wtmod4 != 3)
            return false
        end
    end

    return true
end

"""
    extractkey(P, maxit=100)

Given an X-program as a BitArray, extract its secret key
Returns a tuple of the key and the number of candidate keys tested
"""
function extractkey(P; maxit=100, sysmaxit=1.2)

    keystried = 0
    for i in 1:maxit
        # generate a system of equations that should generate our key w/ 50% prob
        system = gensystem(P, floor(size(P)[1]*sysmaxit))

        # solve for the vectors satisfying that system
        candidate_keys = backsolve(system)

        # check if any of those keys were right
        for key in eachrow(candidate_keys)
            keystried += 1
            key = GF2Array(key)
            if checkkey(P, key)
                return GF2Array(key), keystried
            end
        end

        # if none of those keys were right, we got unlucky with d
        # try again
    end
    throw(ErrorException("max iterations reached"))
end

"""
    vectobin(v)

Convert a BitVector v to a string of 0's and 1's (for printing)
"""
function vectobin(v)
    s = ""
    for x in v
        s *= "$(Int(x))"
    end
    s
end

"""
    vectob64(v)

Convert a BitVector v to a base64 string
"""
function vectob64(v)
    # pad the vector to a multiple of 8 bits
    npad = 7 - ((length(v)-1) % 8)
    v = [falses(npad); v]

    s = Array{UInt8}(undef, (length(v)÷8,))
    for i in 1:8:length(v)
        sidx = (i-1)÷8 + 1
        s[sidx] = 0
        for (j, b) in enumerate(v[i:i+7])
            s[sidx] |= UInt8(b << (8-j))
        end
    end
    base64encode(s)
end

"""
    gensamples(s, nsamples)

Generate nsamples bitstrings such that ~85%
of them are non-orthogonal to a secret string s. Returns
a BitArray for which the columns are the samples
"""
function gensamples(s, nsamples)
    # length of samples
    n = size(s)[1]

    # when we get a sample that is non-orthogonal,
    # always include it. When we get a sample that is orthogonal, include
    # it with the appropriate probability

    # how often to include orthogonal samples so that the resulting fraction is right
    threshold = 1/(cos(pi/8)^2) - 1

    rtn = GF2Array(undef, (n, nsamples))

    count = 0
    while (count < nsamples)
        sample = GF2Array(bitrand(n))
        if (dot(sample, s) || rand()<threshold)
            count += 1
            rtn[:, count] = sample
        end
    end
    rtn
end
