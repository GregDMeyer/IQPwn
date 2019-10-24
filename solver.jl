
using LinearAlgebra
using Random

"""
    readprogram(filename)

Read an X-program from the file filename;
return a BitArray containing *its transpose* (for efficiency).
"""
function readprogram(filename)
    lines = eachline(filename)
    nr = parse(Int, iterate(lines)[1][6:end])
    nc = parse(Int, iterate(lines)[1][6:end])

    # store it transposed for column-major efficiency
    rtn = BitArray(undef, (nc, nr))

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

Given a program P as a BitArray and a desired
number of samples nsamples, generate samples by
the classical method from the original paper, with
a success probability of 0.75 per sample.
"""
function papersoln(P, nsamples)
    n = size(P)[1]
    rtn = falses((n, nsamples))
    for i in 1:nsamples
        d = bitrand(n)
        e = bitrand(n)

        for p in eachcol(P)
            if ((dot(d, p)%2) + (dot(e, p)%2) < 2)
                rtn[:,i] .⊻= p
            end
        end
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
    for rowi in n:-1:2
        # if there is not already a row here, we need to guess
        if !(solvemat[rowi, rowi])
            # put a 1 in this slot
            solvemat[rowi, rowi] = true

            # append guesses to the end of the matrix
            s = copy(solvemat[n+1:end, :])
            s[n+1:end, rowi] .= true
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
    system = falses((n+1, n))

    d = bitrand(n)

    rank = 0
    for i in 1:maxiters

        # generate a sample
        sample = falses(n)
        e = bitrand(n)
        for p in eachcol(P)
            if ((dot(d, p)%2) + (dot(e, p)%2) < 2)
                sample .⊻= p
            end
        end

        # add a 1 to the end
        sample = [sample; [true]]

        # eliminate on it until it can fill a row that is empty
        for rowi in 1:n
            if sample[rowi]
                if system[rowi, rowi] # there is already data here, orthogonalize
                    sample .⊻= system[:, rowi]
                else # we found a new linearly independent row!
                    system[:, rowi] = sample
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
    for i in 1:Ns
        d = bitrand(n)
        tot = 0
        for p in eachcol(P)
            if (dot(s, p)%2) == 1
                tot += dot(d, p)%2
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
"""
function extractkey(P; maxit=100)
    for i in 1:maxit
        # generate a system of equations that should generate our key w/ 50% prob
        system = gensystem(P, size(P)[1]*5)

        # solve for the vectors satisfying that system
        candidate_keys = backsolve(system)

        # check if any of those keys were right
        for key in eachrow(candidate_keys)
            if checkkey(P, key)
                return key
            end
        end

        # if none of those keys were right, we try again
    end
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
    vectohex(v)

Convert a BitVector v to a hex string
"""
function vectohex(v)
    # pad the vector to a multiple of 8 bits
    npad = 7 - ((length(v)-1) % 8)
    v = [falses(npad); v]

    s = ""
    for i in 1:8:length(v)
        char = zero(UInt8)
        for (j, b) in enumerate(v[i:i+7])
            char |= UInt8(b << (8-j))
        end
        s *= "$(string(char, base=16))"
    end
    s
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

    rtn = BitArray(undef, (n, nsamples))

    count = 0
    while (count < nsamples)
        sample = bitrand(n)
        if (dot(sample, s)%2 == 1 || rand()<threshold)
            count += 1
            rtn[:, count] = sample
        end
    end
    rtn
end
