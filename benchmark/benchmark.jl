"""
This file contains functions for timing the solver and
extracting other data about its performance.
"""

include("../solver/solver.jl")
include("genprograms.jl")

using Primes

"""
    timesolve(q, iters)

Get timing data for how long it takes to extract the secret key from X-programs.
Time to generate the X-programs is not included in the timing results.
"""
function timesolve(q, niters; sysmaxit=0)
    times = Array{Float64}(undef, (niters,))
    candkeys = Array{Int64}(undef, (niters,))
    correct_keys = 0
    totkeystried = 0

    for it in 1:niters
        P, correctkey = genprogram(q, sortrows=false)

        if sysmaxit == 0
            (key, keystried, candsize), t, _, _, _ = @timed extractkey(GF2Array(transpose(P)))
        else
            (key, keystried, candsize), t, _, _, _ = @timed extractkey(GF2Array(transpose(P)), sysmaxit=sysmaxit)
        end

        times[it] = t
        candkeys[it] = candsize
        correct_keys += (key == correctkey)
        totkeystried += keystried
    end
    (times, correct_keys/niters, totkeystried/niters, candkeys)
end

"""
Find the next prime that is 7 mod 8
(for use in generating a quadratic residue code)
"""
function nextgoodprime(x)
    x = nextprime(x)
    while (x%8 != 7)
        x = nextprime(x+1)
    end
    x
end

"""
Find the first and third quartile of the input array
"""
function quartiles(ary)
    sort!(ary)
    l = length(ary)
    return (ary[l÷4], ary[3l÷4])
end

"""
    timesolves(maxq, iters; filename="timing.csv")

Record the average time to solve for each system size up to maxq, averaging over
iters iterations. Write the results to filename.
"""
function timesolves(maxq, niters; filename="timing.csv", qstart=47)
    q = nextgoodprime(qstart)
    qs = []
    correctfracs = []

    open(filename, "w") do f
        while (q < maxq)
            print("Running q=$q... ")
            times, correct, keyschecked, _ = timesolve(q, niters)
            println("done")

            push!(qs, q)
            push!(correctfracs, correct)

            quarts = quartiles(times)
            println(f, "$q,$(sum(times)/length(times)),$(quarts[1]),$(quarts[2]),$keyschecked")
            flush(f)

            # spread out the values of q in a nice way
            q = nextgoodprime(11*q÷10)
        end
    end

    # make sure everything was right!
    println()
    println("Fraction of keys correctly extracted:")
    for (q, correct) in zip(qs, correctfracs)
        println("$q : $correct")
    end
end

"""
    gencandkeydist(q, iters)

Test many X-programs, and count how large the set of candidate keys was for each.
Then select those that had large numbers of candidate keys, and rerun them, recording
how many candidate keys they had on the second run. Write the output to a file.
"""
function gendistdata(q, niters; threshold=16, filename="candkeys.csv")
    hard = []

    candkeys = Array{Int64}(undef, (niters,))
    for it in 1:niters
        P, correctkey = genprogram(q, sortrows=false)
        key, solveiters, nkeys = extractkey(GF2Array(transpose(P)), sysmaxit=10)
        if key != correctkey
            println("key incorrect!")
        end

        candkeys[it] = nkeys

        if nkeys >= threshold
            push!(hard, P)
        end
    end

    if !isempty(hard)
        # rerun the hard ones
        subiters = fld(niters, length(hard))
        hardcandkeys = Array{Int64}(undef, (subiters*length(hard),))
        for (n,P) in enumerate(hard)
            for it in 1:subiters
                (_, _, nkeys), _, _, _, _ = @timed extractkey(GF2Array(transpose(P)), sysmaxit=10)
                hardcandkeys[(n-1)*subiters + it] = nkeys
            end
        end
    else
        println("No 'hard' instances found.")
        hardcandkeys = Array{Int64}(undef, (0,))
    end

    # write output to file
    open(filename, "w") do f
        println(f, join((string(x) for x in candkeys), ","))

        if !isempty(hardcandkeys)
            println(f, join((string(x) for x in hardcandkeys), ","))
        else
            println(f, "")
        end
    end
end
