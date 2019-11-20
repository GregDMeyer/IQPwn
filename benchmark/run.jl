#!/usr/bin/env julia

include("benchmark.jl")

qstart = 47
qend   = 1000
distq  = 487
iters  = 1000

println("Generating timing data for q in range $qstart-$qend, with $iters iterations")
timesolves(qend, iters, qstart=qstart)

println()
print("Generating candidate key distribution for q=$distq with $iters iterations... ")
gendistdata(distq, iters)
println("done")
