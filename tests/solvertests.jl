
using Test
using LinearAlgebra
import Random
Random.seed!(0xBEEFCAFE)

include("../solver/solver.jl")

@testset "Solver" begin

    @testset "MatrixFileIO" begin
        s = IOBuffer()
        b = Random.bitrand((10,5))

        arraytofile(b, s)
        seek(s, 0)
        b2 = readprogram(s)

        @test b == b2
    end

    @testset "backsolve" begin
        n = 10

        x = transpose(bitrand(n))
        b = [I; x]
        @test backsolve(b) == x

        b[4,4] = 0
        b[end, 4] = 0

        r = [x; x]
        r[1,4] = 0
        r[2,4] = 1

        @test backsolve(b) == r

        for i in n:-1:2
            if i == 4
                continue
            end

            for j in i-1:-1:1
                if j == 4
                    continue
                end

                if rand(Bool)
                    b[i:end, j] .⊻= b[i:end, i]
                end
            end
        end
        @test backsolve(b) == r
    end

    @testset "checkkey" begin
        P = readprogram("test103.prog")
        s = "01001010010011010001101100111011001001111110110100101"
        v = [parse(Bool, c) for c in s]
        @test checkkey(P, v)

        @testset "bad" for i in 1:50
            @test !checkkey(P, bitrand(length(v)))
        end
    end

    @testset "vectobin" begin
        s = "01001010010011010001101100111011001001111110110100101"
        v = [parse(Bool, c) for c in s]
        @test vectobin(v) == s
    end

    @testset "vectob64" begin
        s = "01001010010011010001101100111011001001111110110100101"
        v = [parse(Bool, c) for c in s]
        c = "CUmjZ2T9pQ=="
        @test vectob64(v) == c
    end

end
