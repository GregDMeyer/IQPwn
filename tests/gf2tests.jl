
using Test
using LinearAlgebra
import Random
Random.seed!(0xBEEFCAFE)

include("../solver/gf2.jl")

@testset "GF2Array" begin

    @testset "constructors" begin
        @test GF2Vector == GF2Array{1}
        @test GF2Matrix == GF2Array{2}

        GF2Array([])
        GF2Array([true])
        GF2Array([true false true])

        @testset "data l=$l $t" for
            l in (1, 64, 100, 128, 1000, (2,2), (64, 64), (100, 128)),
            t in (Random.bitrand, falses, trues)

            x = t(l)
            xg = GF2Array(x)
            @test size(x) == size(xg)
            @test length(x) == length(xg)
            @test all(x .== xg)
        end
    end

    @testset "getsetindex" begin
        @testset "values l=$l $t" for
            l in (1, 64, 100, 128, 1000, (2,2), (64, 64), (100, 128)),
            t in (Random.bitrand, falses, trues)

            x = t(l)
            xg = GF2Array(undef, size(x))

            xg[:] = x[:]
            @test all(x .== xg)

            @test_throws BoundsError xg[0]
            @test_throws BoundsError xg[length(xg)+1]
        end
    end

    @testset "dot" begin

        xg = GF2Array(falses(10))
        yg = GF2Array(falses(11))
        @test_throws DimensionMismatch dot(xg, yg)

        @testset "value l=$l" for
            l in (1, 64, 100, 128, 1000)

            x = Random.bitrand(l)
            y = Random.bitrand(l)

            xg = GF2Array(x)
            yg = GF2Array(y)

            @test sum(x .& y)%2 == dot(xg, yg)
        end
    end

    @testset "dotcol" begin

        xg = GF2Array(falses((10, 2)))
        yg = GF2Array(falses(11))
        @test_throws DimensionMismatch dotcol(yg, xg, 1)

        yg = GF2Array(falses(10))
        @test_throws BoundsError dotcol(yg, xg, 0)
        @test_throws BoundsError dotcol(yg, xg, 3)

        @testset "value l=$l" for
            l in (1, 64, 100, 128, 1000)

            x = Random.bitrand(l)
            y = Random.bitrand((l,2))

            xg = GF2Array(x)
            yg = GF2Array(y)

            @testset "col $i" for i in 1:2
                @test sum(x .& y[:,i])%2 == dotcol(xg, yg, i)
                @test sum(x .& y[:,i])%2 == unsafe_dotcol(xg, yg, i)
            end
        end
    end

    @testset "add" begin

        xg = GF2Array(falses(10))
        yg = GF2Array(falses(11))
        @test_throws DimensionMismatch add!(xg, yg)

        @testset "value l=$l" for
            l in (1, 64, 100, 128, 1000)

            x = Random.bitrand(l)
            y = Random.bitrand(l)

            xg = GF2Array(x)
            yg = GF2Array(y)

            x .⊻= y
            add!(xg, yg)
            @test all(x .== xg)
        end
    end

    @testset "addcol" begin
        xg = GF2Array(falses((10, 2)))
        yg = GF2Array(falses(11))
        @test_throws DimensionMismatch addcol!(yg, xg, 1)

        yg = GF2Array(falses(10))
        @test_throws BoundsError addcol!(yg, xg, 0)
        @test_throws BoundsError addcol!(yg, xg, 3)

        @test_throws BoundsError addcol!(xg, 0, 1)
        @test_throws BoundsError addcol!(xg, 3, 1)
        @test_throws BoundsError addcol!(xg, 1, 0)
        @test_throws BoundsError addcol!(xg, 1, 3)

        @testset "value l=$l" for
            l in (1, 64, 100, 128, 1000)

            x = Random.bitrand(l)
            x2 = copy(x)
            y = Random.bitrand((l,2))
            y2 = copy(y)

            xg = GF2Array(x)
            yg = GF2Array(y)

            @testset "col $i" for i in 1:2
                x .= x2
                y .= y2
                xg .= x2
                yg .= y2

                x .⊻= y[:,i]
                addcol!(xg, yg, i)
                @test all(x .== xg)

                y[:,i] .⊻= y[:,3-i]
                addcol!(yg, i, 3-i)
                @test all(y .== yg)
            end
        end
    end

end
