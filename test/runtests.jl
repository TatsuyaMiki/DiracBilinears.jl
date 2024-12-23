using Revise
using Test
import DiracBilinears as DB
import LinearAlgebra as LA

@testset "DiracBilinears.jl" begin
    # Write your tests here.
end

@testset "Read wfc" begin
    a = 8.422887281
    cona = 1.330207305
    a1 = a*[1.0, 0.0, 0.0]
    a2 = a*[-1.0/2.0, √(3.0)/2.0, 0.0]
    a3 = a*[0.0, 0.0, cona]
    b1ref = 2π*(LA.cross(a2, a3)./LA.dot(a1, LA.cross(a2, a3)))
    b2ref = 2π*(LA.cross(a3, a1)./LA.dot(a2, LA.cross(a3, a1)))
    b3ref = 2π*(LA.cross(a1, a2)./LA.dot(a3, LA.cross(a1, a2)))

    file = "./tmp/wfc1.dat"
    wfc = DB.read_wfc(file)
    @test wfc.ik == 1
    @test wfc.xk == [0.0, 0.0, 0.0]
    @test wfc.npol == 2
    @test wfc.nbnd == 58
    @test isapprox(wfc.b1, b1ref, atol=1e-5)
    @test isapprox(wfc.b2, b2ref, atol=1e-5)
    @test isapprox(wfc.b3, b3ref, atol=1e-5)
end
