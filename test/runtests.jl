using Revise
using Test
import DiracBilinears as DB
import LinearAlgebra as LA

@testset "Read wfc" begin
    a = 8.422887281
    cona = 1.330207305
    a1 = a*[1.0, 0.0, 0.0]
    a2 = a*[-1.0/2.0, √(3.0)/2.0, 0.0]
    a3 = a*[0.0, 0.0, cona]
    b1ref = 2π*(LA.cross(a2, a3)./LA.dot(a1, LA.cross(a2, a3)))
    b2ref = 2π*(LA.cross(a3, a1)./LA.dot(a2, LA.cross(a3, a1)))
    b3ref = 2π*(LA.cross(a1, a2)./LA.dot(a3, LA.cross(a1, a2)))

    wfc = DB.read_wfc("./tmp/wfc1.dat")
    @test wfc.ik == 1
    @test wfc.xk == [0.0, 0.0, 0.0]
    @test wfc.npol == 2
    @test wfc.nbnd == 58
    @test isapprox(wfc.b1, b1ref, atol=1e-5)
    @test isapprox(wfc.b2, b2ref, atol=1e-5)
    @test isapprox(wfc.b3, b3ref, atol=1e-5)
end

@testset "Read xml" begin
    xml = DB.read_xml("./tmp/data-file-schema.xml")
    eref = [-1.098768132555416, -1.098768132555414, -1.098454799660214, -1.098454799660212, -1.097670850862821]
    @test isapprox(xml.e[1:5, 1], eref, atol=1e-12)
    @test isapprox(xml.ef, 3.061205827232420e-001, atol=1e-12)
    @test xml.nbnd == 58
    @test xml.nxk == 64
    a1ref = [8.692739547500000, 0.000000000000000, 0.000000000000000]
    a2ref = [-4.346369773750000, 7.528133276617179, 0.000000000000000]
    a3ref = [0.000000000000000, 0.000000000000000, 1.114938333300000e001]
    @test isapprox(xml.a1, a1ref, atol=1e-12)
    @test isapprox(xml.a2, a2ref, atol=1e-12)
    @test isapprox(xml.a3, a3ref, atol=1e-12)
end

@testset "nrmesh" begin
    xml = DB.read_xml("./tmp/data-file-schema.xml")
    nrmesh = DB.make_nrmesh(xml=xml, nrmesh=(0,0,0))
    @test nrmesh == (24,24,30)
    nrmesh = DB.make_nrmesh(xml=xml, nrmesh=(32,32,32))
    @test nrmesh == (32,32,32)
end


@testset "Make ck" begin
    mill = zeros(Int, (3, 2))
    evc = ones(ComplexF64, (2, 2, 2))
    wfc = DB.Wfc(1, [0.0, 0.0, 0.0], 2, false, 1.0, 2, 2, 2, 2, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], mill, evc)
    ck, ∇ck = DB.make_c_k((4, 4, 4), wfc)
    @test size(ck) == (4, 4, 4, 2, 2)
    @test isapprox(∇ck, zeros(ComplexF64, (4, 4, 4, 2, 2, 3)), atol=1e-12)
end

@testset "R-grid" begin
    n111, degen1 = DB.calc_rgrid(mpmesh=(1, 1, 1))
    @test n111[:, 1] == [0, 0, 0]
    n222, degen2 = DB.calc_rgrid(mpmesh=(2, 2, 2))
    @test n222[:, 1] == [-1, -1, -1]
    @test n222[:, 2] == [-1, -1, 0]
    @test n222[:, end] == [0, 0, 0]
    n333, degen3 = DB.calc_rgrid(mpmesh=(3, 3, 3))
    @test n333[:, 1] == [-1, -1, -1]
    @test n333[:, 2] == [-1, -1, 0]
    @test n333[:, 3] == [-1, -1, 1]
    @test n333[:, end] == [1, 1, 1]
end


# For RESPACK

# @testset "read dat.sample-k" begin
#     rg, rginv, nsymq, nnp = read_symmetry(wfndir*"/dat.symmetry")
#     ks, nxk = DB.read_sample_k("./tmp/dat.sample-k")
#     @test nxk == 3
#     @test maximum(abs, ks[:, 1] - [0.0, 0.0, 0.0]) < 1e-7
#     @test maximum(abs, ks[:, 2] - [0.0, 0.0, 1.0]) < 1e-7
#     @test maximum(abs, ks[:, 3] - [0.0, 0.0, -1.0]) < 1e-7
# end

@testset "lattice" begin
    a1, a2, a3 = DB.read_lattice("./tmp/dat.lattice")
    b1, b2, b3 = DB.calc_b(a1, a2, a3)
    @test maximum(abs, a1 - [8.695296059644800, 0.000000000000000, 0.000000000000000]) < 1e-7
    @test maximum(abs, a2 - [-4.347648029822400, 7.530347281079660, 0.000000000000000]) < 1e-7
    @test maximum(abs, a3 - [0.000000000000000, 0.000000000000000, 11.149501302368339]) < 1e-7
    @test abs(LA.dot(a1, b1) - 2π) < 1e-7
    @test abs(LA.dot(a2, b2) - 2π) < 1e-7
    @test abs(LA.dot(a3, b3) - 2π) < 1e-7
    @test abs(LA.dot(a1, b2)) < 1e-7
    @test abs(LA.dot(a2, b3)) < 1e-7
    @test abs(LA.dot(a3, b1)) < 1e-7
    @test abs(LA.dot(a1, b3)) < 1e-7
    @test abs(LA.dot(a2, b1)) < 1e-7
    @test abs(LA.dot(a3, b2)) < 1e-7
end

# @testset "read kg" begin
#     io = open("./tmp/dat.kg", "r");
#     gs, ng = DB.read_kg(io)
#     close(io)
#     @test ng == 3
#     @test gs[:, 1] == [0, 0, 0]
#     @test gs[:, 2] == [0, 0, -1]
#     @test gs[:, 3] == [0, 0, 1]
# end
