using Revise
using Test
using DiracBilinears
import LinearAlgebra as LA
import FFTW

@testset "Read wfc" begin
    a = 8.422887281
    cona = 1.330207305
    a1 = a*[1.0, 0.0, 0.0]
    a2 = a*[-1.0/2.0, √(3.0)/2.0, 0.0]
    a3 = a*[0.0, 0.0, cona]
    b1ref = 2π*(LA.cross(a2, a3)./LA.dot(a1, LA.cross(a2, a3)))
    b2ref = 2π*(LA.cross(a3, a1)./LA.dot(a2, LA.cross(a3, a1)))
    b3ref = 2π*(LA.cross(a1, a2)./LA.dot(a3, LA.cross(a1, a2)))

    wfc = read_wfc("./tmp/wfc1.dat")
    @test wfc.ik == 1
    @test wfc.xk == [0.0, 0.0, 0.0]
    @test wfc.npol == 2
    @test wfc.nbnd == 58
    @test isapprox(wfc.b1, b1ref, atol=1e-5)
    @test isapprox(wfc.b2, b2ref, atol=1e-5)
    @test isapprox(wfc.b3, b3ref, atol=1e-5)
end

@testset "Read xml" begin
    xml = read_xml("./tmp/data-file-schema.xml")
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
    xml = read_xml("./tmp/data-file-schema.xml")
    nrmesh = make_nrmesh(xml=xml, nrmesh=(0,0,0))
    @test nrmesh == (48,48,60)
    nrmesh = make_nrmesh(xml=xml, nrmesh=(32,32,32))
    @test nrmesh == (32,32,32)
end

@testset "density utils" begin
    # calc_occupation: step
    e = [-1.0, 0.0, 1.0]
    ef = 0.0
    occ = calc_occupation(e; ef=ef, smearing="step", degauss=0.1, δμ=0.0, emin=100000.0)
    @test occ == [1.0, 0.0, 0.0]
    # δμ shifts the threshold upward (in eV)
    occ_shift = calc_occupation(e; ef=ef, smearing="step", degauss=0.1, δμ=DiracBilinears.hartree2ev, emin=100000.0)
    @test occ_shift == [1.0, 1.0, 0.0]
    # emin excludes deep levels (in eV)
    occ_emin = calc_occupation(e; ef=ef, smearing="step", degauss=0.1, δμ=0.0, emin=0.0)
    @test occ_emin == [0.0, 0.0, 0.0]

    # make_zeros_density shapes
    @test size(make_zeros_density("ρ", (2, 3, 4))) == (2, 3, 4)
    @test size(make_zeros_density("ms", (2, 3, 4))) == (3, 2, 3, 4)
    @test size(make_zeros_density("j", (2, 3, 4))) == (3, 2, 3, 4)
    @test size(make_zeros_density("∇ρ", (2, 3, 4))) == (3, 2, 3, 4)
    @test size(make_zeros_density("∇ms", (2, 3, 4))) == (2, 3, 4)
    @test size(make_zeros_density("τz", (2, 3, 4))) == (2, 3, 4)
    @test size(make_zeros_density("ps", (2, 3, 4))) == (3, 2, 3, 4)

    # is_2d
    @test is_2d((1, 4, 4)) == (true, true, false, false)
    @test is_2d((4, 1, 4)) == (true, false, true, false)
    @test is_2d((4, 4, 1)) == (true, false, false, true)
    @test is_2d((4, 4, 4)) == (false, false, false, false)

    # calc_fourier_k behavior for reduced dimensions
    ck = reshape(ComplexF64.(1:4), (1, 2, 2))
    ukn1 = calc_fourier_k((1, 2, 2), ck)
    ukn1_ref = FFTW.bfft(sum(ck, dims=1), [2, 3])
    @test ukn1 == ukn1_ref

    ck2 = reshape(ComplexF64.(1:4), (2, 1, 2))
    ukn2 = calc_fourier_k((2, 1, 2), ck2)
    ukn2_ref = FFTW.bfft(sum(ck2, dims=2), [1, 3])
    @test ukn2 == ukn2_ref

    ck3 = reshape(ComplexF64.(1:4), (2, 2, 1))
    ukn3 = calc_fourier_k((2, 2, 1), ck3)
    ukn3_ref = FFTW.bfft(sum(ck3, dims=3), [1, 2])
    @test ukn3 == ukn3_ref
end

@testset "density calculations" begin
    nrmesh = (1, 1, 1)
    wfc = Wfc(
        1,
        [0.0, 0.0, 0.0],
        1,
        false,
        1.0,
        1,
        1,
        2,
        1,
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        zeros(Int, (3, 1)),
        zeros(ComplexF64, (2, 1, 1)),
    )
    occ = [1.0]

    # ρ and ms with a simple spinor
    ukn1 = reshape(ComplexF64[1.0 + 0.0im, 0.0 + 1.0im], (1, 1, 1, 2, 1))
    rho = calc_density_ρ(copy(ukn1), occ)
    @test isapprox(rho[1, 1, 1], 2.0, atol=1e-12)
    ms = calc_density_ms(copy(ukn1), occ)
    @test size(ms) == (3, 1, 1, 1)
    @test isapprox(ms[:, 1, 1, 1], [0.0, -2.0, 0.0], atol=1e-12)

    # ∇ρ with imaginary gradient along x
    ∇ukn1 = zeros(ComplexF64, (1, 1, 1, 2, 1, 3))
    ∇ukn1[:, :, :, :, :, 1] .= im .* ukn1
    grad_rho = calc_density_∇ρ(wfc, copy(ukn1), ∇ukn1, occ)
    @test isapprox(grad_rho[:, 1, 1, 1], [-4.0, 0.0, 0.0], atol=1e-12)

    # ∇ms with imaginary gradient along y
    ∇ukn2 = zeros(ComplexF64, (1, 1, 1, 2, 1, 3))
    ∇ukn2[:, :, :, :, :, 2] .= im .* ukn1
    grad_ms = calc_density_∇ms(wfc, copy(ukn1), ∇ukn2, occ)
    @test isapprox(grad_ms[1, 1, 1], 4.0, atol=1e-12)

    # j, τz, ps with a real spinor
    ukn2 = reshape(ComplexF64[1.0 + 0.0im, 1.0 + 0.0im], (1, 1, 1, 2, 1))
    ∇ukn3 = zeros(ComplexF64, (1, 1, 1, 2, 1, 3))
    ∇ukn3[:, :, :, :, :, 1] .= ukn2

    j = calc_density_j(nrmesh, wfc, copy(ukn2), ∇ukn3, occ)
    @test isapprox(j[:, 1, 1, 1], [4.0, 0.0, 0.0], atol=1e-12)

    τz = calc_density_τz(nrmesh, wfc, copy(ukn2), ∇ukn3, occ)
    @test isapprox(τz[1, 1, 1], 4.0, atol=1e-12)

    ps = calc_density_ps(nrmesh, wfc, copy(ukn2), ∇ukn3, occ)
    @test isapprox(ps[:, 1, 1, 1], [0.0, 0.0, 0.0], atol=1e-12)
end


@testset "Make ck" begin
    mill = zeros(Int, (3, 2))
    evc = ones(ComplexF64, (2, 2, 2))
    wfc = Wfc(1, [0.0, 0.0, 0.0], 2, false, 1.0, 2, 2, 2, 2, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], mill, evc)
    ck, ∇ck = make_c_k((4, 4, 4), wfc, is∇u=true)
    @test size(ck) == (4, 4, 4, 2, 2)
    @test isapprox(∇ck, zeros(ComplexF64, (4, 4, 4, 2, 2, 3)), atol=1e-12)
end

@testset "R-grid" begin
    n111, degen1 = calc_rgrid(mpmesh=(1, 1, 1))
    @test n111[:, 1] == [0, 0, 0]
    n222, degen2 = calc_rgrid(mpmesh=(2, 2, 2))
    @test n222[:, 1] == [-1, -1, -1]
    @test n222[:, 2] == [-1, -1, 0]
    @test n222[:, end] == [0, 0, 0]
    n333, degen3 = calc_rgrid(mpmesh=(3, 3, 3))
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
    a1, a2, a3 = read_lattice("./tmp/dat.lattice")
    b1, b2, b3 = calc_b(a1, a2, a3)
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
