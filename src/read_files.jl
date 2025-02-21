struct Wfc
    ik::Int # k-point index (1 to number of k-points)
    xk::Vector{Float64} # k-point coordinates
    ispin::Int # spin index for LSDA case: ispin=1 for spin-up, ispin=2 for spin-down for unpolarized or non-colinear cases, ispin=1 always
    Γonly::Bool # write or read only half of the plane waves
    scalef::Float64 # scale factor applied to wavefunctions
    ngw::Int # number of plane waves (PW)
    igwx::Int # max number of PW (may be larger than ngw, not sure why)
    npol::Int # number of spin states for PWs: 2 for non-colinear case, 1 otherwise
    nbnd::Int # number of wavefunctions
    b1::Vector{Float64}
    b2::Vector{Float64}
    b3::Vector{Float64} # primitive reciprocal lattice vectors
    mill::Matrix{Int} # miller indices: h=mill(1,i), k=mill(2,i), l=mill(3,i), the i-th PW has wave vector (k+G)(:)=xk(:)+h*b1(:)+k*b2(:)+ l*b3(:)
    evc::Array{ComplexF64, 3} # wave functions in the PW basis set, The first index runs on PW components, the second index runs on band states. For non-colinear case, each PW has a spin component first igwx components have PW with up spin, second igwx components have PW with down spin
end

struct Xml
    e::Matrix{Float64}
    ef::Float64
    nbnd::Int
    nxk::Int
    a1::Vector{Float64}
    a2::Vector{Float64}
    a3::Vector{Float64}
    fftgrid::Vector{Int}
    atoms::Vector{String}
    atomicpos::Matrix{Float64}
end

function read_wfc(filename::String)
    io = open(filename, "r");
    seek(io, 4)
    ik = Int(read(io, Int32))
    xk = zeros(Float64, 3)
    read!(io, xk)
    ispin = Int(read(io, Int32))
    Γonly = Bool(read(io, Int32))
    scalef = read(io, Float64)
    pos1 = mark(io)

    seek(io, pos1 + 8)
    ngw = Int(read(io, Int32))
    igwx = Int(read(io, Int32))
    npol = Int(read(io, Int32))
    nbnd = Int(read(io, Int32))
    pos2 = mark(io)

    seek(io, pos2 + 8)
    b1 = zeros(Float64, 3)
    read!(io, b1)
    b2 = zeros(Float64, 3)
    read!(io, b2)
    b3 = zeros(Float64, 3)
    read!(io, b3)
    pos3 = mark(io)

    seek(io, pos3 + 8)
    mill = zeros(Int32, (3, igwx))
    mill = Int.(read!(io, mill))
    pos4 = mark(io)

    seek(io, pos4 + 8)
    evc = zeros(ComplexF64, (npol, nbnd, igwx))
    for ib in 1:nbnd
        for is in 1:npol
            evc[is, ib, 1:end] = read!(io, zeros(ComplexF64, igwx))
        end
        pos5 = mark(io)
        seek(io, pos5 + 8)
    end
    close(io)
    return Wfc(ik, xk, ispin, Γonly, scalef, ngw, igwx, npol, nbnd, b1, b2, b3, mill, evc)
end

function read_xml(filename::String)
    doc = EzXML.readxml(filename)
    primates = EzXML.root(doc)
    nbnd = parse(Int, EzXML.nodecontent(findfirst("//nbnd/text()", primates)))
    nxk = parse(Int, EzXML.nodecontent(findfirst("//nks/text()", primates)))
    a1 = parse.(Float64, split(EzXML.nodecontent(findfirst("//a1/text()", primates))))
    a2 = parse.(Float64, split(EzXML.nodecontent(findfirst("//a2/text()", primates))))
    a3 = parse.(Float64, split(EzXML.nodecontent(findfirst("//a3/text()", primates))))

    astruct = findfirst("//atomic_structure", primates)
    apos = findfirst("atomic_positions", astruct)
    atom = findall("atom", apos)
    atoms = [att["name"] for att in atom]
    natom = parse(Int, astruct["nat"])
    atomicpos = zeros(3, natom)
    for ia in 1:natom
        atomicpos[:, ia] = parse.(Float64, split(EzXML.nodecontent(atom[ia]), r"\s+"))
    end

    fft = findfirst("//fft_grid", primates)
    fftgrid = parse.(Int, [fft["nr1"], fft["nr2"], fft["nr3"]])
    species = EzXML.nodecontent.(findall("//eigenvalues/text()", primates))
    e = zeros(Float64, (nbnd, nxk))
    for ik in 1:nxk
        e[1:end, ik] = parse.(Float64, split(species[ik]))
    end
    ef = parse(Float64, EzXML.nodecontent(findfirst("//fermi_energy/text()", primates)))
    return Xml(e, ef, nbnd, nxk, a1, a2, a3, fftgrid, atoms, atomicpos)
end
