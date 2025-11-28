
function calc_wannier_matrix(;calc::String, wfndir::String, wandir::String, rgrid::Matrix{Int}, npol::Int=2)
    rg, rginv, nsymq, nnp = read_symmetry(wfndir*"/dat.symmetry")
    ski, ks, numirr, numrot, trs, rw, nkirr, ntk = read_sample_k(wfndir*"/dat.sample-k", nsymq, rg)
    ngi, ntg = read_nkm(wfndir*"/dat.nkm", nkirr)
    a1, a2, a3 = read_lattice(wfndir*"/dat.lattice")
    b1, b2, b3 = calc_b(a1, a2, a3)
    kgi = read_kg(wfndir*"/dat.kg", ngi, ntg, nkirr)
    ecut4psi = calc_ecut_for_ψ(kgi, ntg, nkirr, ngi, ski, b1, b2, b3)
    ngs, mills = calc_kg0s(ecut4psi, rg, rw, ngi, ntg, trs, ski, b1, b2, b3, ntk, numirr, numrot)
    nxk = length(ngs)
    nr = size(rgrid)[2]
    iowan = open(wandir*"/dat.wan", "r");
    seek(iowan, 4)
    nwfc = Int(read(iowan, Int32))
    o = make_zeros_wannier(calc, nr, nwfc)
    ax = axes(o)
    for ik in 1:nxk
        k = ks[:, ik]
        ng = ngs[ik]
        mill = mills[:, :, ik]
        cs = read_wan(iowan, nwfc, npol, ng)
        otmp = calc_wannier_ok(calc, cs, k, mill, b1, b2, b3, nwfc, ng, nxk)
        for ir in 1:nr
            o[ax[1:end-1]..., ir] += otmp .* exp(-im*2π*LA.dot(k, rgrid[:, ir]))
        end
    end
    close(iowan)
    return o
end

function calc_rgrid(;mpmesh::Tuple=(0,0,0), rfile::String="none")
    @assert rfile != "none" || mpmesh != (0, 0, 0) "No argument specified."
    if mpmesh == (0, 0, 0)
        mpmesh = read_wan_grid(rfile)
    end
    rgrid = zeros(Int, (3, mpmesh[1]*mpmesh[2]*mpmesh[3]))
    ir = 0
    for ir1 in 1:mpmesh[1]
        n1 = ir1 - div(mpmesh[1], 2) - 1
        for ir2 in 1:mpmesh[2]
            n2 = ir2 - div(mpmesh[2], 2) - 1
            for ir3 in 1:mpmesh[3]
                ir += 1
                n3 = ir3 - div(mpmesh[3], 2) - 1
                rgrid[1, ir] = n1
                rgrid[2, ir] = n2
                rgrid[3, ir] = n3
            end
        end
    end
    degen = ones(Int, mpmesh[1]*mpmesh[2]*mpmesh[3])
    return rgrid, degen
end

function read_hrdat_rgrid(hrfile::String)
    open(hrfile, "r") do io
        str = readline(io) # comment
        norb = parse(Int, readline(io)) # norb
        nr = parse(Int, readline(io)) # nr
        degen = zeros(Int, nr)
        rgrid = zeros(Int, (3, nr))
        for i in 1:div(nr, 15)
            d = parse.(Int, split(readline(io)))
            degen[(i-1)*15+1:i*15] = d
        end
        if 1 * (nr%15 != 0) == 1
            d = parse.(Int, split(readline(io)))
            degen[div(nr, 15)*15+1:end] = d
        end
        ir = 0
        for i in 1:nr*norb*norb
            str = readline(io)
            dat = parse.(Int, split(str)[1:5])
            if dat[4] == 1 && dat[5] == 1
                ir += 1
                rgrid[:, ir] = dat[1:3]
            end
        end
        return rgrid, degen
    end
end

function make_zeros_wannier(calc::String, nr::Int, nwfc::Int)
    if calc == "ρ"
        return zeros(ComplexF64, (nwfc, nwfc, nr))
    elseif calc == "ms"
        return zeros(ComplexF64, (3, nwfc, nwfc, nr))
    elseif calc == "j"
        return zeros(ComplexF64, (3, nwfc, nwfc, nr))
    elseif calc == "τz" || calc == "tau_z" || calc == "chirality"
        return zeros(ComplexF64, (nwfc, nwfc, nr))
    elseif calc == "ps"
        return zeros(ComplexF64, (3, nwfc, nwfc, nr))
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function calc_wannier_ok(calc::String, cs::Array{ComplexF64, 3}, k::Vector{Float64}, mill::Matrix{Int}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, nwfc::Int, igwx::Int, nxk::Int)
    if calc == "ρ"
        return calc_wan_ρ(cs, igwx, nwfc, nxk)
    elseif calc == "ms"
        return calc_wan_ms(cs, igwx, nwfc, nxk)
    elseif calc == "j"
        return calc_wan_j(cs, k, mill, b1, b2, b3, igwx, nwfc, nxk)
    elseif calc == "τz" || calc == "tau_z" || calc == "chirality"
        return calc_wan_τz(cs, k, mill, b1, b2, b3, igwx, nwfc, nxk)
    elseif calc == "ps"
        return calc_wan_ps(cs, k, mill, b1, b2, b3, igwx, nwfc, nxk)
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function calc_b(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64})
    b1 = 2π*LA.cross(a2, a3)/LA.dot(a1, LA.cross(a2, a3))
    b2 = 2π*LA.cross(a3, a1)/LA.dot(a2, LA.cross(a3, a1))
    b3 = 2π*LA.cross(a1, a2)/LA.dot(a3, LA.cross(a1, a2))
    return b1, b2, b3
end

function calc_wan_ρ(cs::Array{ComplexF64, 3}, igwx::Int, nwfc::Int, nxk::Int)
    ρk = zeros(ComplexF64, (nwfc, nwfc))
    for ipw in 1:igwx
        cuk = cs[ipw, :, :]
        ρk += ES.ein"sm,sn->mn"(conj.(cuk), cuk) ./ nxk
    end
    return ρk
end

function calc_wan_ms(cs::Array{ComplexF64, 3}, igwx::Int, nwfc::Int, nxk::Int)
    mk = zeros(ComplexF64, (3, nwfc, nwfc))
    for ipw in 1:igwx
        cuk = cs[ipw, :, :]
        mk += -ES.optein"sm,tn,sti->imn"(conj.(cuk), cuk, σ) ./ nxk
    end
    return mk
end

function calc_wan_j(cs::Array{ComplexF64, 3}, k::Vector{Float64}, mill::Matrix{Int}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, igwx::Int, nwfc::Int, nxk::Int)
    psk = zeros(ComplexF64, (3, nwfc, nwfc))
    for ipw in 1:igwx
        kgvec = (k[1] + mill[1,ipw])*b1 .+ (k[2] + mill[2,ipw])*b2 .+ (k[3] + mill[3,ipw])*b3
        cuk = cs[ipw, :, :]
        psk += 2.0 .* ES.optein"i,sm,sn->imn"(kgvec, conj.(cuk), cuk) ./ nxk
    end
    return psk
end

function calc_wan_τz(cs::Array{ComplexF64, 3}, k::Vector{Float64}, mill::Matrix{Int}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, igwx::Int, nwfc::Int, nxk::Int)
    τzk = zeros(ComplexF64, (nwfc, nwfc))
    for ipw in 1:igwx
        kgvec = (k[1] + mill[1,ipw])*b1 .+ (k[2] + mill[2,ipw])*b2 .+ (k[3] + mill[3,ipw])*b3
        cuk = cs[ipw, :, :]
        τzk += 2.0 .* ES.optein"i,sti,sm,tn->mn"(kgvec, σ, conj.(cuk), cuk) ./ nxk
    end
    return τzk
end

function calc_wan_ps(cs::Array{ComplexF64, 3}, k::Vector{Float64}, mill::Matrix{Int}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, igwx::Int, nwfc::Int, nxk::Int)
    psk = zeros(ComplexF64, (3, nwfc, nwfc))
    for ipw in 1:igwx
        kgvec = (k[1] + mill[1,ipw])*b1 .+ (k[2] + mill[2,ipw])*b2 .+ (k[3] + mill[3,ipw])*b3
        cuk = cs[ipw, :, :]
        psk += -2.0 .* ES.optein"ijk,j,stk,sm,tn->imn"(ϵijk, kgvec, σ, conj.(cuk), cuk) ./ nxk
    end
    return psk
end

function write_wannier_matrix(hr::Array{ComplexF64, 3}, rs::Matrix{Int}, degen::Vector{Int}; savefile::String)
    nr::Int = size(rs)[2]
    nwfc::Int = size(hr)[1]
    io = open(savefile,"w")
    PF.@printf(io, "%2s\n", "Written on "*"$(Dates.now())")
    PF.@printf(io, "%10d\n", nwfc)
    PF.@printf(io, "%10d\n", nr)
    for ir in 1:nr-1
        if ir % 15 == 0
            PF.@printf(io, "%5d\n", degen[ir])
        else
            PF.@printf(io, "%5d", degen[ir])
        end
    end
    PF.@printf(io, "%5d\n", degen[end])
    for ir in 1:nr
        for m in 1:nwfc
            for n in 1:nwfc
                PF.@printf(io, "%5d%5d%5d%5d%5d%25.17f%25.17f\n", rs[1, ir], rs[2, ir], rs[3, ir], n, m, real(hr[n, m, ir]), imag(hr[n, m, ir]))
            end
        end
    end
    close(io)
end
