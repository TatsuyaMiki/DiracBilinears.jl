

function calc_wannier_matrix(;calc::String, nrs::Tuple, wfndir::String, wandir::String, npol::Int=2)
    rs = make_rmesh(nrs)

    rg, rginv, nsymq, nnp = read_symmetry(wfndir*"/dat.symmetry")
    ski, ks, numirr, numrot, trs, rw, nkirr, ntk = read_sample_k(wfndir*"/dat.sample-k", nsymq, rg)
    ngi, ntg = read_nkm(wfndir*"/dat.nkm", nkirr)
    a1, a2, a3 = read_lattice(wfndir*"/dat.lattice")
    b1, b2, b3 = calc_b(a1, a2, a3)
    kgi = read_kg(wfndir*"/dat.kg", ngi, ntg, nkirr)
    ecut4psi = calc_ecut_for_ψ(kgi, ntg, nkirr, ngi, ski, b1, b2, b3)
    ngs, mills = calc_kg0s(ecut4psi, rg, rw, ngi, ntg, trs, ski, b1, b2, b3, ntk, numirr, numrot)
    nxk = length(ngs)
    iowan = open(wandir*"/dat.wan", "r");
    seek(iowan, 4)
    nwfc = Int(read(iowan, Int32))
    o = make_zeros_wannier(calc, nrs, nwfc)
    ax = axes(o)
    for ik in 1:nxk
        # if ik % 100 == 0
        #     println("Calculating at ik=$(ik)...")
        # end
        println("Calculating at ik=$(ik)..."); flush(stdout)
        k = ks[:, ik]
        ng = ngs[ik]
        mill = mills[:, :, ik]
        cs = read_wan(iowan, nwfc, npol, ng)
        otmp = calc_wannier_ok(calc, cs, k, mill, b1, b2, b3, nwfc, ng, nxk)
        for ir in 1:(nrs[1]*nrs[2]*nrs[3])
            o[ax[1:end-1]..., ir] += otmp .* exp(-im*2π*LA.dot(k, rs[:, ir]))
        end
    end
    close(iowan)
    return o, rs
end

function make_zeros_wannier(calc::String, nrs::Tuple, nwfc::Int)
    if calc == "ρ"
        return zeros(ComplexF64, (nwfc, nwfc, nrs[1]*nrs[2]*nrs[3]))
    elseif calc == "ms"
        return zeros(ComplexF64, (3, nwfc, nwfc, nrs[1]*nrs[2]*nrs[3]))
    elseif calc == "τz" || calc == "chirality"
        return zeros(ComplexF64, (nwfc, nwfc, nrs[1]*nrs[2]*nrs[3]))
    elseif calc == "ps"
        return zeros(ComplexF64, (3, nwfc, nwfc, nrs[1]*nrs[2]*nrs[3]))
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function calc_wannier_ok(calc::String, cs::Array{ComplexF64, 3}, k::Vector{Float64}, mill::Matrix{Int}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, nwfc::Int, igwx::Int, nxk::Int)
    if calc == "ρ"
        return calc_wan_ρ(cs, igwx, nwfc, nxk)
    elseif calc == "ms"
        return calc_wan_ms(cs, igwx, nwfc, nxk)
    elseif calc == "τz" || calc == "chirality"
        return calc_wan_τz(cs, k, mill, b1, b2, b3, igwx, nwfc, nxk)
    elseif calc == "ps"
        return calc_wan_ps(cs, k, mill, b1, b2, b3, igwx, nwfc, nxk)
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function make_rmesh(nrs)
    nvec = zeros(Int, (3, nrs[1]*nrs[2]*nrs[3]))
    ir = 0
    for ir1 in 1:nrs[1]
        n1 = ir1 - div(nrs[1], 2) - 1
        for ir2 in 1:nrs[2]
            n2 = ir2 - div(nrs[2], 2) - 1
            for ir3 in 1:nrs[3]
                ir += 1
                n3 = ir3 - div(nrs[3], 2) - 1
                nvec[1, ir] = n1
                nvec[2, ir] = n2
                nvec[3, ir] = n3
            end
        end
    end
    return nvec
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
        mk += ES.ein"sm,tn,sti->imn"(conj.(cuk), cuk, σ) ./ nxk
    end
    return mk
end

function calc_wan_τz(cs::Array{ComplexF64, 3}, k::Vector{Float64}, mill::Matrix{Int}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, igwx::Int, nwfc::Int, nxk::Int)
    τzk = zeros(ComplexF64, (nwfc, nwfc))
    for ipw in 1:igwx
        kgvec = (k[1] + mill[1,ipw])*b1 .+ (k[2] + mill[2,ipw])*b2 .+ (k[3] + mill[3,ipw])*b3
        cuk = cs[ipw, :, :]
        τzk += 2.0 .* ES.ein"i,sti,sm,tn->mn"(kgvec, σ, conj.(cuk), cuk) ./ nxk
    end
    return τzk
end

function calc_wan_ps(cs::Array{ComplexF64, 3}, k::Vector{Float64}, mill::Matrix{Int}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, igwx::Int, nwfc::Int, nxk::Int)
    psk = zeros(ComplexF64, (3, nwfc, nwfc))
    for ipw in 1:igwx
        kgvec = (k[1] + mill[1,ipw])*b1 .+ (k[2] + mill[2,ipw])*b2 .+ (k[3] + mill[3,ipw])*b3
        cuk = cs[ipw, :, :]
        psk += 2.0 .* ES.ein"ijk,j,stk,sm,tn->imn"(ϵijk, kgvec, σ, conj.(cuk), cuk) ./ nxk
    end
    return psk
end
