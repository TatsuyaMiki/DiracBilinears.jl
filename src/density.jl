
function calc_density(;calc::String, nrmesh::Tuple, qedir::String, n1::Float64=0.0, n2::Float64=0.0, n3::Float64=0.0, δμ::Float64=0.0)
    o = make_zeros_density(calc, nrmesh)
    xml = read_xml(qedir*"/data-file-schema.xml")
    volume = xml.nxk*abs(LA.dot(xml.a3, LA.cross(xml.a1, xml.a2)))
    for ik in 1:xml.nxk
        wfcfile = qedir*"/wfc$(ik).dat"
        wfc = read_wfc(wfcfile)

        idx = findall(x -> x < xml.ef + δμ/au2ev, xml.e[:, ik])
        occ = zeros(wfc.nbnd)
        occ[idx] .= 1.0

        ck, ∇ck = make_c_k(nrmesh, wfc; n1=n1, n2=n2, n3=n3)
        ukn = calc_fourier_k(nrmesh, ck)/√(volume)
        ∇ukn = calc_fourier_k(nrmesh, ∇ck)/√(volume)
        ok = calc_density_ok(calc, nrmesh, wfc, ukn, ∇ukn, occ)
        o += ok
    end
    return o
end

function make_zeros_density(calc::String, nrmesh::Tuple)
    if calc == "ρ"
        return zeros(nrmesh...)
    elseif calc == "ms"
        return zeros(3, nrmesh...)
    elseif calc == "τz" || calc == "chirality"
        return zeros(nrmesh...)
    elseif calc == "ps"
        return zeros(3, nrmesh...)
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function calc_density_ok(calc::String, nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    if calc == "ρ"
        return calc_density_ρ(ukn, occ)
    elseif calc == "∇ρ"
        return calc_density_∇ρ(wfc, ukn, ∇ukn, occ)
    elseif calc == "∇ms"
        return calc_density_∇ms(wfc, ukn, ∇ukn, occ)
    elseif calc == "ms"
        return calc_density_ms(ukn, occ)
    elseif calc == "τz" || calc == "chirality"
        return calc_density_τz(nrmesh, wfc, ukn, ∇ukn, occ)
    elseif calc == "ps"
        return calc_density_ps(nrmesh, wfc, ukn, ∇ukn, occ)
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function calc_density_ρ(ukn, occ::Vector{Float64})
    return ES.ein"xyzbs,b->xyz"(abs2.(ukn), occ)
end

function calc_density_ms(ukn, occ::Vector{Float64})
    return ES.ein"xyzbs,sti,xyzbt,b->ixyz"(conj.(ukn), σ, ukn, occ)
end

function calc_density_∇ρ(wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, ib, :] .*= occ[ib]
    end
    return 2.0*imag.(ES.ein"ixyzbs,xyzbs->xyzsti"(conj.(∇ukn), ukn))
end

function calc_density_∇ms(wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, ib, :] .*= occ[ib]
    end
    return 2.0*imag.(ES.ein"ixyzbs,ist,xyzbt->xyz"(conj.(∇ukn), σ, ukn))
end

function calc_density_τz(nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @assert wfc.npol == 2
    kudu = zeros(ComplexF64, (3, nrmesh..., wfc.nbnd, wfc.npol))
    @inbounds for ix in 1:3
        kudu[ix, :, :, :, :, :] = wfc.xk[ix]*ukn[:, :, :, :, :] .+ ∇ukn[:, :, :, :, ix, :]
    end
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, ib, :] .*= occ[ib]
    end
    tmp = ES.ein"xyzbs,ixyzbt->xyzsti"(conj.(ukn), kudu)
    return 2.0*real.(ES.ein"xyzsti,sti->xyz"(tmp, σ))
end

function calc_density_ps(nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @assert wfc.npol == 2
    kudu = zeros(ComplexF64, (3, nrmesh..., wfc.nbnd, wfc.npol))
    for ix in 1:3
        kudu[ix, :, :, :, :, :] = wfc.xk[ix]*ukn[:, :, :, :, :] .+ ∇ukn[:, :, :, :, ix, :]
    end
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, ib, :] .*= occ[ib]
    end
    tmp = ES.ein"xyzbs,ixyzbt->xyzsti"(conj.(ukn), kudu)
    return 2.0*real.(ES.ein"ijk,stj,xyzstk->ixyz"(ϵijk, σ, tmp))
end

function calc_fourier_k(nrmesh::Tuple, ck)
    if nrmesh[1] == 1
        cktmp = sum(ck, dims=1)
        ukn = FFTW.bfft(cktmp, [2,3])
        return ukn
    elseif nrmesh[2] == 1
        cktmp = sum(ck, dims=2)
        ukn = FFTW.bfft(cktmp, [1,3])
        return ukn
    elseif nrmesh[3] == 1
        cktmp = sum(ck, dims=3)
        ukn = FFTW.bfft(cktmp, [1,2])
        return ukn
    else
        ukn = FFTW.bfft(ck, [1,2,3])
        return ukn
    end
end

function make_c_k(nrmesh::Tuple, wfc::Wfc; n1::Float64=0.0, n2::Float64=0.0, n3::Float64=0.0)
    hmax = maximum(wfc.mill[1,:])
    hmin = minimum(wfc.mill[1,:])
    kmax = maximum(wfc.mill[2,:])
    kmin = minimum(wfc.mill[2,:])
    lmax = maximum(wfc.mill[3,:])
    lmin = minimum(wfc.mill[3,:])
    nh = hmax - hmin + 1
    nk = kmax - kmin + 1
    nl = lmax - lmin + 1
    ns1 = nrmesh[1]
    ns2 = nrmesh[2]
    ns3 = nrmesh[3]
    i1 = abs(-div(nrmesh[1], 2) - hmin)
    i2 = abs(-div(nrmesh[2], 2) - kmin)
    i3 = abs(-div(nrmesh[3], 2) - lmin)
    if nrmesh[1] == 1
        ns1 = nh
        i1 = 0
    else
        @assert nh < nrmesh[1]
        ns1 = nrmesh[1]
    end
    if nrmesh[2] == 1
        ns2 = nk
        i2 = 0
    else
        @assert nk < nrmesh[2]
        ns2 = nrmesh[2]
    end
    if nrmesh[3] == 1
        ns3 = nl
        i3 = 0
    else
        @assert nl < nrmesh[3]
        ns3 = nrmesh[3]
    end
    cktmp = zeros(ComplexF64, (ns1, ns2, ns3, wfc.nbnd, wfc.npol))
    ∇cktmp = zeros(ComplexF64, (ns1, ns2, ns3, wfc.nbnd, 3, wfc.npol))
    for ipw in 1:wfc.igwx
        gvec = wfc.mill[1,ipw]*wfc.b1 + wfc.mill[2,ipw]*wfc.b2 + wfc.mill[3,ipw]*wfc.b3
        ih = wfc.mill[1,ipw] - hmin + 1
        ik = wfc.mill[2,ipw] - kmin + 1
        il = wfc.mill[3,ipw] - lmin + 1
        epn = 1.0
        if nrmesh[1] == 1
            epn = exp(2π*im*n1*wfc.mill[1,ipw])*epn
        end
        if nrmesh[2] == 1
            epn = exp(2π*im*n2*wfc.mill[2,ipw])*epn
        end
        if nrmesh[3] == 1
            epn = exp(2π*im*n3*wfc.mill[3,ipw])*epn
        end
        for is in 1:wfc.npol
            evci = wfc.evc[is, :, ipw]
            cktmp[i1 + ih, i2 + ik, i3 + il, :, is] = evci*epn
            for ii in 1:3
                ∇cktmp[i1 + ih, i2 + ik, i3 + il, :, ii, is] = gvec[ii]*evci*epn
            end
        end
    end
    ck = zeros(ComplexF64, (ns1, ns2, ns3, wfc.nbnd, wfc.npol))
    ∇ck = zeros(ComplexF64, (ns1, ns2, ns3, wfc.nbnd, 3, wfc.npol))
    for is in 1:wfc.npol
        for ib in 1:wfc.nbnd
            ck[:, :, :, ib, is] = FFTW.ifftshift(cktmp[:, :, :, ib, is])
            for ii in 1:3
                ∇ck[:, :, :, ib, ii, is] = FFTW.ifftshift(∇cktmp[:, :, :, ib, ii, is])
            end
        end
    end
    return ck, ∇ck
end

function write_vesta_plot(f0::Array{Float64, 3}; qedir::String, savefile::String, atoms::Vector{String}=["none"], atomic_positions::Matrix{Float64}=zeros(3,2), r1::Int=1, r2::Int=1, r3::Int=1)
    xml = read_xml(qedir*"/data-file-schema.xml")
    na1, na2, na3 = size(f0)
    fplot = zeros(Float64, (na1*r1, na2*r2, na3*r3))
    for i3 in 1:r3
        for i2 in 1:r2
            for i1 in 1:r1
                fplot[na1*(i1-1)+1:na1*i1, na2*(i2-1)+1:na2*i2, na3*(i3-1)+1:na3*i3] = copy(f0)
            end
        end
    end
    bohr2ang = 0.5291772083
    a1ang = xml.a1.*bohr2ang
    a2ang = xml.a2.*bohr2ang
    a3ang = xml.a3.*bohr2ang
    io = open(savefile, "w")
    PF.@printf(io, "%2s\n", "# Written on "*"$(Dates.now())")
    PF.@printf(io, "%2s\n", "CRYSTAL")
    PF.@printf(io, "%2s\n", "PRIMVEC")
    PF.@printf(io, "%15f%10f%10f\n", a1ang[1], a1ang[2], a1ang[3])
    PF.@printf(io, "%15f%10f%10f\n", a2ang[1], a2ang[2], a2ang[3])
    PF.@printf(io, "%15f%10f%10f\n", a3ang[1], a3ang[2], a3ang[3])
    PF.@printf(io, "%2s\n", "CONVVEC")
    PF.@printf(io, "%15f%10f%10f\n", a1ang[1], a1ang[2], a1ang[3])
    PF.@printf(io, "%15f%10f%10f\n", a2ang[1], a2ang[2], a2ang[3])
    PF.@printf(io, "%15f%10f%10f\n", a3ang[1], a3ang[2], a3ang[3])
    if atoms[1] != "none"
        PF.@printf(io, "%2s\n", "PRIMCOORD")
        PF.@printf(io, "%15f%10f\n", length(atoms), length(unique(atoms)))
        for ia in 1:length(atoms)
            ra = atomic_positions[1, ia]*a1ang + atomic_positions[2, ia]*a2ang + atomic_positions[3, ia]*a3ang
            PF.@printf(io, "%2s%10f%10f%10f\n", atoms[ia], ra[1], ra[2], ra[3])
        end
    end
    PF.@printf(io, "%2s\n", "")
    PF.@printf(io, "%2s\n", "")
    PF.@printf(io, "%2s\n", "")
    PF.@printf(io, "%2s\n", "BEGIN_BLOCK_DATAGRID_3D")
    PF.@printf(io, "%2s\n", "3D_field")
    PF.@printf(io, "%2s\n", "BEGIN_DATAGRID_3D_UNKNOWN")
    PF.@printf(io, "%15d%10d%10d\n", na1*r1, na2*r2, na3*r3)
    PF.@printf(io, "%15f%10f%10f\n", 0.0, 0.0, 0.0)
    PF.@printf(io, "%15f%10f%10f\n", a1ang[1]*r1, a1ang[2]*r1, a1ang[3]*r1)
    PF.@printf(io, "%15f%10f%10f\n", a2ang[1]*r2, a2ang[2]*r2, a2ang[3]*r2)
    PF.@printf(io, "%15f%10f%10f\n", a3ang[1]*r3, a3ang[2]*r3, a3ang[3]*r3)
    ia = 0
    for ia3 in 1:na3*r3
        for ia2 in 1:na2*r2
            for ia1 in 1:na1*r1
                ia += 1
                if ia%6 == 0
                    PF.@printf(io, "%10f\n", fplot[ia1, ia2, ia3])
                else
                    PF.@printf(io, "%10f", fplot[ia1, ia2, ia3])
                end 
            end
        end
    end
    PF.@printf(io, "%2s\n", "")
    PF.@printf(io, "%2s\n", "END_DATAGRID_3D")
    PF.@printf(io, "%2s\n", "END_BLOCK_DATAGRID_3D")
    close(io)
end
