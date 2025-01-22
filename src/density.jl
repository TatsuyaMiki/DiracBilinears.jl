
function calc_density(;calc::String, qedir::String, n1::Float64=0.0, n2::Float64=0.0, n3::Float64=0.0, nrmesh::Tuple=(0,0,0), δμ::Float64=0.0)
    xml = read_xml(qedir*"/data-file-schema.xml")
    nrmesh_ = make_nrmesh(xml=xml, nrmesh=nrmesh)
    o = make_zeros_density(calc, nrmesh_)
    volume = xml.nxk*abs(LA.dot(xml.a3, LA.cross(xml.a1, xml.a2)))
    for ik in 1:xml.nxk
        wfcfile = qedir*"/wfc$(ik).dat"
        wfc = read_wfc(wfcfile)

        idx = findall(x -> x < xml.ef + δμ/au2ev, xml.e[:, ik])
        occ = zeros(wfc.nbnd)
        occ[idx] .= 1.0
        ck, ∇ck = make_c_k(nrmesh_, wfc; n1=n1, n2=n2, n3=n3)
        ukn = calc_fourier_k(nrmesh_, ck)/√(volume)
        ∇ukn = calc_fourier_k(nrmesh_, ∇ck)/√(volume)
        ok = calc_density_ok(calc, nrmesh_, wfc, ukn, ∇ukn, occ)
        o += ok
    end
    return o
end

function make_nrmesh(;xml::Xml, nrmesh::Tuple=(0,0,0))
    if nrmesh == (0,0,0)
        return (div(xml.fftgrid[1], 2), div(xml.fftgrid[2], 2), div(xml.fftgrid[3], 2))
    else
        return nrmesh
    end
end


function make_zeros_density(calc::String, nrmesh::Tuple)
    if calc == "ρ" || calc == "rho"
        return zeros(nrmesh...)
    elseif calc == "ms"
        return zeros(3, nrmesh...)
    elseif calc == "j"
        return zeros(3, nrmesh...)
    elseif calc == "∇ρ" || calc == "nabla_rho"
        return zeros(3, nrmesh...)
    elseif calc == "∇ms" || calc == "nabla_ms"
        return zeros(nrmesh...)
    elseif calc == "τz" || calc == "tau_z" || calc == "chirality"
        return zeros(nrmesh...)
    elseif calc == "ps"
        return zeros(3, nrmesh...)
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function calc_density_ok(calc::String, nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    if calc == "ρ" || calc == "rho"
        return calc_density_ρ(ukn, occ)
    elseif calc == "ms"
        return calc_density_ms(ukn, occ)
    elseif calc == "j"
        return return calc_density_j(nrmesh, wfc, ukn, ∇ukn, occ)
    elseif calc == "∇ρ" || calc == "nabla_rho"
        return calc_density_∇ρ(wfc, ukn, ∇ukn, occ)
    elseif calc == "∇ms" || calc == "nabla_ms"
        return calc_density_∇ms(wfc, ukn, ∇ukn, occ)
    elseif calc == "τz" || calc == "tau_z" || calc == "chirality"
        return calc_density_τz(nrmesh, wfc, ukn, ∇ukn, occ)
    elseif calc == "ps"
        return calc_density_ps(nrmesh, wfc, ukn, ∇ukn, occ)
    else
        @assert false "Invalid value assigned to 'calc'."
    end
end

function calc_density_ρ(ukn, occ::Vector{Float64})
    return ES.ein"xyzsb,b->xyz"(abs2.(ukn), occ)
end

function calc_density_ms(ukn, occ::Vector{Float64})
    return ES.ein"xyzsb,sti,xyztb,b->ixyz"(conj.(ukn), σ, ukn, occ)
end

function calc_density_∇ρ(wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, :, ib] .*= occ[ib]
    end
    return 2.0*imag.(ES.ein"xyzsbi,xyzsb->ixyz"(conj.(∇ukn), ukn))
end

function calc_density_∇ms(wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, :, ib] .*= occ[ib]
    end
    return 2.0*imag.(ES.ein"xyzsbi,sti,xyztb->xyz"(conj.(∇ukn), σ, ukn))
end

function calc_density_j(nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @assert wfc.npol == 2
    kudu = zeros(ComplexF64, (3, nrmesh..., wfc.npol, wfc.nbnd))
    for ix in 1:3
        kudu[ix, :, :, :, :, :] = wfc.xk[ix]*ukn[:, :, :, :, :] .+ ∇ukn[:, :, :, :, :, ix]
    end
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, :, ib] .*= occ[ib]
    end
    return 2.0*real.(ES.ein"xyzsb,ixyzsb->ixyz"(conj.(ukn), kudu))
end

function calc_density_τz(nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @assert wfc.npol == 2
    kudu = zeros(ComplexF64, (3, nrmesh..., wfc.npol, wfc.nbnd))
    @inbounds for ix in 1:3
        kudu[ix, :, :, :, :, :] = wfc.xk[ix]*ukn[:, :, :, :, :] .+ ∇ukn[:, :, :, :, :, ix]
    end
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, :, ib] .*= occ[ib]
    end
    tmp = ES.ein"xyzsb,ixyztb->xyzsti"(conj.(ukn), kudu)
    return 2.0*real.(ES.ein"xyzsti,sti->xyz"(tmp, σ))
end

function calc_density_ps(nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @assert wfc.npol == 2
    kudu = zeros(ComplexF64, (3, nrmesh..., wfc.npol, wfc.nbnd))
    for ix in 1:3
        kudu[ix, :, :, :, :, :] = wfc.xk[ix]*ukn[:, :, :, :, :] .+ ∇ukn[:, :, :, :, :, ix]
    end
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, :, ib] .*= occ[ib]
    end
    tmp = ES.ein"xyzsb,ixyztb->xyzsti"(conj.(ukn), kudu)
    return -2.0*real.(ES.ein"ijk,stj,xyzstk->ixyz"(ϵijk, σ, tmp))
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
    hmin, hmax = extrema(wfc.mill[1,:])
    kmin, kmax = extrema(wfc.mill[2,:])
    lmin, lmax = extrema(wfc.mill[3,:])
    @assert hmax - hmin < nrmesh[1] || nrmesh[1] == 1 "'nrmesh' must be nrmesh >= FFT grid."
    @assert kmax - kmin < nrmesh[2] || nrmesh[2] == 1 "'nrmesh' must be nrmesh >= FFT grid."
    @assert lmax - lmin < nrmesh[3] || nrmesh[3] == 1 "'nrmesh' must be nrmesh >= FFT grid."
    ns1, ns2, ns3 = nrmesh
    i1 = abs(-div(nrmesh[1], 2) - hmin)
    i2 = abs(-div(nrmesh[2], 2) - kmin)
    i3 = abs(-div(nrmesh[3], 2) - lmin)
    ns1 = (hmax - hmin) * (nrmesh[1] == 1) + nrmesh[1] * (nrmesh[1] != 1)
    ns2 = (kmax - kmin) * (nrmesh[2] == 1) + nrmesh[2] * (nrmesh[2] != 1)
    ns3 = (lmax - lmin) * (nrmesh[3] == 1) + nrmesh[3] * (nrmesh[3] != 1)
    epn1 = exp(im*2π*n1)
    epn2 = exp(im*2π*n2)
    epn3 = exp(im*2π*n3)
    ihs = wfc.mill[1, :] .- hmin .+ 1
    iks = wfc.mill[2, :] .- kmin .+ 1
    ils = wfc.mill[3, :] .- lmin .+ 1
    cktmp = zeros(ComplexF64, (ns1, ns2, ns3, wfc.npol, wfc.nbnd))
    ∇cktmp = zeros(ComplexF64, (ns1, ns2, ns3, wfc.npol, wfc.nbnd, 3))
    @inbounds for ipw in 1:wfc.igwx
        gvec = wfc.mill[1, ipw] * wfc.b1 + wfc.mill[2, ipw] * wfc.b2 + wfc.mill[3, ipw] * wfc.b3
        epn = epn1^(wfc.mill[1, ipw] * (nrmesh[1] == 1)) * epn2^(wfc.mill[2, ipw] * (nrmesh[2] == 1)) * epn3^(wfc.mill[3, ipw] * (nrmesh[3] == 1))
        cktmp[i1 + ihs[ipw], i2 + iks[ipw], i3 + ils[ipw], :, :] .= wfc.evc[:, :, ipw] * epn
        @inbounds for ii in 1:3
            ∇cktmp[i1 + ihs[ipw], i2 + iks[ipw], i3 + ils[ipw], :, :, ii] .= gvec[ii] * wfc.evc[:, :, ipw] * epn
        end
    end
    return FFTW.ifftshift(cktmp, 1:3), FFTW.ifftshift(∇cktmp, 1:3)
end

function write_density(f0::Array{Float64, 3}; qedir::String="manual", savefile::String, atoms::Vector{String}=["none"], atomic_positions::Matrix{Float64}=zeros(3,2), a1::Vector{Float64}=zeros(Float64, 3), a2::Vector{Float64}=zeros(Float64, 3), a3::Vector{Float64}=zeros(Float64, 3))
    bohr2ang = 0.5291772083
    a1ang = zeros(Float64, 3)
    a2ang = zeros(Float64, 3)
    a3ang = zeros(Float64, 3)
    if qedir == "manual"
        a1ang = a1.*bohr2ang
        a2ang = a2.*bohr2ang
        a3ang = a3.*bohr2ang
    else
        @assert a1 == zeros(Float64, 3) && a2 == zeros(Float64, 3) && a3 == zeros(Float64, 3)
        xml = read_xml(qedir*"/data-file-schema.xml")
        a1ang = xml.a1.*bohr2ang
        a2ang = xml.a2.*bohr2ang
        a3ang = xml.a3.*bohr2ang
    end
    na1, na2, na3 = size(f0)
    fplot = zeros(Float64, (na1+1, na2+1, na3+1))
    fplot[1:end-1, 1:end-1, 1:end-1] = copy(f0)
    fplot[end, 1:end-1, 1:end-1] = copy(f0[1,:,:])
    fplot[1:end-1, end, 1:end-1] = copy(f0[:,1,:])
    fplot[1:end-1, 1:end-1, end] = copy(f0[:,:,1])
    fplot[1:end-1, end, end] = copy(f0[:,1,1])
    fplot[end, 1:end-1, end] = copy(f0[1,:,1])
    fplot[end, end, 1:end-1] = copy(f0[1,1,:])
    fplot[end, end, end] = copy(f0[1,1,1])
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
    natom = length(atoms)
    if atoms[1] != "none"
        PF.@printf(io, "%2s\n", "PRIMCOORD")
        PF.@printf(io, "%15f%10f\n", natom, length(unique(atoms)))
        for ia in 1:natom
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
    PF.@printf(io, "%15d%10d%10d\n", na1+1, na2+1, na3+1)
    PF.@printf(io, "%15f%10f%10f\n", 0.0, 0.0, 0.0)
    PF.@printf(io, "%15f%10f%10f\n", a1ang[1], a1ang[2], a1ang[3])
    PF.@printf(io, "%15f%10f%10f\n", a2ang[1], a2ang[2], a2ang[3])
    PF.@printf(io, "%15f%10f%10f\n", a3ang[1], a3ang[2], a3ang[3])
    ia = 0
    for ia3 in 1:na3+1
        for ia2 in 1:na2+1
            for ia1 in 1:na1+1
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