
function calc_density(;calc::String, qedir::String, n1::Float64=0.0, n2::Float64=0.0, n3::Float64=0.0, nrmesh::Tuple=(0,0,0), δμ::Float64=0.0, emin::Float64=100000.0, smearing::String="step", degauss::Float64=0.01)
    ## - n1, n2, n3 are parameters which can be used for calculations in a 2D plane.
    ##   The calculations are performed on a plane perpendicular to ai (i=1,2,3) that passes through ni*ai.
    ## - The chemical potential can be shifted using δμ (eV).
    ## - The core levels can be excluded by using emin (eV).
    is∇u = is∇ukn(calc)
    xml = read_xml(qedir*"/data-file-schema.xml")
    nrmesh_ = make_nrmesh(xml=xml, nrmesh=nrmesh)
    o = make_zeros_density(calc, nrmesh_)
    volume = xml.nxk*abs(LA.dot(xml.a3, LA.cross(xml.a1, xml.a2)))
    for ik in 1:xml.nxk
        wfc = qewfc(ik, xml, qedir)
        occ = calc_occupation(xml.e[:, ik]; ef=xml.ef, smearing=smearing, degauss=degauss, δμ=δμ, emin=emin)
        ck, ∇ck = make_c_k(nrmesh_, wfc; n1=n1, n2=n2, n3=n3, is∇u=is∇u)
        ukn = calc_fourier_k(nrmesh_, ck)/√(volume)
        ∇ukn = is∇u == true ? calc_fourier_k(nrmesh_, ∇ck)/√(volume) : ∇ck
        o += calc_density_ok(calc, nrmesh_, wfc, ukn, ∇ukn, occ)
    end
    return o
end

is∇ukn(calc::String) = calc in ["j", "∇ρ", "nabla_rho", "∇ms", "nabla_ms", "τz", "tau_z", "chirality", "ps"] 

function qewfc(ik::Int, xml::Xml, qedir::String)
    if xml.slda == "false"
        wfc = read_wfc(qedir*"/wfc$(ik).dat")
    elseif xml.slda == "true" && xml.noncolin == "false"
        wfcup = read_wfc(qedir*"/wfcup$(ik).dat")
        wfcdw = read_wfc(qedir*"/wfcdw$(ik).dat")
        npol = 2
        nbnd = 2wfcup.nbnd
        igwx = max(wfcup.igwx, wfcdw.igwx)
        evc = zeros(ComplexF64, (npol, nbnd, igwx))
        evc[1, 1:wfcup.nbnd, :] = wfcup.evc
        evc[2, wfcup.nbnd+1:end, :] = wfcdw.evc
        wfc = Wfc(wfcup.ik, wfcup.xk, wfcup.ispin, wfcup.Γonly, wfcup.scalef, wfcup.ngw, igwx, npol, nbnd, wfcup.b1, wfcup.b2, wfcup.b3, wfcup.mill, evc)
    end
    return wfc
end

function calc_occupation(e::Vector{Float64}; ef::Float64, smearing::String="step", degauss::Float64=1.0, δμ::Float64=0.0, emin::Float64=100000.0)
    if (smearing == "m-p") || (smearing == "mp")
        occ = methfessel_paxton_step.((e .- (ef + δμ/hartree2ev)) ./ degauss; n=1)
    elseif smearing == "step"
        idx = findall(x -> (x < ef + δμ/hartree2ev), e)
        occ = zeros(Float64, size(e))
        occ[idx] .= 1.0
    else
        @assert false "Invalid value assigned to 'smearing'."
    end
    if emin < 100000.0
        idx = findall(x -> (x < ef + emin/hartree2ev), e)
        occ[idx] .= 0.0
    end
    return occ
end

function make_nrmesh(;xml::Xml, nrmesh::Tuple=(0,0,0))
    if nrmesh == (0,0,0)
        return (xml.fftgrid[1], xml.fftgrid[2], xml.fftgrid[3])
    else
        return nrmesh
    end
end


function make_zeros_density(calc::String, nrmesh::Tuple)
    if calc == "ρ" || calc == "rho"
        return zeros(Float64, nrmesh)
    elseif calc == "ms"
        return zeros(Float64, (3, nrmesh...))
    elseif calc == "j"
        return zeros(Float64, (3, nrmesh...))
    elseif calc == "∇ρ" || calc == "nabla_rho"
        return zeros(Float64, (3, nrmesh...))
    elseif calc == "∇ms" || calc == "nabla_ms"
        return zeros(Float64, nrmesh)
    elseif calc == "τz" || calc == "tau_z" || calc == "chirality"
        return zeros(Float64, nrmesh)
    elseif calc == "ps"
        return zeros(Float64, (3, nrmesh...))
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
    return -real.(ES.optein"xyzsb,sti,xyztb,b->ixyz"(conj.(ukn), σ, ukn, occ))
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
    return -2.0*imag.(ES.optein"xyzsbi,sti,xyztb->xyz"(conj.(∇ukn), σ, ukn))
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
    return -2.0*real.(ES.optein"ijk,stj,xyzstk->ixyz"(ϵijk, σ, tmp))
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

function is_2d(nrmesh::Tuple)
    is1 = nrmesh[1] == 1
    is2 = nrmesh[2] == 1
    is3 = nrmesh[3] == 1
    return is1 || is2 || is3, is1, is2, is3
end

function make_c_k(nrmesh::Tuple, wfc::Wfc; n1::Float64=0.0, n2::Float64=0.0, n3::Float64=0.0, is∇u=false)
    is2d, is1, is2, is3 = is_2d(nrmesh)
    hmin, hmax = extrema(@view wfc.mill[1, :])
    kmin, kmax = extrema(@view wfc.mill[2, :])
    lmin, lmax = extrema(@view wfc.mill[3, :])
    @assert hmax - hmin + 1 <= nrmesh[1] || is1 "'nrmesh' are too small."
    @assert kmax - kmin + 1 <= nrmesh[2] || is2 "'nrmesh' are too small."
    @assert lmax - lmin + 1 <= nrmesh[3] || is3 "'nrmesh' are too small."
    i1 = abs(-div(nrmesh[1], 2) - hmin)
    i2 = abs(-div(nrmesh[2], 2) - kmin)
    i3 = abs(-div(nrmesh[3], 2) - lmin)
    ns1 = (hmax - hmin) * is1 + nrmesh[1] * (!is1)
    ns2 = (kmax - kmin) * is2 + nrmesh[2] * (!is2)
    ns3 = (lmax - lmin) * is3 + nrmesh[3] * (!is3)
    ihs = i1 .+ (@view wfc.mill[1, :]) .- hmin .+ 1
    iks = i2 .+ (@view wfc.mill[2, :]) .- kmin .+ 1
    ils = i3 .+ (@view wfc.mill[3, :]) .- lmin .+ 1
    if is2d == true
        epn1 = exp(im*2π*n1)
        epn2 = exp(im*2π*n2)
        epn3 = exp(im*2π*n3)
    end

    cktmp = zeros(ComplexF64, (ns1, ns2, ns3, wfc.npol, wfc.nbnd))
    ∇cktmp = is∇u == true ? zeros(ComplexF64, (ns1, ns2, ns3, wfc.npol, wfc.nbnd, 3)) : Array{ComplexF64}(undef, 1,1,1,1,1,1)
    @inbounds for ipw in 1:wfc.igwx
        h, k, l = wfc.mill[:, ipw]
        if is2d == true
            epn = epn1^(is1 ? h : 0) * epn2^(is2 ? k : 0) * epn3^(is3 ? l : 0)
        else
            epn = 1.0
        end
        i1p, i2p, i3p = ihs[ipw], iks[ipw], ils[ipw]
        evcepn = wfc.evc[:, :, ipw] .* epn
        cktmp[i1p, i2p, i3p, :, :] .= evcepn
        if is∇u == true
            gvec = h * wfc.b1 + k * wfc.b2 + l * wfc.b3 
            for ii in 1:3
                ∇cktmp[i1p, i2p, i3p, :, :, ii] .= evcepn .* gvec[ii]
            end
        end
    end
    ck = FFTW.ifftshift(cktmp, 1:3)
    ∇ck = is∇u == true ? FFTW.ifftshift(∇cktmp, 1:3) : ∇cktmp
    return ck, ∇ck
end

function write_density(f0::Array{Float64, 3}; qedir::String="manual", savefile::String, atoms::Vector{String}=["none"], atomicpos::Matrix{Float64}=zeros(3,2), a1::Vector{Float64}=zeros(Float64, 3), a2::Vector{Float64}=zeros(Float64, 3), a3::Vector{Float64}=zeros(Float64, 3), comment::String="# Written on "*"$(Dates.now())", format="xsf")
    if format == "xsf"
        write_xsf(f0; qedir=qedir, savefile=savefile, atoms=atoms, atomicpos=atomicpos, a1=a1, a2=a2, a3=a3, comment=comment)
    elseif format == "grd"
        write_grd(f0; qedir=qedir, savefile=savefile, comment=comment)
    end
end

function write_xsf(f0::Array{Float64, 3}; qedir::String, savefile::String, atoms::Vector{String}, atomicpos::Matrix{Float64}, a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, comment::String)
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
        atoms = xml.atoms
        atomicpos = xml.atomicpos .* bohr2ang
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
    PF.@printf(io, "%2s\n", comment)
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
            PF.@printf(io, "%2s%10f%10f%10f\n", atoms[ia], atomicpos[1, ia], atomicpos[2, ia], atomicpos[3, ia])
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

angle_vec(a, b) = LA.atand(LA.norm(LA.cross(a,b)), LA.dot(a,b))

function write_grd(f0::Array{Float64, 3}; qedir::String, savefile::String, comment::String)
    xml = read_xml(qedir*"/data-file-schema.xml")
    a1ang = xml.a1.*bohr2ang
    a2ang = xml.a2.*bohr2ang
    a3ang = xml.a3.*bohr2ang
    α = round(angle_vec(xml.a2, xml.a3))
    β = round(angle_vec(xml.a3, xml.a1))
    γ = round(angle_vec(xml.a1, xml.a2))
    
    na1, na2, na3 = size(f0)
    io = open(savefile, "w")
    PF.@printf(io, "%2s\n", "# "*comment)
    PF.@printf(io, "%10f%12f%12f%15f%12f%12f\n", LA.norm(a1ang, 2), LA.norm(a2ang, 2), LA.norm(a3ang, 2), α, β, γ)
    PF.@printf(io, "%10d%10d%10d\n", na1, na2, na3)
    ia = 0
    for ia1 in 1:na1
        for ia2 in 1:na2
            for ia3 in 1:na3
                ia += 1
                if ia%6 == 0
                    PF.@printf(io, "%12f\n", f0[ia1, ia2, ia3])
                else
                    PF.@printf(io, "%12f", f0[ia1, ia2, ia3])
                end 
            end
        end
    end
    PF.@printf(io, "%2s\n", "")
    close(io)
end
