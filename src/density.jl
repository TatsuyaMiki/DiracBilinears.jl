
function calc_density(;calc::String, qedir::String, n1::Float64=0.0, n2::Float64=0.0, n3::Float64=0.0, nrmesh::Tuple=(0,0,0), δμ::Float64=0.0, emin::Float64=100000.0, smearing::String="step", degauss::Float64=0.01, ikst::Int=1, ikend::Union{Int,Nothing}=nothing)
    ## - n1, n2, n3 are parameters which can be used for calculations in a 2D plane
    ##   The calculations are performed on a plane perpendicular to ai (i=1,2,3) that passes through ni*ai
    ## - The chemical potential can be shifted using δμ (eV)
    ## - The core levels can be excluded by using emin (eV)
    ## - ikst and ikend specify the first and last k-point indices
    is∇u = is∇ukn(calc)
    xml = read_xml(qedir*"/data-file-schema.xml")
    ikend_ = isnothing(ikend) ? xml.nxk : ikend
    1 ≤ ikst ≤ ikend_ ≤ xml.nxk || error("k-point indices must satisfy 1 ≤ ikst ≤ ikend ≤ $(xml.nxk)")
    nrmesh_ = make_nrmesh(xml=xml, nrmesh=nrmesh)
    o = make_zeros_density(calc, nrmesh_)
    volume = xml.nxk * abs(LA.dot(xml.a3, LA.cross(xml.a1, xml.a2)))
    for ik in ikst:ikend_
        wfc = qewfc(ik, xml, qedir)
        occ = calc_occupation(xml.e[:, ik]; ef=xml.ef, smearing=smearing, degauss=degauss, δμ=δμ, emin=emin)
        ukn, ∇ukn = make_c_k(nrmesh_, wfc; n1=n1, n2=n2, n3=n3, is∇u=is∇u)
        calc_fourier_k!(ukn; nrmesh=nrmesh_)
        ukn ./= √(volume)
        is∇u == true ? calc_fourier_k!(∇ukn; nrmesh=nrmesh_) : ∇ukn
        ∇ukn ./= √(volume)
        calc_density_ok!(o; calc=calc, nrmesh=nrmesh_, wfc=wfc, ukn=ukn, ∇ukn=∇ukn, occ=occ)
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
        evc = zeros(ComplexF64, npol, nbnd, igwx)
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
        error("Invalid value assigned to 'smearing'")
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
        return zeros(Float64, nrmesh[1], nrmesh[2], nrmesh[3])
    elseif calc == "ms"
        return zeros(Float64, 3, nrmesh[1], nrmesh[2], nrmesh[3])
    elseif calc == "j"
        return zeros(Float64, 3, nrmesh[1], nrmesh[2], nrmesh[3])
    elseif calc == "∇ρ" || calc == "nabla_rho"
        return zeros(Float64, 3, nrmesh[1], nrmesh[2], nrmesh[3])
    elseif calc == "∇ms" || calc == "nabla_ms"
        return zeros(Float64, nrmesh[1], nrmesh[2], nrmesh[3])
    elseif calc == "τz" || calc == "tau_z" || calc == "chirality"
        return zeros(Float64, nrmesh[1], nrmesh[2], nrmesh[3])
    elseif calc == "ps"
        return zeros(Float64, 3, nrmesh[1], nrmesh[2], nrmesh[3])
    else
        error("Invalid value assigned to 'calc'")
    end
end


function calc_density_ok!(out; calc::String, nrmesh::Tuple, wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    if calc == "ρ" || calc == "rho"
        calc_density_ρ!(out; ukn=ukn, occ=occ)
        return out
    elseif calc == "ms"
        calc_density_ms!(out; ukn=ukn, occ=occ)
        return out
    elseif calc == "j"
        calc_density_j!(out; wfc=wfc, ukn=ukn, ∇ukn=∇ukn, occ=occ)
        return out
    elseif calc == "∇ρ" || calc == "nabla_rho"
        calc_density_∇ρ!(out; ukn=ukn, ∇ukn=∇ukn, occ=occ)
        return out
    elseif calc == "∇ms" || calc == "nabla_ms"
        calc_density_∇ms!(out; wfc=wfc, ukn=ukn, ∇ukn=∇ukn, occ=occ)
        return out
    elseif calc == "τz" || calc == "tau_z" || calc == "chirality"
        calc_density_τz!(out; wfc=wfc, ukn=ukn, ∇ukn=∇ukn, occ=occ)
        return out
    elseif calc == "ps"
        calc_density_ps!(out; wfc=wfc, ukn=ukn, ∇ukn=∇ukn, occ=occ)
        return out
    else
        error("Invalid value assigned to 'calc'")
    end
end

function calc_density_ρ!(out; ukn, occ::Vector{Float64})
    nx, ny, nz, ns, nb = size(ukn)
    @inbounds for ib in 1:nb
        occib = occ[ib]
        occib == 0.0 && continue
        for is in 1:ns, iz in 1:nz, iy in 1:ny, ix in 1:nx
            out[ix, iy, iz] += abs2(ukn[ix, iy, iz, is, ib]) * occib
        end
    end
    return out
end

function calc_density_ms!(out; ukn, occ::Vector{Float64})
    nx, ny, nz, ns, nb = size(ukn)
    @inbounds for ib in 1:nb
        occib = occ[ib]
        occib == 0.0 && continue
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            u1 = ukn[ix, iy, iz, 1, ib]
            u2 = ukn[ix, iy, iz, 2, ib]
            out[1, ix, iy, iz] -= 2.0 * real(conj(u1) * u2) * occib
            out[2, ix, iy, iz] -= 2.0 * imag(conj(u1) * u2) * occib
            out[3, ix, iy, iz] -= (abs2(u1) - abs2(u2)) * occib
        end
    end
    return out
end


function calc_density_∇ms(wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    @inbounds for ib in 1:wfc.nbnd
        ukn[:, :, :, :, ib] .*= occ[ib]
    end
    return -2.0*imag.(ES.optein"xyzsbi,sti,xyztb->xyz"(conj.(∇ukn), σ, ukn))
end

function calc_density_j!(out; wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    nx, ny, nz, ns, nb = size(ukn)
    @inbounds for ib in 1:nb
        occib = occ[ib]
        occib == 0.0 && continue
        for is in 1:ns, iz in 1:nz, iy in 1:ny, ix in 1:nx
            u = ukn[ix, iy, iz, is, ib]
            uc = conj(u) * occib
            out[1, ix, iy, iz] += 2.0 * real(uc * (wfc.xk[1] * u + ∇ukn[ix, iy, iz, is, ib, 1]))
            out[2, ix, iy, iz] += 2.0 * real(uc * (wfc.xk[2] * u + ∇ukn[ix, iy, iz, is, ib, 2]))
            out[3, ix, iy, iz] += 2.0 * real(uc * (wfc.xk[3] * u + ∇ukn[ix, iy, iz, is, ib, 3]))
        end
    end
    return out
end

function calc_density_∇ρ!(out; ukn, ∇ukn, occ::Vector{Float64})
    nx, ny, nz, ns, nb = size(ukn)
    @inbounds for ib in 1:nb
        occib = occ[ib]
        occib == 0.0 && continue
        for is in 1:ns, iz in 1:nz, iy in 1:ny, ix in 1:nx
            u = ukn[ix, iy, iz, is, ib]
            out[1, ix, iy, iz] += 2.0 * imag(conj(∇ukn[ix, iy, iz, is, ib, 1]) * u) * occib
            out[2, ix, iy, iz] += 2.0 * imag(conj(∇ukn[ix, iy, iz, is, ib, 2]) * u) * occib
            out[3, ix, iy, iz] += 2.0 * imag(conj(∇ukn[ix, iy, iz, is, ib, 3]) * u) * occib
        end
    end
    return out
end

function calc_density_∇ms!(out; wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    wfc.npol == 2 || error("npol must be 2")
    nx, ny, nz, ns, nb = size(ukn)
    ns == 2 || error("ns must be 2")
    @inbounds for ib in 1:nb
        occib = occ[ib]
        occib == 0.0 && continue
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            u1 = ukn[ix, iy, iz, 1, ib]
            u2 = ukn[ix, iy, iz, 2, ib]
            dux1 = ∇ukn[ix, iy, iz, 1, ib, 1]
            dux2 = ∇ukn[ix, iy, iz, 2, ib, 1]
            duy1 = ∇ukn[ix, iy, iz, 1, ib, 2]
            duy2 = ∇ukn[ix, iy, iz, 2, ib, 2]
            duz1 = ∇ukn[ix, iy, iz, 1, ib, 3]
            duz2 = ∇ukn[ix, iy, iz, 2, ib, 3]
            val =
                conj(dux1) * u2 + conj(dux2) * u1 +
                (-im) * conj(duy1) * u2 + (im) * conj(duy2) * u1 +
                conj(duz1) * u1 - conj(duz2) * u2
            out[ix, iy, iz] += -2.0 * imag(val) * occib
        end
    end
    return out
end


function calc_density_τz!(out; wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    wfc.npol == 2 || error("npol must be 2")
    nx, ny, nz, ns, nb = size(ukn)
    ns == 2 || error("ns must be 2")
    @inbounds for ib in 1:nb
        occib = occ[ib]
        occib == 0.0 && continue
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            u1 = ukn[ix, iy, iz, 1, ib]
            u2 = ukn[ix, iy, iz, 2, ib]
            c1 = conj(u1) * occib
            c2 = conj(u2) * occib
            k1u1 = wfc.xk[1] * u1 + ∇ukn[ix, iy, iz, 1, ib, 1]
            k1u2 = wfc.xk[1] * u2 + ∇ukn[ix, iy, iz, 2, ib, 1]
            k2u1 = wfc.xk[2] * u1 + ∇ukn[ix, iy, iz, 1, ib, 2]
            k2u2 = wfc.xk[2] * u2 + ∇ukn[ix, iy, iz, 2, ib, 2]
            k3u1 = wfc.xk[3] * u1 + ∇ukn[ix, iy, iz, 1, ib, 3]
            k3u2 = wfc.xk[3] * u2 + ∇ukn[ix, iy, iz, 2, ib, 3]
            val =
                c1 * k1u2 + c2 * k1u1 +
                (-im) * c1 * k2u2 + (im) * c2 * k2u1 +
                c1 * k3u1 - c2 * k3u2
            out[ix, iy, iz] += 2.0 * real(val)
        end
    end
    return out
end

function calc_density_ps!(out; wfc::Wfc, ukn, ∇ukn, occ::Vector{Float64})
    wfc.npol == 2 || error("npol must be 2")
    nx, ny, nz, ns, nb = size(ukn)
    ns == 2 || error("ns must be 2")
    @inbounds for ib in 1:nb
        occib = occ[ib]
        occib == 0.0 && continue
        for iz in 1:nz, iy in 1:ny, ix in 1:nx
            u1 = ukn[ix, iy, iz, 1, ib]
            u2 = ukn[ix, iy, iz, 2, ib]
            c1 = conj(u1) * occib
            c2 = conj(u2) * occib

            k1u1 = wfc.xk[1] * u1 + ∇ukn[ix, iy, iz, 1, ib, 1]
            k1u2 = wfc.xk[1] * u2 + ∇ukn[ix, iy, iz, 2, ib, 1]
            k2u1 = wfc.xk[2] * u1 + ∇ukn[ix, iy, iz, 1, ib, 2]
            k2u2 = wfc.xk[2] * u2 + ∇ukn[ix, iy, iz, 2, ib, 2]
            k3u1 = wfc.xk[3] * u1 + ∇ukn[ix, iy, iz, 1, ib, 3]
            k3u2 = wfc.xk[3] * u2 + ∇ukn[ix, iy, iz, 2, ib, 3]

            # tmp_jk = sum_{s,t} conj(u_s) * (k_j u)_t * sigma_k[s,t]
            tmp12 = (-im) * c1 * k1u2 + (im) * c2 * k1u1
            tmp13 = c1 * k1u1 - c2 * k1u2
            tmp21 = c1 * k2u2 + c2 * k2u1
            tmp23 = c1 * k2u1 - c2 * k2u2
            tmp31 = c1 * k3u2 + c2 * k3u1
            tmp32 = (-im) * c1 * k3u2 + (im) * c2 * k3u1
            out[1, ix, iy, iz] += -2.0 * real(tmp23 - tmp32)
            out[2, ix, iy, iz] += -2.0 * real(tmp31 - tmp13)
            out[3, ix, iy, iz] += -2.0 * real(tmp12 - tmp21)
        end
    end
    return out
end

function calc_fourier_k!(ukn; nrmesh::Tuple)
    if nrmesh[1] == 1
        ukn .= sum(ukn, dims=1)
        FFTW.bfft!(ukn, [2,3])
        return ukn
    elseif nrmesh[2] == 1
        ukn .= sum(ukn, dims=2)
        FFTW.bfft!(ukn, [1,3])
        return ukn
    elseif nrmesh[3] == 1
        ukn .= sum(ukn, dims=3)
        FFTW.bfft!(ukn, [1,2])
        return ukn
    else
        FFTW.bfft!(ukn, [1,2,3])
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
    (hmax - hmin + 1 ≤ nrmesh[1] || is1) || error("nrmesh are too small")
    (kmax - kmin + 1 ≤ nrmesh[2] || is2) || error("nrmesh are too small")
    (lmax - lmin + 1 ≤ nrmesh[3] || is3) || error("nrmesh are too small")
    ns1 = is1 ? hmax - hmin + 1 : nrmesh[1]
    ns2 = is2 ? kmax - kmin + 1 : nrmesh[2]
    ns3 = is3 ? lmax - lmin + 1 : nrmesh[3]
    if is2d == true
        epn1 = exp(im*2π*n1)
        epn2 = exp(im*2π*n2)
        epn3 = exp(im*2π*n3)
    end

    ck = zeros(ComplexF64, ns1, ns2, ns3, wfc.npol, wfc.nbnd)
    ∇ck = is∇u == true ? zeros(ComplexF64, ns1, ns2, ns3, wfc.npol, wfc.nbnd, 3) : Array{ComplexF64}(undef, 1,1,1,1,1,1)
    @inbounds for ipw in 1:wfc.igwx
        h = wfc.mill[1, ipw]
        k = wfc.mill[2, ipw]
        l = wfc.mill[3, ipw]
        if is2d == true
            epn = epn1^(is1 ? h : 0) * epn2^(is2 ? k : 0) * epn3^(is3 ? l : 0)
        else
            epn = 1.0
        end
        i1p = is1 ? (h - hmin + 1) : (mod(h, nrmesh[1]) + 1)
        i2p = is2 ? (k - kmin + 1) : (mod(k, nrmesh[2]) + 1)
        i3p = is3 ? (l - lmin + 1) : (mod(l, nrmesh[3]) + 1)
        gx = h * wfc.b1[1] + k * wfc.b2[1] + l * wfc.b3[1]
        gy = h * wfc.b1[2] + k * wfc.b2[2] + l * wfc.b3[2]
        gz = h * wfc.b1[3] + k * wfc.b2[3] + l * wfc.b3[3]
        @inbounds for ib in 1:wfc.nbnd, is in 1:wfc.npol
            v = wfc.evc[is, ib, ipw] * epn
            ck[i1p, i2p, i3p, is, ib] = v
            if is∇u
                ∇ck[i1p, i2p, i3p, is, ib, 1] = v * gx
                ∇ck[i1p, i2p, i3p, is, ib, 2] = v * gy
                ∇ck[i1p, i2p, i3p, is, ib, 3] = v * gz
            end
        end
    end
    return ck, ∇ck
end

function write_density(f0::Array{Float64, 3}; qedir::String="manual", savefile::String, atoms::Vector{String}=["none"], atomicpos::Matrix{Float64}=zeros(3,2), a1::Vector{Float64}=zeros(Float64, 3), a2::Vector{Float64}=zeros(Float64, 3), a3::Vector{Float64}=zeros(Float64, 3), comment::String="# Written on "*"$(Dates.now())", format::String="xsf", unit::String="bohr")
    if qedir != "manual"
        a1 == zeros(Float64, 3) && a2 == zeros(Float64, 3) && a3 == zeros(Float64, 3) || error("a1, a2, a3 can only be specified when qedir='manual'")
        xml = read_xml(qedir*"/data-file-schema.xml")
        a1 = xml.a1
        a2 = xml.a2
        a3 = xml.a3
        atoms = xml.atoms
        atomicpos = xml.atomicpos
    end
    f0, a1, a2, a3, atomicpos = convert_units(f0, a1, a2, a3, atomicpos; unit=unit)
    if format == "xsf"
        write_xsf(f0; savefile=savefile, a1=a1, a2=a2, a3=a3, atoms=atoms, atomicpos=atomicpos, comment=comment)
    elseif format == "grd"
        write_grd(f0; savefile=savefile, a1=a1, a2=a2, a3=a3, comment=comment)
    end
end

function convert_units(f0::Array{Float64, 3}, a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, atomicpos::Matrix{Float64}; unit::String)
    if unit=="ang"
        f0 ./= bohr2ang^3
        a1 .*= bohr2ang
        a2 .*= bohr2ang
        a3 .*= bohr2ang
        atomicpos .*= bohr2ang
    end
    return f0, a1, a2, a3, atomicpos
end

function write_xsf(f0::Array{Float64, 3}; savefile::String, a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, atoms::Vector{String}, atomicpos::Matrix{Float64}, comment::String)
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
    PF.@printf(io, "%15f%10f%10f\n", a1[1], a1[2], a1[3])
    PF.@printf(io, "%15f%10f%10f\n", a2[1], a2[2], a2[3])
    PF.@printf(io, "%15f%10f%10f\n", a3[1], a3[2], a3[3])
    PF.@printf(io, "%2s\n", "CONVVEC")
    PF.@printf(io, "%15f%10f%10f\n", a1[1], a1[2], a1[3])
    PF.@printf(io, "%15f%10f%10f\n", a2[1], a2[2], a2[3])
    PF.@printf(io, "%15f%10f%10f\n", a3[1], a3[2], a3[3])
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
    PF.@printf(io, "%15f%10f%10f\n", a1[1], a1[2], a1[3])
    PF.@printf(io, "%15f%10f%10f\n", a2[1], a2[2], a2[3])
    PF.@printf(io, "%15f%10f%10f\n", a3[1], a3[2], a3[3])
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

function write_grd(f0::Array{Float64, 3}; savefile::String, a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, comment::String)
    α = round(angle_vec(a2, a3))
    β = round(angle_vec(a3, a1))
    γ = round(angle_vec(a1, a2))
    
    na1, na2, na3 = size(f0)
    io = open(savefile, "w")
    PF.@printf(io, "%2s\n", "# "*comment)
    PF.@printf(io, "%10f%12f%12f%15f%12f%12f\n", LA.norm(a1, 2), LA.norm(a2, 2), LA.norm(a3, 2), α, β, γ)
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
