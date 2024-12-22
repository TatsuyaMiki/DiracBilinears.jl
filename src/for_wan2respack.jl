
function read_symmetry(filename::String)
    io = open(filename, "r");
    nsymq = parse(Int, readline(io))
    nnp = parse(Int, readline(io))
    rg = zeros(Int, (3, 3, nsymq))
    pg = zeros(Int, (3, nsymq))
    rginv = zeros(Int, (3, 3, nsymq))
    for iop in 1:nsymq
        rgtmp = parse.(Int, split(readline(io)))
        pg = parse.(Int, split(readline(io)))
        it = 0
        for j in 1:3
            for i in 1:3
                it += 1
                rg[i, j, iop] = rgtmp[it]
            end
        end
    end
    close(io)
    rginv = zeros(Int, (3, 3, nsymq))
    for iop in 1:nsymq
        rginv[:, :, iop] = inv(rg[:, :, iop])
    end
    return rg, rginv, nsymq, nnp
end


function read_nkm(filename::String, Nk_irr::Int)
    ngi = zeros(Int, Nk_irr)
    io = open(filename, "r");
    for ik in 1:Nk_irr
        ngi[ik] = parse(Int, readline(io))
    end
    close(io)
    ntg = maximum(ngi) 
    return ngi, ntg
end

function make_kg0(NTG::Int, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, Gcut::Float64, q::Vector{Float64}, NGL::Int)
    KG0 = zeros(Int, (3, NTG))
    qgL = zeros(3)
    igL = 0
    for igL1 in -NGL:NGL
        for igL2 in -NGL:NGL
            for igL3 in -NGL:NGL
                qgL = (q[1] + igL1)*b1 + (q[2] + igL2)*b2 + (q[3] + igL3)*b3    
                qgL2 = qgL[1]^2 + qgL[2]^2 + qgL[3]^2
                if qgL2 <= Gcut
                    igL=igL+1
                    KG0[1,igL] = igL1
                    KG0[2,igL] = igL2
                    KG0[3,igL] = igL3 
                end
            end
        end
    end
    NG = igL
    return KG0, NG
end

function read_kg(filename::String, ngi::Vector{Int}, ntg::Int, Nk_irr::Int)
    io = open(filename, "r");
    kgi = zeros(Int, 3, ntg, Nk_irr)
    for ik in 1:Nk_irr
        ngiref = parse(Int, readline(io))
        @assert ngiref == ngi[ik]
        for ig in 1:ngiref
            kgi[:, ig, ik] = parse.(Int, split(readline(io)))
        end
    end
    close(io)
    return kgi
end

function calc_ecut_for_ψ(kgi::Array{Int, 3}, ntg::Int, Nk_irr::Int, ngi::Vector{Int}, ski::Matrix{Float64}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64})
    LKGI = zeros(ntg, Nk_irr)
    ecut = 0.0
    for ik in 1:Nk_irr
        for ig in 1:ngi[ik]
            ktmp = (ski[1,ik] + kgi[1,ig,ik])*b1 + (ski[2,ik] + kgi[2,ig,ik])*b2 + (ski[3,ik] + kgi[3,ig,ik])*b3
            LKGI[ig,ik] = ktmp[1]^2 + ktmp[2]^2 + ktmp[3]^2
        end
        ecut = maximum(LKGI) + 1e-8
    end
    return ecut
end

function calc_kg0s(Gcut::Float64, rg::Array{Int, 3}, RW::Matrix{Int}, ngi::Vector{Int}, ntg::Int, trs::Vector{Int}, ski::Matrix{Float64}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, NTK::Int, numirr::Vector{Int}, numrot::Vector{Int})
    bmin = minimum([LA.norm(b1, 2), LA.norm(b2, 2), LA.norm(b3, 2)])
    NGL = round(Int, √(Gcut) / bmin) + 10
    ng0 = zeros(Int, NTK)
    kg0 = zeros(Int, (3, ntg, NTK))
    KGtmp = zeros(Int, (3,ntg))
    for jk in 1:NTK
        ik = numirr[jk]
        iop = numrot[jk]
        ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik] + RW[:,jk]
        KGtmp, NG_for_psi = make_kg0(ntg,b1,b2,b3, Gcut, ktmp, NGL)
        @assert NG_for_psi == ngi[ik] "ERROR; NG_for_psi should be NGI[ik]"
        ng0[jk] = NG_for_psi
        if trs[jk]==1
            kg0[:,:,jk] = KGtmp
        elseif trs[jk] == -1
            kg0[:,:,jk] = -KGtmp # notice on - sign 
        end
    end
    return ng0, kg0
end


function kcheck(ktmp::Vector{Float64})
    RWtmp = zeros(Int, 3)
    δbz = 1.0e-6
    for i in 1:3
        if ktmp[i] > 1.50 + δbz
            ktmp[i] += - 2.0
            RWtmp[i] = -2
        end 
        if ktmp[i] > 0.50 + δbz
            ktmp[i] += - 1.0
            RWtmp[i] = -1
        end
        if ktmp[i] <= -1.50 + δbz
            ktmp[i] += 2.0
            RWtmp[i] = 2 
        end
        if ktmp[i] <= -0.5 + δbz
            ktmp[i] += 1.0
            RWtmp[i] = 1
        end
    end
    return ktmp, RWtmp
end

function kcheck_trs(ktmp::Vector{Float64})
    RWtmp = zeros(Int, 3)
    δbz = -1.0e-6
    for i in 1:3
        if ktmp[i] >= 1.50 + δbz
            ktmp[i] += -2.0
            RWtmp[i] = -2
        end 
        if ktmp[i] >= 0.50 + δbz 
            ktmp[i] += -1.0
            RWtmp[i] = -1
        end 
        if ktmp[i] < -1.50 + δbz 
            ktmp[i] = ktmp[1]+2.0
            RWtmp[i] = 2 
        end
        if ktmp[i] < -0.50 + δbz 
            ktmp[i] += 1.0
            RWtmp[i] = 1
        end
    end
    return ktmp, RWtmp
end


function est_NTK(Nk_irr::Int, Nsymq::Int, ski::Matrix{Float64}, rg::Array{Int, 3})
    N=Nk_irr*Nsymq*2
    # SK0
    SK0 = zeros(3, N)
    jk=0
    for ik in 1:Nk_irr
        for iop in 1:Nsymq
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, RWtmp = kcheck(ktmp) # rewind check 
            flag1000 = false
            for iik in 1:jk
                if abs(SK0[1,iik] - ktmp[1]) < 1e-4 && abs(SK0[2,iik] - ktmp[2])<1e-4 && abs(SK0[3,iik] - ktmp[3]) < 1e-4
                    flag1000 = true
                    break
                end
            end
            if flag1000 == false
                jk = jk + 1
                SK0[:,jk] = -ktmp
            end
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, RWtmp = kcheck_trs(ktmp) # rewind check 
            flag2000 = false
            for iik in 1:jk
                if maximum(abs, SK0[:,iik] - ktmp) < 1e-4
                    flag2000 = true
                    break
                end
            end
            if flag2000 == false
                jk = jk + 1
                SK0[:,jk] = ktmp
            end
        end
    end
    NTK = jk
    if NTK > N
        @assert false "Estimated NTK is too large; stop"
    end
    return SK0, NTK
end

function est_nkbi(N::Int, SK::Matrix{Float64}) 
    nkb = zeros(Int, 3)
    for b in 1:3
        x = 1.0
        for i in 1:N
            if abs(SK[b,i]) < 1e-7
                continue
            end
            if abs(SK[b,i]) < x
                x = abs(SK[b,i])  
            end
        end
        nkb[b] = round(Int, 1.0 / x)
    end
    NTK = nkb[1] * nkb[2] * nkb[3]
    return NTK
end 

function read_sample_k(filename::String, nsymq::Int, rg::Array{Int, 3})
    io = open(filename, "r");
    Nk_irr = parse(Int, readline(io))
    ski = zeros(3, Nk_irr)
    for ik in 1:Nk_irr
        ski[:, ik] = parse.(Float64, split(readline(io)))
    end
    close(io)

    SK0, NTK = est_NTK(Nk_irr,nsymq,ski,rg)
    #--
    # SK0,numirr,numrot,trs,RW
    SK0 = zeros(3,NTK)
    numirr = zeros(Int, NTK)
    numrot = zeros(Int, NTK)
    trs = zeros(Int, NTK)
    RW = zeros(Int, (3,NTK))
    numMK = zeros(Int, Nk_irr)
    jk = 0
    for ik in 1:Nk_irr
        initial_flg = 0
        for iop in 1:nsymq
            # sym
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, RWtmp = kcheck(ktmp) # rewind check
            flag1000 = false
            for iik in 1:jk
                if maximum(abs, SK0[:,iik] - ktmp) < 1e-4
                    flag1000 = true
                    break
                end
            end
            if flag1000 == false
                jk += 1
                SK0[:,jk] = ktmp
                numirr[jk] = ik
                numrot[jk] = iop
                trs[jk] = 1
                RW[:,jk] = RWtmp
                if initial_flg == 0
                    numMK[ik] = jk
                    initial_flg = 1
                end
            end
            # time-reversal
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, RWtmp = kcheck_trs(ktmp) # rewind check
            flag2000 = false
            for iik = 1:jk
                if maximum(abs, SK0[:,iik] + ktmp) < 1e-4
                    flag2000 = true
                    break
                end
            end
            if flag2000 == false
                jk += 1
                SK0[:,jk] = -ktmp 
                numirr[jk] = ik
                numrot[jk] = iop
                trs[jk] = -1
                RW[:,jk] = RWtmp 
            end
        end
    end
    NTK = est_nkbi(NTK, SK0)
    @assert NTK == jk "ERROR; NTK(=$(NTK)) should be jk(=$(jk))"
    return ski, SK0, numirr, numrot, trs, RW, Nk_irr, NTK
end





function read_lattice(filename)
    io = open(filename, "r");
    a1 = parse.(Float64, split(readline(io)))
    a2 = parse.(Float64, split(readline(io)))
    a3 = parse.(Float64, split(readline(io)))
    close(io)
    return a1, a2, a3
end

function read_wan(io, nwfc::Int, npol::Int, ng::Int)
    cs = zeros(ComplexF64, (ng, npol, nwfc))
    skip(io, 8)
    for is in 1:npol
        for ib in 1:nwfc
            for ig in 1:ng
                cs[ig, is, ib] = read(io, ComplexF64)
            end
        end
    end
    return cs
end
