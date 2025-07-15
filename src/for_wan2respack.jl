
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


function read_nkm(filename::String, nkirr::Int)
    ngi = zeros(Int, nkirr)
    io = open(filename, "r");
    for ik in 1:nkirr
        ngi[ik] = parse(Int, readline(io))
    end
    close(io)
    ntg = maximum(ngi) 
    return ngi, ntg
end

function make_kg0(ntg::Int, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, gcut::Float64, q::Vector{Float64}, ngl::Int)
    kg0 = zeros(Int, (3, ntg))
    qgL = zeros(3)
    igl = 0
    for igl1 in -ngl:ngl
        for igl2 in -ngl:ngl
            for igl3 in -ngl:ngl
                qgL = (q[1] + igl1)*b1 + (q[2] + igl2)*b2 + (q[3] + igl3)*b3    
                qgL2 = qgL[1]^2 + qgL[2]^2 + qgL[3]^2
                if qgL2 <= gcut
                    igl=igl+1
                    kg0[1,igl] = igl1
                    kg0[2,igl] = igl2
                    kg0[3,igl] = igl3 
                end
            end
        end
    end
    ng = igl
    return kg0, ng
end

function read_kg(filename::String, ngi::Vector{Int}, ntg::Int, nkirr::Int)
    io = open(filename, "r");
    kgi = zeros(Int, 3, ntg, nkirr)
    for ik in 1:nkirr
        ngiref = parse(Int, readline(io))
        @assert ngiref == ngi[ik]
        for ig in 1:ngiref
            kgi[:, ig, ik] = parse.(Int, split(readline(io)))
        end
    end
    close(io)
    return kgi
end

function calc_ecut_for_ψ(kgi::Array{Int, 3}, ntg::Int, nkirr::Int, ngi::Vector{Int}, ski::Matrix{Float64}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64})
    lkgi = zeros(ntg, nkirr)
    ecut = 0.0
    for ik in 1:nkirr
        for ig in 1:ngi[ik]
            ktmp = (ski[1,ik] + kgi[1,ig,ik])*b1 + (ski[2,ik] + kgi[2,ig,ik])*b2 + (ski[3,ik] + kgi[3,ig,ik])*b3
            lkgi[ig,ik] = ktmp[1]^2 + ktmp[2]^2 + ktmp[3]^2
        end
        ecut = maximum(lkgi) + 1e-8
    end
    return ecut
end

function calc_kg0s(gcut::Float64, rg::Array{Int, 3}, rw::Matrix{Int}, ngi::Vector{Int}, ntg::Int, trs::Vector{Int}, ski::Matrix{Float64}, b1::Vector{Float64}, b2::Vector{Float64}, b3::Vector{Float64}, ntk::Int, numirr::Vector{Int}, numrot::Vector{Int})
    bmin = minimum([LA.norm(b1, 2), LA.norm(b2, 2), LA.norm(b3, 2)])
    ngl = round(Int, √(gcut) / bmin) + 10
    ng0 = zeros(Int, ntk)
    kg0 = zeros(Int, (3, ntg, ntk))
    kgtmp = zeros(Int, (3,ntg))
    for jk in 1:ntk
        ik = numirr[jk]
        iop = numrot[jk]
        ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik] + rw[:,jk]
        kgtmp, ngψ = make_kg0(ntg,b1,b2,b3, gcut, ktmp, ngl)
        @assert ngψ == ngi[ik] "ERROR; ngψ should be ngI[ik]"
        ng0[jk] = ngψ
        if trs[jk]==1
            kg0[:,:,jk] = kgtmp
        elseif trs[jk] == -1
            kg0[:,:,jk] = -kgtmp # notice on - sign 
        end
    end
    return ng0, kg0
end


function kcheck(ktmp::Vector{Float64})
    rwtmp = zeros(Int, 3)
    δbz = 1.0e-6
    for i in 1:3
        if ktmp[i] > 1.50 + δbz
            ktmp[i] += - 2.0
            rwtmp[i] = -2
        end 
        if ktmp[i] > 0.50 + δbz
            ktmp[i] += - 1.0
            rwtmp[i] = -1
        end
        if ktmp[i] <= -1.50 + δbz
            ktmp[i] += 2.0
            rwtmp[i] = 2 
        end
        if ktmp[i] <= -0.5 + δbz
            ktmp[i] += 1.0
            rwtmp[i] = 1
        end
    end
    return ktmp, rwtmp
end

function kcheck_trs(ktmp::Vector{Float64})
    rwtmp = zeros(Int, 3)
    δbz = -1.0e-6
    for i in 1:3
        if ktmp[i] >= 1.50 + δbz
            ktmp[i] += -2.0
            rwtmp[i] = -2
        end 
        if ktmp[i] >= 0.50 + δbz 
            ktmp[i] += -1.0
            rwtmp[i] = -1
        end 
        if ktmp[i] < -1.50 + δbz 
            ktmp[i] = ktmp[1]+2.0
            rwtmp[i] = 2 
        end
        if ktmp[i] < -0.50 + δbz 
            ktmp[i] += 1.0
            rwtmp[i] = 1
        end
    end
    return ktmp, rwtmp
end


function est_ntk(nkirr::Int, nsymq::Int, ski::Matrix{Float64}, rg::Array{Int, 3})
    n=nkirr*nsymq*2
    # sk0
    sk0 = zeros(3, n)
    jk=0
    for ik in 1:nkirr
        for iop in 1:nsymq
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, rwtmp = kcheck(ktmp) # rewind check 
            flag1000 = false
            for iik in 1:jk
                if abs(sk0[1,iik] - ktmp[1]) < 1e-4 && abs(sk0[2,iik] - ktmp[2])<1e-4 && abs(sk0[3,iik] - ktmp[3]) < 1e-4
                    flag1000 = true
                    break
                end
            end
            if flag1000 == false
                jk = jk + 1
                sk0[:,jk] = -ktmp
            end
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, rwtmp = kcheck_trs(ktmp) # rewind check 
            flag2000 = false
            for iik in 1:jk
                if maximum(abs, sk0[:,iik] - ktmp) < 1e-4
                    flag2000 = true
                    break
                end
            end
            if flag2000 == false
                jk = jk + 1
                sk0[:,jk] = ktmp
            end
        end
    end
    ntk = jk
    if ntk > n
        @assert false "Estimated ntk is too large; stop"
    end
    return sk0, ntk
end

function est_nkbi(n::Int, sk::Matrix{Float64}) 
    nkb = zeros(Int, 3)
    for b in 1:3
        x = 1.0
        for i in 1:n
            if abs(sk[b,i]) < 1e-7
                continue
            end
            if abs(sk[b,i]) < x
                x = abs(sk[b,i])  
            end
        end
        nkb[b] = round(Int, 1.0 / x)
    end
    ntk = nkb[1] * nkb[2] * nkb[3]
    return ntk
end 

function read_sample_k(filename::String, nsymq::Int, rg::Array{Int, 3})
    io = open(filename, "r");
    nkirr = parse(Int, readline(io))
    ski = zeros(3, nkirr)
    for ik in 1:nkirr
        ski[:, ik] = parse.(Float64, split(readline(io)))
    end
    close(io)

    sk0, ntk = est_ntk(nkirr,nsymq,ski,rg)
    #--
    # sk0,numirr,numrot,trs,rw
    sk0 = zeros(3,ntk)
    numirr = zeros(Int, ntk)
    numrot = zeros(Int, ntk)
    trs = zeros(Int, ntk)
    rw = zeros(Int, (3,ntk))
    nummk = zeros(Int, nkirr)
    jk = 0
    for ik in 1:nkirr
        initial_flg = 0
        for iop in 1:nsymq
            # sym
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, rwtmp = kcheck(ktmp) # rewind check
            flag1000 = false
            for iik in 1:jk
                if maximum(abs, sk0[:,iik] - ktmp) < 1e-4
                    flag1000 = true
                    break
                end
            end
            if flag1000 == false
                jk += 1
                sk0[:,jk] = ktmp
                numirr[jk] = ik
                numrot[jk] = iop
                trs[jk] = 1
                rw[:,jk] = rwtmp
                if initial_flg == 0
                    nummk[ik] = jk
                    initial_flg = 1
                end
            end
            # time-reversal
            ktmp = rg[:,1,iop]*ski[1,ik] + rg[:,2,iop]*ski[2,ik] + rg[:,3,iop]*ski[3,ik]
            ktmp, rwtmp = kcheck_trs(ktmp) # rewind check
            flag2000 = false
            for iik = 1:jk
                if maximum(abs, sk0[:,iik] + ktmp) < 1e-4
                    flag2000 = true
                    break
                end
            end
            if flag2000 == false
                jk += 1
                sk0[:,jk] = -ktmp 
                numirr[jk] = ik
                numrot[jk] = iop
                trs[jk] = -1
                rw[:,jk] = rwtmp 
            end
        end
    end
    ntk = est_nkbi(ntk, sk0)
    @assert ntk == jk "ERROR; ntk(=$(ntk)) should be jk(=$(jk))"
    return ski, sk0, numirr, numrot, trs, rw, nkirr, ntk
end





function read_lattice(filename::String)
    io = open(filename, "r");
    a1 = parse.(Float64, split(readline(io)))
    a2 = parse.(Float64, split(readline(io)))
    a3 = parse.(Float64, split(readline(io)))
    close(io)
    return a1, a2, a3
end

function read_wan(io, nwfc::Int, npol::Int, ng::Int)
    cs = zeros(ComplexF64, (ng, npol, nwfc))
    seek(io, 8)
    for is in 1:npol
        for ib in 1:nwfc
            for ig in 1:ng
                cs[ig, is, ib] = read(io, ComplexF64)
            end
        end
    end
    return cs
end
