
stepc(x) = step(-x)

function methfessel_paxton_δ(x::Float64; n::Int=1)
    n ≥ 1 || error("n must be n ≥ 1")
    d = 1.0 / √(π)
    hold = 1.0
    hodd = 2 * x
    a = 1.0 / √(π)
    for i in 1:n
        a *= -1.0 / (i * 4.0)
        hnew = 2.0 * x * hodd - 2.0 * (2.0*i - 1.0) * hold
        d += a * hnew
        hodd = 2.0 * x * hnew - 2.0 * 2.0*i * hodd
        hold = hnew
    end
    d *= exp(-x^2)
    return d
end

function methfessel_paxton_step(x::Float64; n::Int=1)
    n ≥ 1 || error("n must be n ≥ 1")
    s0 = 0.5 * (1.0 - SF.erf(x))
    hold = 1.0
    hodd = 2 * x
    sn = 0.0 
    a = 1.0 / √(π)
    for i in 1:n
        a *= -1.0 / (i * 4.0)
        hnew = 2.0 * x * hodd - 2.0 * (2.0*i - 1.0) * hold
        sn += a * hodd
        hodd = 2.0 * x * hnew - 2.0 * 2.0*i * hodd
        hold = hnew
    end
    s = s0 + sn * exp(-x^2)
    return s
end

function calc_ikrange(nk::Int, nproc::Int, rank::Int)
    @assert 1 ≤ nproc
    @assert 0 ≤ rank < nproc
    @assert 1 ≤ nk

    q, r = divrem(nk, nproc)
    nlocal = q + (rank < r ? 1 : 0)
    ikst = rank * q + min(rank, r) + 1
    ikend = ikst + nlocal - 1
    return ikst, ikend
end
