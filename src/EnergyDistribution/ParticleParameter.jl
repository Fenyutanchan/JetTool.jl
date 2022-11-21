transverse_momentum(a::AbstractParticle)::Real  =   norm(a.Momentum[1:2])

function pseudo_rapidity(a::AbstractParticle)::Real
    pa  =   norm(a.Momentum)
    pza =   last(a.Momentum)

    result  =   (1/2) * log(
        (pa + pza) / (pa - pza)
    )
    return  result
end