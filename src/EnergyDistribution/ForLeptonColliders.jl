function energy_within_cone_r(j::Jet, r::Real)::Real
    result  =   0
    for particle ∈ j.Particles
        Δr  =   calc_Δθ(j, particle)

        if Δr ≤ r
            result  +=  particle.Energy
        end
    end
    return  result
end

function calc_Δθ(a::AbstractParticle, b::AbstractParticle)
    cosΔθ   =   transpose(a.Momentum) * b.Momentum / (
        norm(a.Momentum) * norm(b.Momentum)
    )
    cosΔθ   =   if cosΔθ > 1
        printstyled(
            "Warning: $cosΔθ is larger than 1! We have set it to 1. Please check it!\n",
            color=:red
        )
        1
    elseif cosΔθ < -1
        printstyled(
            "Warning: $cosΔθ is less than -1! We have set it to -1. Please check it!\n",
            color=:red
        )
        -1
    else
        cosΔθ
    end
    
    return  acos(cosΔθ)
end