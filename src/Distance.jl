calc_distance(a::AbstractParticle, b::AbstractParticle, p::Real)::Real   =   min(
    transverse_momentum(a)^(2 * p), transverse_momentum(b)^(2 * p)
) * calc_ΔR(a, b)^2

calc_ΔR(a::AbstractParticle, b::AbstractParticle)::Real   =   sqrt(
    calc_Δy(a, b)^2 + calc_Δφ(a, b)^2
)

calc_Δy(a::AbstractParticle, b::AbstractParticle)::Real =   abs(
    pseudo_rapidity(a) - pseudo_rapidity(b)
)

function calc_Δφ(a::AbstractParticle, b::AbstractParticle)::Real
    pT_inner_product    =   transpose(a.Momentum[1:2]) * b.Momentum[1:2]
    cosΔφ               =   pT_inner_product / (
        transverse_momentum(a) * transverse_momentum(b)
    )
    if abs(cosΔφ) > 1
        return  acos(cosΔφ / abs(cosΔφ))
    else
        return  acos(cosΔφ)
    end
end