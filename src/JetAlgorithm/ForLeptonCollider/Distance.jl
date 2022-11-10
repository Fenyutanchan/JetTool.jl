calc_cosΔθ(
    a::AbstractParticle,
    b::AbstractParticle
)::Real = transpose(a.Momentum) * b.Momentum / (
    abs(a.Momentum) * abs(b.Momentum)
)

distance_Spherical_Generalized_kT(
    a::AbstractParticle,
    b::AbstractParticle,
    R::Real;
    p::Real=-1
)   =   min(
    a.Energy^(2p),
    b.Energy^(2p)
) * (1 - calc_cosΔθ(a, b)) / (1 - cos(R))
distance_Spherical_Generalized_kT(
    a::AbstractParticle;
    p::Real=-1
)   =   a.Energy^(2p)

distance_Spherical_kT(
    a::AbstractParticle,
    b::AbstractParticle
)   =   2 * min(a.Energy^2, b.Energy^2) * (1 - calc_cosΔθ(a, b))