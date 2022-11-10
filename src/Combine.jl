function combine_particles(a::Jet, b::Jet)::Jet
    new_momentum    =   a.Momentum + b.Momentum
    new_energy      =   a.Energy + b.Energy
    new_mass_sqr    =   new_energy^2 - sum(new_momentum.^2)
    new_mass        =   if new_mass_sqr ≥   0
        sqrt(new_mass_sqr)
    else
        printstyled(
            "Warning: The mass sqruared $new_mass_sqr is less than 0! We have set it to 0. Please check it!\n",
            color=:red
        )
        0
    end

    return  Jet(
        nothing,
        new_momentum,
        new_energy,
        new_mass,
        vcat(a.Particles, b.Particles)
    )
end
function combine_particles(a::AbstractParticle, b::AbstractParticle)::Jet
    return  combine_particles(
        Jet(a), Jet(b)
    )
end