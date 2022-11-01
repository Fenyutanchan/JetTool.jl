function transverse_momentum_within_cone_r(j::Jet, r::Real)::Real
    result  =   0
    for particle ∈ j.Particles
        ΔR  =   calc_ΔR(j, particle)

        if ΔR ≤ r
            result  +=  transverse_momentum(particle)
        end
    end
    return  result
end

function energy_momentum_within_cone_r(j::Jet, r::Real)::Real
    result  =   0
    for particle ∈ j.Particles
        ΔR  =   calc_ΔR(j, particle)

        if ΔR ≤ r
            result  +=  particle.Energy
        end
    end
    return  result
end