function combine_partons(event::Event, jet_radius::Real, p::Real)::Event
    in_event_particles  =   Jet.(event.Particles)
    out_event_particles =   Jet[]

    while !isempty(in_event_particles)
        min_distance    =   typemax(Float64)
        delete_flag     =   [0, 0, 0]
        for ii ∈ eachindex(in_event_particles)
            p_i             =   in_event_particles[ii]
            beam_distance_i =   transverse_momentum(p_i)^(2 * p) * jet_radius^2
            if beam_distance_i ≤ min_distance
                delete_flag =   [ii, 0, 0]
                min_distance    =   beam_distance_i
            end
            for jj ∈ 1:(ii-1)
                p_j             =   in_event_particles[jj]
                distance_ii_jj  =   calc_distance(p_i, p_j, p)
                if distance_ii_jj ≤ min_distance
                    delete_flag     =   [0, jj, ii]
                    min_distance    =   distance_ii_jj
                end
            end
        end
        @assert sum(delete_flag) != 0

        if (iszero ∘ first)(delete_flag)
            @assert !iszero(delete_flag[2]) && !iszero(delete_flag[3])
            p_i =   in_event_particles[delete_flag[2]]
            p_j =   in_event_particles[delete_flag[3]]
            deleteat!(in_event_particles, delete_flag[2:end])
            push!(
                in_event_particles,
                combine_partons(p_i, p_j)
            )
        else
            @assert iszero(delete_flag[2:end])
            push!(
                out_event_particles,
                Jet(
                    in_event_particles[
                        first(delete_flag)
                    ]
                )
            )
            deleteat!(
                in_event_particles,
                first(delete_flag)
            )
        end
    end

    return  Event(
        event.Event_Weight,
        out_event_particles
    )
end

function combine_partons(a::Jet, b::Jet)::Jet
    # println("Making parton combination.")
    
    new_PDGID   =   if b.PDGID == 21
        a.PDGID
    elseif a.PDGID == 21
        b.PDGID
    elseif (a.PDGID + b.PDGID) == 0
        21
    else
        println("Parton combination exception!")
        1
    end

    new_momentum    =   a.Momentum + b.Momentum
    new_energy      =   a.Energy + b.Energy
    new_mass        =   sqrt(new_energy^2 - norm(new_momentum)^2)

    # if typeof(a) == Jet
    #     if typeof(b) == Jet
    #         return  Jet(
    #             new_PDGID,
    #             new_momentum,
    #             new_energy,
    #             new_mass,
    #             vcat(a.Particles, b.Particles)
    #         )
    #     elseif typeof(b) == Particle
    #         return  Jet(
    #             new_PDGID,
    #             new_momentum,
    #             new_energy,
    #             new_mass,
    #             vcat(a.Particles, [b])
    #         )
    #     else
    #         throw("Exception for the type of AbstractParticle $(typeof(b)).")
    #     end
    # elseif typeof(a) == Particle
    #     if typeof(b) == Jet
    #         return  Jet(
    #             new_PDGID,
    #             new_momentum,
    #             new_energy,
    #             new_mass,
    #             vcat([a], b.Particles)
    #         )
    #     elseif typeof(b) == Particle
    #         return  Jet(
    #             new_PDGID,
    #             new_momentum,
    #             new_energy,
    #             new_mass,
    #             [a, b]
    #         )
    #     else
    #         throw("Exception for the type of AbstractParticle $(typeof(b)).")
    #     end
    # else
    #     throw("Exception for the type of AbstractParticle $(typeof(a)).")
    # end

    return  Jet(
        new_PDGID,
        new_momentum,
        new_energy,
        new_mass,
        vcat(a.Particles, b.Particles)
    )
end
function combine_partons(a::AbstractParticle, b::AbstractParticle)::Jet
    return  combine_partons(
        Jet(a), Jet(b)
    )
end

function combine_particles_for_events(
    event::Event, jet_radius::Real;
    alg::Symbol=:kT, p::Real=-1 # default for anti-kT
)::Event
    if alg == :kT
        return  combine_particles_for_events_kT(event, jet_radius, p)
    elseif alg == :cone
        return  combine_particles_for_events_cone(event, jet_radius)
    else
        throw("Jet algorithm not defined for $alg.")
    end
end

function combine_particles_for_events_kT(event::Event, jet_radius::Real, p::Real)::Event
    in_event_particles  =   Jet.(event.Particles)
    out_event_particles =   Jet[]

    while !isempty(in_event_particles)
        min_distance    =   typemax(Float64)
        delete_flag     =   [0, 0, 0]
        for ii ∈ eachindex(in_event_particles)
            p_i             =   in_event_particles[ii]
            beam_distance_i =   transverse_momentum(p_i)^(2 * p) * jet_radius^2
            if beam_distance_i ≤ min_distance
                delete_flag     =   [ii, 0, 0]
                min_distance    =   beam_distance_i
            end
            for jj ∈ 1:(ii-1)
                p_j             =   in_event_particles[jj]
                distance_ii_jj  =   calc_distance(p_i, p_j, p)
                if distance_ii_jj ≤ min_distance
                    delete_flag     =   [0, jj, ii]
                    min_distance    =   distance_ii_jj
                end
            end
        end
        @assert sum(delete_flag) != 0

        if (iszero ∘ first)(delete_flag)
            @assert !iszero(delete_flag[2]) && !iszero(delete_flag[3])
            p_i =   in_event_particles[delete_flag[2]]
            p_j =   in_event_particles[delete_flag[3]]
            deleteat!(in_event_particles, delete_flag[2:end])
            push!(
                in_event_particles,
                combine_particles(p_i, p_j)
            )
        else
            @assert iszero(delete_flag[2:end])
            push!(
                out_event_particles,
                Jet(
                    in_event_particles[
                        first(delete_flag)
                    ]
                )
            )
            deleteat!(
                in_event_particles,
                first(delete_flag)
            )
        end
    end

    return  Event(
        event.Event_Weight,
        out_event_particles
    )
end

function combine_particles_for_events_cone(event::Event, jet_radius::Real)::Event
    out_event_particles =   Jet.(event.Particles)
    dead_flags          =   [false for particle ∈ out_event_particles]

    while true
        combine_flag    =   false

        for ii ∈ eachindex(out_event_particles)
            if dead_flags[ii]
                continue
            end
            for jj ∈ 1:(ii-1)
                if dead_flags[jj]
                    continue
                end

                ΔR  =   calc_ΔR(
                    out_event_particles[ii],
                    out_event_particles[jj]
                )
                if ΔR ≤ jet_radius
                    out_event_particles[jj] =   combine_particles(
                        out_event_particles[jj],
                        out_event_particles[ii]
                    )

                    dead_flags[ii]  =   true
                    combine_flag    =   true
                    break
                end
            end
        end

        if !combine_flag
            break
        end
    end

    deleteat!(out_event_particles, dead_flags)

    return  Event(
        event.Event_Weight,
        out_event_particles
    )
end

function combine_particles(a::Jet, b::Jet)::Jet
    new_momentum    =   a.Momentum + b.Momentum
    new_energy      =   a.Energy + b.Energy
    new_mass        =   sqrt(new_energy^2 - norm(new_momentum)^2)

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