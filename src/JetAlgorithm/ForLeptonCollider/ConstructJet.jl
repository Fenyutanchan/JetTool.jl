function construct_jets_for_lepton_collision(
    event::Union{Event, Jet},
    jet_radius::Real;
    alg::Symbol=:Spherical_Generalized_kT,
    kT_parameter::Real=-1
)::Union{Event, Jet}
    @assert alg ∈ [
        :Spherical_Generalized_kT,
        # :Spherical_kT,
        # :Spherical_Cambridge,
        # :JADE,
        # :Spherical_SISCone,
    ]

    combination_function  =   (eval ∘ Symbol)(
        "combination_" * String(alg)
    )
    
    parameters_list =   deleteat!(
        [jet_radius, kT_parameter],
        .![
            combination_function(jet_radius_flag),
            combination_function(kT_parameter_flag)
        ]
    )
    return combination_function(event, parameters_list...)
end

function combination_Spherical_Generalized_kT(
    event::Union{Event, Jet},
    jet_radius::Real,
    kT_parameter::Real=-1
)::Union{Event, Jet}
    in_event_particles  =   Jet.(event.Particles)
    out_event_particles =   Jet[]

    while !isempty(in_event_particles)
        min_distance    =   typemax(Float64)
        delete_flag     =   [0, 0, 0]
        for ii ∈ eachindex(in_event_particles)
            p_i             =   in_event_particles[ii]
            beam_distance_i =   distance_Spherical_Generalized_kT(
                p_i;
                p=kT_parameter
            )
            if beam_distance_i ≤ min_distance
                delete_flag     =   [ii, 0, 0]
                min_distance    =   beam_distance_i
            end
            for jj ∈ 1:(ii-1)
                p_j             =   in_event_particles[jj]
                distance_ii_jj  =   distance_Spherical_Generalized_kT(
                    p_i, p_j, jet_radius;
                    p=kT_parameter
                )
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

    return if typeof(event) == Event
        Event(
            event.Event_Weight,
            out_event_particles
        )
    elseif typeof(event) == Jet
        Jet(
            event.PDGID,
            event.Momentum,
            event.Energy,
            event.Mass,
            out_event_particles
        )
    end
end
combination_Spherical_Generalized_kT(::typeof(jet_radius_flag))::Bool   =   true
combination_Spherical_Generalized_kT(::typeof(kT_parameter_flag))::Bool =   true