module  JetTool

    using   LinearAlgebra
    using   StructParticle

    # export  combine_particles_for_events
    # export  energy_momentum_within_cone_r
    # export  transverse_momentum

    # include("Distance.jl")
    # include("EergyDistribution.jl")
    # include("ParticleParameter.jl")

    export construct_jets_for_lepton_collision

    include("Combine.jl")

    (include ∘ joinpath)(
        "JetAlgorithm",
        "JetAlgorithm.jl"
    )

end # module JetTool
