module  JetTool

    using   LinearAlgebra
    using   StructParticle

    # export  combine_particles_for_events
    # export  energy_momentum_within_cone_r
    # export  transverse_momentum

    # include("Combine.jl")
    # include("Distance.jl")
    # include("EergyDistribution.jl")
    # include("ParticleParameter.jl")

    export construct_jets_for_lepton_collision

    (include ∘ joinpath)(
        "JetAlgorithm",
        "ForLeptonCollider",
        "ForLeptonCollider.jl"
    )

end # module JetTool
