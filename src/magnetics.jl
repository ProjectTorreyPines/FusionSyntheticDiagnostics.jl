import Interpolations: AbstractInterpolation

export add_magnetics!, compute_magnetics!

const magnetic_probe_types = [
    :b_field_pol_probe,
    :b_field_tor_probe,
    :flux_loop,
    :rogowski_coil,
    :shunt,
]

"""
    add_magnetics!(
        config::Union{String, Dict{Symbol, Any}},
        @nospecialize(ids::IMAS.dd)=IMAS.dd();
        overwrite::Bool=false,
        kwargs...,
    )::IMAS.dd

Add magnetic diagnostics to IMAS structure using a `JSON` file or Julia `Dict` and
compute their outputs `kwargs` are passed to [`compute_magnetics!`](@ref).
"""
function add_magnetics!(
    config::Union{String, Dict{Symbol, Any}},
    @nospecialize(ids::IMAS.dd)=IMAS.dd();
    overwrite::Bool=false,
    kwargs...,
)::IMAS.dd
    add_diagnostic!(
        config,
        :magnetics,
        ids;
        overwrite=overwrite,
        channel=magnetic_probe_types,
    )
    compute_magnetics!(ids; kwargs...)
    return ids
end

"""
    compute_magnetics!(@nospecialize(ids::IMAS.dd)=IMAS.dd(); kwargs...)

Convinience function to compute all magnetic probes if required arguments are provided
in keywords.

For `b_field_pol_probe` and `flux_loops`,
`PSI_interpolant::AbstractInterpolation` is required which can be
obtained by calling
[IMAS.Ψ_interpolant](https://projecttorreypines.github.io/IMAS.jl/dev/api/#IMAS.%CF%88_interpolant)
"""
function compute_magnetics!(@nospecialize(ids::IMAS.dd)=IMAS.dd(); kwargs...)
    for probe ∈ magnetic_probe_types
        try
            compute_magnetics!(getproperty(ids.magnetics, probe); kwargs...)
        catch
            continue
        end
    end
end

"""
    compute_magnetics!(
        probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}};
        PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
    ) where {T <: Real}

Synthetic diagnostics for the poloidal magnetic probes saved in the data structure.
Returns a vector with, for each probe, the measured component of Bpol aligned with the
probe. Also, updates the data structure with time and synthetic data. Data is saved at
global time.
"""
function compute_magnetics!(
    probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}};
    PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
) where {T <: Real}
    #
    if isempty(probes) || isnothing(PSI_interpolant)
        return T[]
    end

    # get location of the probes and their poloidal angle
    r = [probe.position.r for probe ∈ probes]
    z = [probe.position.z for probe ∈ probes]
    poloidal_angle = [probe.poloidal_angle for probe ∈ probes] # clock-wise angle from the horizontal

    # local Br and Bz at probe location
    Br, Bz = IMAS.Br_Bz(PSI_interpolant, r, z)

    #component of polodial field measured by the probe 
    B = Br .* cos.(poloidal_angle) .- Bz .* sin.(poloidal_angle)

    # update dd with synthetic data
    for (k, b) ∈ enumerate(B)
        IMAS.@ddtime(probes[k].field.data = b)
    end

    return B
end

"""
    compute_magnetics!(
        loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}};
        PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
    ) where {T <: Real}

Synthetic diagnostics for the flux loops saved in the data structure.
Returns a vector with, for each loop, the measured poloidal flux.
Also, updates the data structure with time and synthetic data. Data is saved at global time.
"""
function compute_magnetics!(
    loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}};
    PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
) where {T <: Real}
    #
    if isempty(loops) || isnothing(PSI_interpolant)
        return T[]
    end

    # get location of the loops
    r = [loop.position.r for loop ∈ loops]
    z = [loop.position.z for loop ∈ loops]

    # psi at the loop positions
    psi = PSI_interpolant.(r, z)

    # update dd with synthetic data
    for (k, p) ∈ enumerate(psi)
        IMAS.@ddtime(loops[k].flux.data = p[1] - p)
    end

    return psi[1] .- psi
end

"""
To do:

Add following methods for compute_magnetics! that are specific to other probes:

compute_magnetics!(probes::IMAS.IDSvector{IMAS.magnetics__b_field_tor_probe{T}})
compute_magnetics!(coils::IMAS.IDSvector{IMAS.magnetics__rogowski_coil{T}})
compute_magnetics!(shunts::IMAS.IDSvector{IMAS.magnetics__shunt{T}})
"""
