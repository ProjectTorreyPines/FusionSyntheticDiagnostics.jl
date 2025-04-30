import Interpolations: AbstractInterpolation

export add_magnetics!, compute_magnetics!

const magnetic_probe_types = [
    :b_field_pol_probe,
    :b_field_tor_probe,
    :flux_loop,
    :rogowski_coil,
    :shunt,
]

const all_mag_probes = Union{
    IMAS.magnetics__b_field_pol_probe{T},
    IMAS.magnetics__b_field_tor_probe{T},
    IMAS.magnetics__flux_loop{T},
    IMAS.magnetics__rogowski_coil{T},
    IMAS.magnetics__shunt{T},
} where {T}

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
    compute_magnetics!(@nospecialize(ids::IMAS.dd)=IMAS.dd();
        equilibrium::Union{Nothing, IMAS.equilibrium}=nothing,
        kwargs...,
    )

Convinience function to compute all magnetic probes if required arguments are provided
in keywords.

For `b_field_pol_probe` and `flux_loops`,
`PSI_interpolant::AbstractInterpolation` is required which can be
obtained by calling
[`IMAS.ψ_interpolant`](https://projecttorreypines.github.io/IMAS.jl/dev/api/#IMAS.%CF%88_interpolant)

Alternatively, if an entire equilibrium is provided in keyword argument, then this
function will compute `PSI_interpolant` for each time slice and compute the probes
for each time step.
"""
function compute_magnetics!(@nospecialize(ids::IMAS.dd)=IMAS.dd();
    equilibrium::Union{Nothing, IMAS.equilibrium}=nothing,
    kwargs...,
)
    if !isnothing(equilibrium)
        for ts ∈ equilibrium.time_slice
            _, _, ψ_interp = IMAS.ψ_interpolant(ts.profiles_2d)
            compute_magnetics!(
                ids.magnetics.b_field_pol_probe;
                PSI_interpolant=ψ_interp,
                time=ts.time,
            )
            compute_magnetics!(
                ids.magnetics.flux_loop;
                PSI_interpolant=ψ_interp,
                time=ts.time,
            )
        end
    else
        for probe ∈ magnetic_probe_types
            try
                compute_magnetics!(getproperty(ids.magnetics, probe); kwargs...)
            catch
                continue
            end
        end
    end
end

"""
    compute_magnetics!(
        probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}};
        PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
        time::Union{Nothing, T}=nothing,
    ) where {T <: Float64}

Synthetic diagnostics for the poloidal magnetic probes saved in the data structure.
Returns a vector with, for each probe, the measured component of Bpol aligned with the
probe. Also, updates the data structure with time and synthetic data. Data is saved at
global time if `time` is not provided.
"""
function compute_magnetics!(
    probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}};
    PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
    time::Union{Nothing, T}=nothing,
) where {T <: Float64}
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
    _write_with_time(probes, :field, time, B)

    return B
end

"""
    compute_magnetics!(
        loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}};
        PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
        time::Union{Nothing, T}=nothing,
    ) where {T <: Float64}

Synthetic diagnostics for the flux loops saved in the data structure.
Returns a vector with, for each loop, the measured poloidal flux.
Also, updates the data structure with time and synthetic data.
Data is saved at global time if `time` is not provided.
"""
function compute_magnetics!(
    loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}};
    PSI_interpolant::Union{Nothing, AbstractInterpolation}=nothing,
    time::Union{Nothing, T}=nothing,
) where {T <: Float64}
    #
    if isempty(loops) || isnothing(PSI_interpolant)
        return T[]
    end

    # get location of the loops
    r = [mean([pos.r for pos ∈ loop.position]) for loop ∈ loops]
    z = [mean([pos.z for pos ∈ loop.position]) for loop ∈ loops]

    # psi at the loop positions
    psi = PSI_interpolant.(r, z)

    # update dd with synthetic data
    _write_with_time(loops, :flux, time, psi[1] .- psi)

    return psi[1] .- psi
end

"""
    _write_with_time(
        probes::IMAS.IDSvector{<:all_mag_probes},
        value_filed::Symbol,
        time::Union{Nothing, T},
        data::Vector{T},
    ) where {T <: Float64}

Internal function to write probe value data. If time is not provided, global time is
used but if time is provided, data is written along with time.
"""
function _write_with_time(
    probes::IMAS.IDSvector{<:all_mag_probes},
    value_filed::Symbol,
    time::Union{Nothing, T},
    data::Vector{T},
) where {T <: Float64}
    if isnothing(time)
        for (probe, value) ∈ zip(probes, data)
            pv = getfield(probe, value_filed)
            IMAS.@ddtime(pv.data = value)
        end
    else
        for (probe, value) ∈ zip(probes, data)
            pv = getfield(probe, value_filed)
            if ismissing(pv, :time)
                pv.time = T[]
            end
            if ismissing(pv, :data)
                pv.data = T[]
            end
            append!(pv.time, time)
            append!(pv.data, value)
        end
    end
end

"""
To do:

Add following methods for compute_magnetics! that are specific to other probes:

compute_magnetics!(probes::IMAS.IDSvector{IMAS.magnetics__b_field_tor_probe{T}})
compute_magnetics!(coils::IMAS.IDSvector{IMAS.magnetics__rogowski_coil{T}})
compute_magnetics!(shunts::IMAS.IDSvector{IMAS.magnetics__shunt{T}})
"""
