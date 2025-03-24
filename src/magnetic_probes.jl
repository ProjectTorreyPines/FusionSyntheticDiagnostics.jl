"""
    compute_magnetic_probes!(
        PSI_interpolant::Interpolations.AbstractInterpolation,
        probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}},
    ) where {T<:Real}

Synthetic diagnostics for the poloidal magnetic probes saved in the data structure.
Returns a vector with, for each probe, the measured component of Bpol aligned with the probe.
Also, updates the data structure with time and synthetic data. Data is saved at global time.
"""
function compute_magnetic_probes!(
    PSI_interpolant::Interpolations.AbstractInterpolation,
    probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}},
) where {T <: Real}

    if isempty(probes)
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
        @ddtime(probes[k].field.data = b)
    end

    return B
end
