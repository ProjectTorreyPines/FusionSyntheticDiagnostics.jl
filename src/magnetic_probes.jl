"""
    compute_magnetic_probes(PSI_interpolant::Interpolations.AbstractInterpolation,
                                 probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}}) where T

Synthetic diagnostics for the poloidal magnetic probes saved in the data structure.
Returns a vector with, for each probe, the measured component of Bpol aligned with the probe.
Also, updates the data structure with time and synthetic data. Data is saved at global time.
"""

function compute_magnetic_probes(PSI_interpolant::Interpolations.AbstractInterpolation,
                                 probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}}) where T
    @assert !isempty(probes) "Magnetic probes not inside the data structure. Impossible to create synthetic diagnostics"    

    # get location of the probes and their poloidal angle
    r=[probe.position.r for probe in probes] 
    z=[probe.position.z for probe in probes]
    poloidal_angle = [probe.poloidal_angle for probe in probes] # clock-wise angle from the horizontal

    Br, Bz = IMAS.Br_Bz(PSI_interpolant,r,z) # local Br and Bz at probe location

    #component of polodial field measured by the probe 
    # B = n ̇Bp = Br cos(angle) - Bz sin(angle)
    B = Br.*cos.(poloidal_angle) .- Bz.*sin.(poloidal_angle)
    # update dd with synthetic data
    for (k,b) in enumerate(B)
        @ddtime(probes[k].field.data = b)
    end
    return B

end


"""
    location_magnetic_probes(probe::IMAS.magnetics__b_field_pol_probe{T}) where T

Returns the location of a poloidal magnetic probe saved in the data structure.
"""

function location_magnetic_probes(probe::IMAS.magnetics__b_field_pol_probe{T}) where T
    @assert !isempty(probe) "Magnetic probe not inside the data structure. Impossible to retrieve position of probe"    
    return (probe.position.r,probe.position.z)
end

function location_magnetic_probes(probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}}) where T
    @assert !isempty(probes) "Magnetic probes not inside the data structure. Impossible to retrieve position of probes"    
    return location_magnetic_probes.(probes)
end


"""
    name_magnetic_probes(probes::IMAS.magnetics__b_field_pol_probe{T}) where T

Returns the name of a poloidal magnetic probe saved in the data structure.
"""
function name_magnetic_probes(probe::IMAS.magnetics__b_field_pol_probe{T}) where T
    @assert !isempty(probe) "Magnetic probes not inside the data structure. Impossible to retrieve name of probes"  
    return probe.name 
end

function name_magnetic_probes(probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}}) where T
    @assert !isempty(probes) "Magnetic probes not inside the data structure. Impossible to retrieve name of probes"  
    return name_magnetic_probes.(probes)
end

"""
    data_magnetic_probes(probe::IMAS.magnetics__b_field_pol_probe{T}) where T

Returns the (time, data) of a poloidal magnetic probe saved in the data structure.
"""
function data_magnetic_probes(probe::IMAS.magnetics__b_field_pol_probe{T}) where T
    @assert !isempty(probe) "Magnetic probe not inside the data structure. Impossible to get data of probe"  
    @assert !isempty(probe.field.data) "No data is saved for the magentic probe in the data structure"  
    return (time = probe.field.time, data = probe.field.data)
end

function data_magnetic_probes(probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}}) where T
    @assert !isempty(probes) "Magnetic probe not inside the data structure. Impossible to get data of probe"  
    return data_magnetic_probes.(probes)
end
