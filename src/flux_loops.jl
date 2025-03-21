"""
    compute_flux_loops(PSI_interpolant::Interpolations.AbstractInterpolation,
                                 loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}}) where T

Synthetic diagnostics for the flux loops saved in the data structure.
Returns a vector with, for each loop, the measured poloidal flux.
Also, updates the data structure with time and synthetic data. Data is saved at global time.
"""

function compute_flux_loops(PSI_interpolant::Interpolations.AbstractInterpolation,
                                 loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}}) where T
    @assert !isempty(loops) "Magnetic loops not inside the data structure. Impossible to create synthetic diagnostics"    

    # get location of the loops 
    psi = [PSI_interpolant(loop.position[1].r,loop.position[1].z) for loop in loops]

    # update dd with synthetic data
    for (k,p) in enumerate(psi)
        @ddtime(loops[k].flux.data = p)
    end
    return psi[1].-psi

end


"""
    location_flux_loops(loop::IMAS.magnetics__flux_loop{T}) where T

Returns the location of a flux loop saved in the data structure.
"""

function location_flux_loops(loop::IMAS.magnetics__flux_loop{T}) where T
    @assert !isempty(loop) "Flux loop not inside the data structure. Impossible to retrieve position of loop"    
    return (loop.position[1].r,loop.position[1].z)
end

function location_flux_loops(loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}}) where T
    @assert !isempty(loops) "Flux loops not inside the data structure. Impossible to retrieve position of loops"    
    return location_flux_loops.(loops)
end


"""
    name_flux_loops(loop::IMAS.magnetics__flux_loop{T}) where T

Returns the name of a flux loop saved in the data structure.
"""
function name_flux_loops(loop::IMAS.magnetics__flux_loop{T}) where T
    @assert !isempty(loop) "Flux loops not inside the data structure. Impossible to retrieve name of loops"  
    return loop.name
end

function name_flux_loops(loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}}) where T
    @assert !isempty(loops) "Flux loops not inside the data structure. Impossible to retrieve name of loops"  
    return name_flux_loops.(loops)
end


"""
    data_flux_loops(loop::IMAS.magnetics__flux_loop{T}) where T

Returns the (time, data) of a flux loop saved in the data structure.
"""
function data_flux_loops(loop::IMAS.magnetics__flux_loop{T}) where T
    @assert !isempty(loop) "Flux loop not inside the data structure. Impossible to get data of loop"  
    @assert !isempty(loop.flux.data) "No data is saved for the flux loop in the data structure"  
    return (time = loop.flux.time, data = loop.flux.data)
end

function data_flux_loops(loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}}) where T
    @assert !isempty(loops) "Flux loops not inside the data structure. Impossible to get data of loop"  
    return data_flux_loops.(loops)
end
