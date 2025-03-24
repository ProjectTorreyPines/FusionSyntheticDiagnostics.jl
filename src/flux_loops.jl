"""
    compute_flux_loops!(
        PSI_interpolant::Interpolations.AbstractInterpolation,
        loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}},
    ) where {T <: Real}

Synthetic diagnostics for the flux loops saved in the data structure.
Returns a vector with, for each loop, the measured poloidal flux.
Also, updates the data structure with time and synthetic data. Data is saved at global time.
"""
function compute_flux_loops!(
    PSI_interpolant::Interpolations.AbstractInterpolation,
    loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}},
) where {T <: Real}
    #
    if isempty(loops)
        return T[]
    end

    # get location of the loops
    r = [loop.position.r for loop ∈ loops]
    z = [loop.position.z for loop ∈ loops]

    # psi at the loop positions
    psi = PSI_interpolant.(r, z)

    # update dd with synthetic data
    for (k, p) ∈ enumerate(psi)
        @ddtime(loops[k].flux.data = p[1] - p)
    end

    return psi[1] .- psi
end
