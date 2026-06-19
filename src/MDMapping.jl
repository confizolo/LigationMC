module MDMapping

export map_k1_to_density

"""
    map_k1_to_density(k1::Float64, ntot::Int) -> Float64

Maps the Gillespie SSA merge rate constant `k1` to the Particle DSMC `density` parameter.
In Particle DSMC, the expected time increment for a specific pair (i,j) is:
    Δt = ntot / ( n_pairs * density * K_ij )
In Gillespie SSA, the rate for a specific unordered pair is:
    Rate = k1 * K_ij
Thus, the expected time to merge that pair is 1 / Rate.
Equating these gives:
    density = ntot * k1
"""
function map_k1_to_density(k1::Float64, ntot::Int)::Float64
    return k1 * Float64(ntot)
end

end
