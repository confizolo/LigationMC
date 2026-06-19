module ParticleDSMC

using Random
using ..PolymerUtils
using ..GillespieSSA: Event, MergeEvent, CyclisationEvent

export run_particle_dsmc!

function run_particle_dsmc!(
    linear_lengths::Vector{Int},
    density::Float64,
    k2::Float64;
    alpha::Float64=1.0,
    nu::Float64=0.5,
    seed::Union{Int, Nothing}=nothing,
    max_steps::Int=500000
)
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    
    ntot = length(linear_lengths)
    masses = copy(linear_lengths)
    
    # We maintain a list of active indices for faster sampling of `i` and `k`
    # However, standard DSMC picks `j` uniformly over all slots.
    
    n_chains = ntot
    time = 0.0
    
    # Initial max rates estimates (must be > 0)
    # We can use the rates of the initial species
    k_max = smoluchowski_kernel(linear_lengths[1], linear_lengths[1]; alpha=alpha, nu=nu)
    r_max = cyclisation_rate(linear_lengths[1], k2; nu=nu)
    if k_max == 0.0; k_max = 1e-4; end
    if r_max == 0.0; r_max = 1e-4; end
    
    Alpha = 1.0 # time splitting parameter
    
    events = Event[]
    steps = 0
    
    while n_chains > 0
        if steps >= max_steps
            error("Exceeded max_steps=$max_steps before exhausting linears.")
        end
        steps += 1
        
        # In the original python script, it stops when n_chains > 1, but we need to cyclise the last chain too!
        # If n_chains == 1, p_ann should be 0 because (n_chains-1) is 0.
        p_ann = 0.0
        if n_chains > 1
            p_ann = 1.0 / (1.0 + (2.0 * ntot * r_max) / ((n_chains - 1) * density * k_max))
        end
        
        if rand(rng) < p_ann
            # Attempt annealing
            active_indices = findall(x -> x > 0, masses)
            i = rand(rng, active_indices)
            j = rand(rng, 1:ntot)
            
            # In the Python script, j is picked uniformly from all slots, then rejected if 0 or == i
            while masses[j] == 0 || j == i
                j = rand(rng, 1:ntot)
            end
            
            mi = masses[i]
            mj = masses[j]
            k_ij = smoluchowski_kernel(mi, mj; alpha=alpha, nu=nu)
            
            if k_ij > k_max
                k_max = k_ij
            else
                if rand(rng) < k_ij / k_max
                    time += 2.0 * Alpha * ntot / (n_chains * (n_chains - 1) * density * k_ij)
                    masses[j] = mi + mj
                    masses[i] = 0
                    n_chains -= 1
                    push!(events, MergeEvent(time, mi, mj, mi + mj))
                end
            end
        else
            # Attempt cyclisation
            active_indices = findall(x -> x > 0, masses)
            k = rand(rng, active_indices)
            mk = masses[k]
            rmk = cyclisation_rate(mk, k2; nu=nu)
            
            if rmk > r_max
                r_max = rmk
            else
                if rand(rng) < rmk / r_max
                    time += (1.0 - Alpha) / (n_chains * rmk)
                    masses[k] = 0
                    n_chains -= 1
                    push!(events, CyclisationEvent(time, mk, mk))
                end
            end
        end
    end
    
    return events
end

end
