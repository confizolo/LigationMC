# DSMC — Particle DSMC Engine for Polymer Growth
#
# Acceptance-rejection particle DSMC for the irreversible Smoluchowski
# equation with cyclisation.  Each linear chain is tracked individually
# (not as a population bucket).  At each step the algorithm picks between
# an annealing (merge) attempt and a cyclisation attempt with probability
# proportional to the current majorant rates, then accepts or rejects
# based on the true rate / majorant ratio.

module DSMC

using Random
using ..PolymerUtils

export Event, MergeEvent, CyclisationEvent, run_dsmc!

# ── Shared event types ────────────────────────────────────────────────

abstract type Event end

struct MergeEvent <: Event
    time::Float64
    length_i::Int
    length_j::Int
    new_length::Int
end

mutable struct CyclisationEvent <: Event
    time::Float64
    linear_length::Int
    ring_length::Int
    links_formed::Int
    linked_ring_ids::Vector{Int}
    ring_id::Union{Int, Nothing}
end

function CyclisationEvent(time::Float64, linear_length::Int, ring_length::Int)
    return CyclisationEvent(time, linear_length, ring_length, 0, Int[], nothing)
end

# ── Engine ────────────────────────────────────────────────────────────

"""
    run_dsmc!(linear_lengths, k1, k2; kwargs...) -> Vector{Event}

Run the particle DSMC until all linear chains have either merged or
cyclised.  The `density` parameter used internally is derived from
`k1` via `density = k1 * ntot` where `ntot = length(linear_lengths)`.

Returns a vector of `MergeEvent` / `CyclisationEvent` records.
"""
function run_dsmc!(
    linear_lengths::Vector{Int},
    k1::Float64,
    k2::Float64;
    alpha::Float64=1.0,
    nu::Float64=0.5,
    seed::Union{Int, Nothing}=nothing,
    max_steps::Int=500000,
    initial_time::Float64=0.0
)
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)

    ntot = length(linear_lengths)
    density = k1 * Float64(ntot)          # k1 → density mapping
    masses = copy(linear_lengths)

    n_chains = ntot
    time = initial_time

    # Initial majorant rate estimates (must be > 0)
    k_max = smoluchowski_kernel(linear_lengths[1], linear_lengths[1]; alpha=alpha, nu=nu)
    r_max = cyclisation_rate(linear_lengths[1], k2; nu=nu)
    if k_max == 0.0; k_max = 1e-4; end
    if r_max == 0.0; r_max = 1e-4; end

    Alpha = 1.0   # time-splitting parameter

    events = Event[]
    steps = 0

    while n_chains > 0
        if steps >= max_steps
            error("Exceeded max_steps=$max_steps before exhausting linears.")
        end
        steps += 1

        # Branching probability: merge vs cyclise
        p_ann = 0.0
        if n_chains > 1
            p_ann = 1.0 / (1.0 + (2.0 * ntot * r_max) / ((n_chains - 1) * density * k_max))
        end

        if rand(rng) < p_ann
            # ── Attempt merge ─────────────────────────────────────
            active_indices = findall(x -> x > 0, masses)
            i = rand(rng, active_indices)
            j = rand(rng, 1:ntot)
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
            # ── Attempt cyclisation ───────────────────────────────
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

    return events, time
end

end
