module DSMC

using Random
using ..PolymerUtils

export Event, MergeEvent, CyclisationEvent, DSMCEngine, step!, run_until_exhausted!

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

mutable struct DSMCEngine
    k1::Float64
    k2::Float64
    alpha::Float64
    nu::Float64
    time::Float64
    
    linear_population::Dict{Int, Int}
    ring_population::Vector{Int}
    rng::AbstractRNG
    n_linear::Int
end

function DSMCEngine(linear_lengths::Vector{Int}, k1::Float64, k2::Float64; alpha::Float64=1.0, nu::Float64=0.5, seed::Union{Int, Nothing}=nothing)
    linear_pop = Dict{Int, Int}()
    for len in linear_lengths
        linear_pop[len] = get(linear_pop, len, 0) + 1
    end
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    n_linear = length(linear_lengths)
    return DSMCEngine(k1, k2, alpha, nu, 0.0, linear_pop, Int[], rng, n_linear)
end

function reaction_channels(engine::DSMCEngine)
    # Build all currently allowed reactions and their propensities.
    # We enumerate cyclisation channels (single length) and merge channels (length pairs).
    lengths = sort(collect(keys(engine.linear_population)))
    channels = Tuple{Symbol, Int, Int}[]
    propensities = Float64[]
    
    for (idx, i) in enumerate(lengths)
        ni = engine.linear_population[i]
        if ni <= 0
            continue
        end
        
        a_cyc = ni * cyclisation_rate(i, engine.k2; nu=engine.nu)
        if a_cyc > 0.0
            push!(channels, (:cyc, i, i))
            push!(propensities, a_cyc)
        end
        
        for j in lengths[idx:end]
            nj = engine.linear_population[j]
            if nj <= 0
                continue
            end
            if i == j && ni < 2
                continue
            end
            
            factor = ni * (nj - (i == j ? 1 : 0))
            if factor <= 0
                continue
            end
            
            a_merge = engine.k1 * factor * smoluchowski_kernel(i, j; alpha=engine.alpha, nu=engine.nu)
            if a_merge > 0.0
                push!(channels, (:merge, i, j))
                push!(propensities, a_merge)
            end
        end
    end
    return channels, propensities
end

function step!(engine::DSMCEngine)::Event
    if engine.n_linear <= 0
        error("No linear polymers left to evolve.")
    end
    
    channels, prop = reaction_channels(engine)
    if isempty(prop)
        error("No available reaction channels while linears remain.")
    end
    
    a0 = sum(prop)
    # Standard Gillespie waiting time draw.
    tau = randexp(engine.rng) / a0
    engine.time += tau
    
    choice = rand(engine.rng) * a0
    cumulative = 0.0
    picked = 0
    for (idx, a) in enumerate(prop)
        cumulative += a
        if choice <= cumulative
            picked = idx
            break
        end
    end
    # Fallback to last index due to precision issues
    if picked == 0
        picked = length(prop)
    end
    
    channel, i, j = channels[picked]
    
    if channel == :merge
        # Consumes two linears and creates one longer linear.
        engine.linear_population[i] -= 1
        if engine.linear_population[i] == 0
            delete!(engine.linear_population, i)
        end
        
        engine.linear_population[j] = get(engine.linear_population, j, 0) - 1
        if engine.linear_population[j] == 0
            delete!(engine.linear_population, j)
        end
        
        new_length = i + j
        engine.linear_population[new_length] = get(engine.linear_population, new_length, 0) + 1
        engine.n_linear -= 1
        
        return MergeEvent(engine.time, i, j, new_length)
    else
        # Cyclisation consumes one linear and creates one new ring entry.
        engine.linear_population[i] -= 1
        if engine.linear_population[i] == 0
            delete!(engine.linear_population, i)
        end
        engine.n_linear -= 1
        
        push!(engine.ring_population, i)
        return CyclisationEvent(engine.time, i, i)
    end
end

function run_until_exhausted!(engine::DSMCEngine; max_steps::Int=50000)
    # Run SSA events until there are no linear chains left in the stage.
    events = Event[]
    steps = 0
    while engine.n_linear > 0
        if steps >= max_steps
            error("Exceeded max_steps=$max_steps before exhausting linears.")
        end
        push!(events, step!(engine))
        steps += 1
    end
    return events
end

end
