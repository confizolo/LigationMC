module Network

using Graphs
using Random
using ..GillespieSSA
using ..PolymerUtils

export NetworkBuilder, add_ring!, process_cyclisation!, largest_component_fraction, degree_distribution

mutable struct NetworkBuilder
    graph::SimpleGraph{Int}
    ring_lengths::Vector{Int}
    rng::AbstractRNG
end

function NetworkBuilder(initial_ring_lengths::Vector{Int}; seed::Union{Int, Nothing}=nothing)
    graph = SimpleGraph{Int}(length(initial_ring_lengths))
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    # The nodes are pre-added by passing N to SimpleGraph
    return NetworkBuilder(graph, copy(initial_ring_lengths), rng)
end

function add_ring!(builder::NetworkBuilder, length::Int)::Int
    add_vertex!(builder.graph)
    push!(builder.ring_lengths, length)
    return nv(builder.graph) # The ID of the newly added vertex
end

function process_cyclisation!(builder::NetworkBuilder, event::CyclisationEvent, A::Float64, box_volume::Float64=1.0)
    # Add the freshly cyclised ring as a new graph node.
    new_ring_id = add_ring!(builder, event.ring_length)
    l_cyc = event.ring_length
    
    n_total = nv(builder.graph) - 1 # excluding the new ring
    if n_total == 0
        event.ring_id = new_ring_id
        return event
    end
    
    targets = Int[]
    for target in 1:n_total
        nring_target = builder.ring_lengths[target]
        if nring_target <= 0
            continue
        end
        
        # Bernoulli-per-target approximation from Poisson mean mu.
        mu = A * Float64(nring_target) * Float64(l_cyc) / box_volume
        if mu > 0.0
            p = 1.0 - exp(-mu)
            if rand(builder.rng) < p
                add_edge!(builder.graph, new_ring_id, target)
                push!(targets, target)
            end
        end
    end
    
    event.links_formed = length(targets)
    event.linked_ring_ids = targets
    event.ring_id = new_ring_id
    
    return event
end

function largest_component_fraction(builder::NetworkBuilder)::Float64
    n = nv(builder.graph)
    if n == 0
        return 0.0
    end
    cc = connected_components(builder.graph)
    if isempty(cc)
        return 0.0
    end
    max_len = maximum(length, cc)
    return Float64(max_len) / Float64(n)
end

function degree_distribution(builder::NetworkBuilder)::Dict{Int, Int}
    deg_dist = Dict{Int, Int}()
    for v in vertices(builder.graph)
        d = degree(builder.graph, v)
        deg_dist[d] = get(deg_dist, d, 0) + 1
    end
    return deg_dist
end

end
