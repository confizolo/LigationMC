module SimulationMC

using Random
using JSON3
using Base.Threads
using ..PolymerUtils
using ..DSMC
using ..Network

export run_single_trial, run_trials_for_system

const FITTED_K1_DEFAULT = 1.0
const FITTED_K2_DEFAULT = 7928.458898884091
const FITTED_A_DEFAULT = 0.20927677484111143

function first_stage_to_half(fracs::Vector{Float64})::Union{Int, Nothing}
    for (idx, val) in enumerate(fracs)
        if val >= 0.5
            return idx
        end
    end
    return nothing
end

function run_stage!(network::NetworkBuilder, config::Dict, stage::Int, phi_ref::Float64, rng::AbstractRNG)
    # Keep the reference monomer density fixed while the ring population grows.
    total_ring_monomers = Float64(sum(network.ring_lengths))
    total_monomers = total_ring_monomers + Float64(config["mlin"] * config["nlin"])
    box_volume = total_monomers / phi_ref
    stage_L = box_volume^(1.0 / 3.0)
    
    linears = fill(config["nlin"], config["mlin"])
    stage_seed = abs(rand(rng, Int))
    events = run_dsmc!(linears, Float64(config["k1"]), Float64(config["k2"]);
                       alpha=Float64(config["alpha"]), nu=Float64(config["nu"]),
                       seed=stage_seed, max_steps=config["max_steps"])
    n_events = length(events)
    event_records = Dict{String, Any}[]
    
    for event in events
        if event isa CyclisationEvent
            event = process_cyclisation!(network, event, Float64(config["val_A"]), box_volume)
            push!(event_records, Dict{String, Any}(
                "stage" => stage,
                "L" => stage_L,
                "event_type" => "cyclisation",
                "time" => event.time,
                "linear_length" => event.linear_length,
                "ring_length" => event.ring_length,
                "links_formed" => event.links_formed,
                "linked_ring_ids" => event.linked_ring_ids,
                "ring_id" => event.ring_id
            ))
        else
            push!(event_records, Dict{String, Any}(
                "stage" => stage,
                "L" => stage_L,
                "event_type" => "merge",
                "time" => event.time,
                "length_i" => event.length_i,
                "length_j" => event.length_j,
                "new_length" => event.new_length
            ))
        end
    end
    
    frac = largest_component_fraction(network)
    return frac, stage_L, n_events, event_records
end

function run_single_trial(config::Dict{String, Any})::Dict{String, Any}
    rng = Random.Xoshiro(config["seed"])
    
    initial_ring_lengths = fill(config["nring"], config["mring"])
    network = NetworkBuilder(initial_ring_lengths, seed=config["seed"])
    network.rng = rng
    
    stage_fractions = Float64[]
    stage_events = Int[]
    stage_box_lengths = Float64[]
    event_timeline = Dict{String, Any}[]
    
    initial_monomers = Float64(config["mring"] * config["nring"] + config["mlin"] * config["nlin"])
    if initial_monomers <= 0.0
        error("Total monomers must be positive.")
    end
    phi_ref = initial_monomers / (Float64(config["L"])^3)
    
    # We use 1-based indexing for Julia (stage = 1,2,...,n_stages).
    for stage in 1:config["n_stages"]
        frac, stage_L, n_events, event_records = run_stage!(network, config, stage, phi_ref, rng)
        push!(stage_fractions, frac)
        push!(stage_box_lengths, stage_L)
        push!(stage_events, n_events)
        append!(event_timeline, event_records)
    end
    
    half_stage = first_stage_to_half(stage_fractions)
    
    return Dict{String, Any}(
        "seed" => config["seed"],
        "largest_component_fraction" => stage_fractions,
        "stages_to_half" => half_stage,
        "box_lengths" => stage_box_lengths,
        "events_per_stage" => stage_events
        # "degree_distribution" => degree_distribution(network) # uncomment if needed, but omitted to save space
    )
end

function run_trials_for_system(sys_cfg::Dict{String, Any}; trials::Int=100, n_stages::Int=100, k1::Float64=FITTED_K1_DEFAULT, k2::Float64=FITTED_K2_DEFAULT, alpha::Float64=1.0, nu::Float64=0.5, val_A::Float64=FITTED_A_DEFAULT, max_steps::Int=1000000)::Vector{Dict{String, Any}}
    
    config_base = Dict{String, Any}(
        "L" => sys_cfg["L"],
        "mring" => sys_cfg["mring"],
        "nring" => sys_cfg["nring"],
        "mlin" => sys_cfg["mlin"],
        "nlin" => sys_cfg["nlin"],
        "k1" => k1,
        "k2" => k2,
        "alpha" => alpha,
        "nu" => nu,
        "val_A" => val_A,
        "max_steps" => max_steps,
        "n_stages" => n_stages
    )
    
    results = Vector{Dict{String, Any}}(undef, trials)
    
    # Dynamic scheduling helps when some trials gel much earlier than others.
    @threads :dynamic for t in 1:trials
        config = copy(config_base)
        config["seed"] = 42 + t * 1337
        results[t] = run_single_trial(config)
    end
    
    return results
end

end
