include("PolymerUtils.jl")
include("GillespieSSA.jl")
include("ParticleDSMC.jl")
include("Network.jl")
include("MDMapping.jl")

module CompareDSMC

using ArgParse
using JSON3
using Statistics
using Random
using ..PolymerUtils
using ..GillespieSSA
using ..ParticleDSMC
using ..Network
using ..MDMapping

const FITTED_K1_DEFAULT = 1.0
const FITTED_K2_DEFAULT = 12840.849325710129
const FITTED_A_DEFAULT = 0.20927677484111143

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--L"
            arg_type = Float64
            default = 80.0
        "--mring"
            arg_type = Int
            default = 27
        "--nring"
            arg_type = Int
            default = 1024
        "--mlin"
            arg_type = Int
            default = 6
        "--nlin"
            arg_type = Int
            default = 64
        "--trials"
            arg_type = Int
            default = 100
        "--out_dir"
            arg_type = String
            default = "./results/compare"
        "--n_stages"
            arg_type = Int
            default = 100
    end
    return parse_args(s)
end

function run_gillespie_stage!(network::NetworkBuilder, config::Dict, stage::Int, phi_ref::Float64, rng::AbstractRNG)
    total_ring_monomers = Float64(sum(network.ring_lengths))
    total_monomers = total_ring_monomers + Float64(config["mlin"] * config["nlin"])
    box_volume = total_monomers / phi_ref
    stage_L = box_volume^(1.0 / 3.0)
    
    linears = fill(config["nlin"], config["mlin"])
    engine = DSMCEngine(linears, config["k1"], config["k2"]; alpha=config["alpha"], nu=config["nu"], seed=nothing)
    engine.rng = rng
    
    events = run_until_exhausted!(engine, max_steps=config["max_steps"])
    
    for event in events
        if event isa CyclisationEvent
            process_cyclisation!(network, event, config["val_A"], box_volume)
        end
    end
    
    return largest_component_fraction(network), events
end

function run_particle_stage!(network::NetworkBuilder, config::Dict, stage::Int, phi_ref::Float64, rng::AbstractRNG)
    total_ring_monomers = Float64(sum(network.ring_lengths))
    total_monomers = total_ring_monomers + Float64(config["mlin"] * config["nlin"])
    box_volume = total_monomers / phi_ref
    stage_L = box_volume^(1.0 / 3.0)
    
    linears = fill(config["nlin"], config["mlin"])
    density = map_k1_to_density(config["k1"], config["mlin"])
    
    # run_particle_dsmc! already accepts density. We need a way to pass the RNG.
    # But wait, run_particle_dsmc! accepts seed, not rng. Let's patch ParticleDSMC to accept rng!
    # For now, we pass seed and it creates its own rng, or we pass a seed derived from the trial seed + stage.
    stage_seed = abs(rand(rng, Int))
    events = run_particle_dsmc!(linears, density, config["k2"]; alpha=config["alpha"], nu=config["nu"], seed=stage_seed, max_steps=config["max_steps"])
    
    for event in events
        if event isa CyclisationEvent
            process_cyclisation!(network, event, config["val_A"], box_volume)
        end
    end
    
    return largest_component_fraction(network), events
end

function run_comparison_trial(config::Dict)
    # Gillespie
    rng_g = Random.Xoshiro(config["seed"])
    net_g = NetworkBuilder(fill(config["nring"], config["mring"]), seed=config["seed"])
    net_g.rng = rng_g
    
    # Particle DSMC
    rng_p = Random.Xoshiro(config["seed"])
    net_p = NetworkBuilder(fill(config["nring"], config["mring"]), seed=config["seed"])
    net_p.rng = rng_p
    
    phi_ref = Float64(config["mring"] * config["nring"] + config["mlin"] * config["nlin"]) / (config["L"]^3)
    
    fracs_g = Float64[]
    fracs_p = Float64[]
    
    events_g = []
    events_p = []
    
    for stage in 1:config["n_stages"]
        fg, eg = run_gillespie_stage!(net_g, config, stage, phi_ref, rng_g)
        push!(fracs_g, fg)
        append!(events_g, eg)
        
        fp, ep = run_particle_stage!(net_p, config, stage, phi_ref, rng_p)
        push!(fracs_p, fp)
        append!(events_p, ep)
    end
    
    half_g = findfirst(x -> x >= 0.5, fracs_g)
    half_p = findfirst(x -> x >= 0.5, fracs_p)
    
    # Collect lengths of cyclised rings
    rings_g = [e.ring_length for e in events_g if e isa CyclisationEvent]
    rings_p = [e.ring_length for e in events_p if e isa CyclisationEvent]
    
    # Valence
    val_g = [e.links_formed for e in events_g if e isa CyclisationEvent]
    val_p = [e.links_formed for e in events_p if e isa CyclisationEvent]
    
    return Dict(
        "gillespie" => Dict("half" => half_g, "rings" => rings_g, "valence" => val_g),
        "particle" => Dict("half" => half_p, "rings" => rings_p, "valence" => val_p)
    )
end

function main()
    args = parse_commandline()
    config = Dict{String, Any}(
        "L" => args["L"], "mring" => args["mring"], "nring" => args["nring"],
        "mlin" => args["mlin"], "nlin" => args["nlin"], "trials" => args["trials"],
        "n_stages" => args["n_stages"], "k1" => FITTED_K1_DEFAULT, "k2" => FITTED_K2_DEFAULT,
        "alpha" => 1.0, "nu" => 0.5, "val_A" => FITTED_A_DEFAULT, "max_steps" => 1000000
    )
    
    out_dir = args["out_dir"]
    mkpath(out_dir)
    
    results = []
    for t in 1:config["trials"]
        config["seed"] = 42 + t * 1337
        push!(results, run_comparison_trial(config))
    end
    
    # Summarize
    halfs_g = [r["gillespie"]["half"] for r in results if r["gillespie"]["half"] !== nothing]
    halfs_p = [r["particle"]["half"] for r in results if r["particle"]["half"] !== nothing]
    
    mg = isempty(halfs_g) ? nothing : mean(halfs_g)
    mp = isempty(halfs_p) ? nothing : mean(halfs_p)
    
    println("Gillespie mean stages-to-half: $mg ($(length(halfs_g)) gelled)")
    println("Particle  mean stages-to-half: $mp ($(length(halfs_p)) gelled)")
    
    # PMF of ring lengths
    all_rings_g = Int[]
    for r in results; append!(all_rings_g, r["gillespie"]["rings"]); end
    
    all_rings_p = Int[]
    for r in results; append!(all_rings_p, r["particle"]["rings"]); end
    
    pmf_g = Dict{Int, Float64}()
    for r in all_rings_g; pmf_g[r] = get(pmf_g, r, 0.0) + 1.0; end
    for (k,v) in pmf_g; pmf_g[k] = v / length(all_rings_g); end
    
    pmf_p = Dict{Int, Float64}()
    for r in all_rings_p; pmf_p[r] = get(pmf_p, r, 0.0) + 1.0; end
    for (k,v) in pmf_p; pmf_p[k] = v / length(all_rings_p); end
    
    open(joinpath(out_dir, "summary.json"), "w") do f
        JSON3.write(f, Dict("gillespie" => Dict("gel_stages" => mg, "pmf" => pmf_g), "particle" => Dict("gel_stages" => mp, "pmf" => pmf_p)))
    end
    println("Saved summary to $(joinpath(out_dir, "summary.json"))")
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CompareDSMC.main()
end
