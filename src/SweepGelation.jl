println("Initializing simulation modules (this may take a moment on first run)...")
include("PolymerUtils.jl")
include("DSMC.jl")
include("Network.jl")
include("Simulation.jl")
println("All modules loaded successfully.")

module SweepGelation

using ArgParse
using JSON3
using Statistics
using ProgressMeter
using ..SimulationMC
using ..PolymerUtils

function parse_commandline()
    s = ArgParseSettings()
    
    @add_arg_table! s begin
        "--nring_min"
            help = "Minimum ring length"
            arg_type = Int
            default = 256
        "--nring_max"
            help = "Maximum ring length"
            arg_type = Int
            default = 4096
        "--nring_step"
            help = "Step size for ring length"
            arg_type = Int
            default = 16
        "--nlin_min"
            help = "Minimum linear length"
            arg_type = Int
            default = 16
        "--nlin_max"
            help = "Maximum linear length"
            arg_type = Int
            default = 512
        "--nlin_step"
            help = "Step size for linear length"
            arg_type = Int
            default = 16
        "--ccsr"
            help = "Concentration scalar for rings"
            arg_type = Float64
            default = 5.0
        "--ccsl"
            help = "Concentration scalar for linears"
            arg_type = Float64
            default = 0.05
        "--L"
            help = "Box size"
            arg_type = Float64
            default = 200.0
        "--trials"
            help = "Number of trials per system"
            arg_type = Int
            default = 1000
        "--n_stages"
            help = "Number of stages"
            arg_type = Int
            default = 100
        "--k1"
            help = "Merge rate constant"
            arg_type = Float64
            default = SimulationMC.FITTED_K1_DEFAULT
        "--k2"
            help = "Cyclisation rate constant"
            arg_type = Float64
            default = SimulationMC.FITTED_K2_DEFAULT
        "--alpha"
            help = "Smoluchowski alpha parameter"
            arg_type = Float64
            default = 1.0
        "--nu"
            help = "Flory scaling exponent"
            arg_type = Float64
            default = 0.5
        "--val_A"
            help = "Valence model parameter A"
            arg_type = Float64
            default = SimulationMC.FITTED_A_DEFAULT
        "--max_steps"
            help = "Max steps per stage"
            arg_type = Int
            default = 1000000
        "--out_dir"
            help = "Output directory"
            arg_type = String
            default = PolymerUtils.RESULTS_DIR
        "--resume"
            help = "Resume from existing summary CSV"
            action = :store_true
    end
    
    return parse_args(s)
end

function build_system_grid(nrings, nlins, ccsr, ccsl, L)
    # Construct physically valid systems where both ring and linear counts are positive.
    systems = Dict{String, Any}[]
    for nring in nrings
        for nlin in nlins
            mring, mlin = calculate_polymer_numbers(ccsr, ccsl, L, RG_RING_128, RG_LINEAR_128, nring, nlin)
            if mring > 0 && mlin > 0
                push!(systems, Dict{String, Any}(
                    "L" => L,
                    "mring" => mring,
                    "nring" => nring,
                    "mlin" => mlin,
                    "nlin" => nlin
                ))
            end
        end
    end
    return systems
end

function main()
    args = parse_commandline()
    
    nrings = range(args["nring_min"], step=args["nring_step"], stop=args["nring_max"]-1)
    nlins = range(args["nlin_min"], step=args["nlin_step"], stop=args["nlin_max"]-1)
    
    systems = build_system_grid(nrings, nlins, args["ccsr"], args["ccsl"], args["L"])
    println("Built grid of $(length(systems)) systems.")
    
    mkpath(args["out_dir"])
    csv_path = joinpath(args["out_dir"], "sweep_summary.csv")
    
    completed_tags = Set{String}()
    if args["resume"] && isfile(csv_path)
        lines = readlines(csv_path)
        if length(lines) > 1
            for line in lines[2:end]
                parts = split(strip(line), ",")
                if length(parts) >= 5
                    try
                        L_val = parse(Float64, parts[1])
                        mring_val = round(Int, parse(Float64, parts[2]))
                        nring_val = round(Int, parse(Float64, parts[3]))
                        mlin_val = round(Int, parse(Float64, parts[4]))
                        nlin_val = round(Int, parse(Float64, parts[5]))
                        # Tag format must match the one used for fresh outputs.
                        tag = "L$(L_val)_mring$(mring_val)_nring$(nring_val)_mlin$(mlin_val)_nlin$(nlin_val)"
                        push!(completed_tags, tag)
                    catch
                        continue
                    end
                end
            end
        end
        println("Resuming: found $(length(completed_tags)) systems already completed.")
    elseif !isfile(csv_path)
        open(csv_path, "w") do f
            write(f, "L,mring,nring,mlin,nlin,mean_stages,std_stages,n_trials\n")
        end
    end
    
    sweep_results = Dict{String, Any}()
    
    println("Starting sweep over $(length(systems)) systems...")
    p = Progress(length(systems), dt=1.0, barglyphs=BarGlyphs("[=> ]"), barlen=40, color=:cyan)
    
    for (idx, sys_cfg) in enumerate(systems)
        tag = "L$(sys_cfg["L"])_mring$(sys_cfg["mring"])_nring$(sys_cfg["nring"])_mlin$(sys_cfg["mlin"])_nlin$(sys_cfg["nlin"])"
        
        if args["resume"] && in(tag, completed_tags)
            next!(p)
            continue
        end
        
        # Immediate feedback for SSH sessions where ProgressMeter bar might not render
        if idx % 10 == 0 || idx == 1
            println("Working on system $idx/$(length(systems)): $tag")
        end
        
        trial_results = run_trials_for_system(
            sys_cfg, 
            trials=args["trials"],
            n_stages=args["n_stages"],
            k1=args["k1"],
            k2=args["k2"],
            alpha=args["alpha"],
            nu=args["nu"],
            val_A=args["val_A"],
            max_steps=args["max_steps"]
        )
        
        stages_to_half = Float64[]
        for r in trial_results
            if r["stages_to_half"] !== nothing
                push!(stages_to_half, Float64(r["stages_to_half"]))
            end
        end
        
        if !isempty(stages_to_half)
            mean_s = mean(stages_to_half)
            std_s = length(stages_to_half) > 1 ? std(stages_to_half) : 0.0
        else
            mean_s = NaN
            std_s = NaN
        end
        
        sys_result = copy(sys_cfg)
        sys_result["mean_stages"] = mean_s
        sys_result["std_stages"] = std_s
        sys_result["n_trials_successful"] = length(trial_results)
        sys_result["n_gelled"] = length(stages_to_half)
        
        sweep_results[tag] = sys_result
        
        open(csv_path, "a") do f
            write(f, "$(sys_cfg["L"]),$(sys_cfg["mring"]),$(sys_cfg["nring"]),$(sys_cfg["mlin"]),$(sys_cfg["nlin"]),$(round(mean_s, digits=4)),$(round(std_s, digits=4)),$(length(stages_to_half))\n")
        end
        
        # Save trial results for this specific system
        results_all_json = joinpath(args["out_dir"], tag, "results_all.json")
        mkpath(joinpath(args["out_dir"], tag))
        open(results_all_json, "w") do f
            # wrap in system tag as root key if needed to match python dict output
            JSON3.write(f, Dict(tag => trial_results))
        end
        
        next!(p; showvalues = [(:system, tag)])
    end
    
    println("Sweep complete.")
end

end

# Entry point
if abspath(PROGRAM_FILE) == @__FILE__
    SweepGelation.main()
end
