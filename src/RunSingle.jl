println("Initializing simulation modules (this may take a moment on first run)...")
include("PolymerUtils.jl")
include("GillespieSSA.jl")
include("Network.jl")
include("Simulation.jl")
println("All modules loaded successfully.")

module RunSingle

using ArgParse
using JSON3
using Statistics
using ..SimulationMC
using ..PolymerUtils

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--L"
            help = "Box side length"
            arg_type = Float64
            default = 200.0
        "--mring"
            help = "Number of ring polymers"
            arg_type = Int
            required = true
        "--nring"
            help = "Monomers per ring polymer"
            arg_type = Int
            required = true
        "--mlin"
            help = "Number of linear polymers"
            arg_type = Int
            required = true
        "--nlin"
            help = "Monomers per linear polymer"
            arg_type = Int
            required = true
        "--trials"
            help = "Number of independent trials"
            arg_type = Int
            default = 1000
        "--n_stages"
            help = "Number of reaction stages per trial"
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
            help = "Max DSMC steps per stage"
            arg_type = Int
            default = 1000000
        "--out_dir"
            help = "Output directory for results"
            arg_type = String
            default = PolymerUtils.RESULTS_DIR
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    sys_cfg = Dict{String, Any}(
        "L"     => args["L"],
        "mring" => args["mring"],
        "nring" => args["nring"],
        "mlin"  => args["mlin"],
        "nlin"  => args["nlin"]
    )

    tag = "L$(sys_cfg["L"])_mring$(sys_cfg["mring"])_nring$(sys_cfg["nring"])_mlin$(sys_cfg["mlin"])_nlin$(sys_cfg["nlin"])"
    println("Running single system: $tag")
    println("  trials=$(args["trials"]), n_stages=$(args["n_stages"])")
    println("  k1=$(args["k1"]), k2=$(args["k2"]), alpha=$(args["alpha"]), nu=$(args["nu"]), val_A=$(args["val_A"])")

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

    # ── Aggregate stages_to_half statistics ──────────────────────────
    stages_to_half = Float64[]
    for r in trial_results
        if r["stages_to_half"] !== nothing
            push!(stages_to_half, Float64(r["stages_to_half"]))
        end
    end

    if !isempty(stages_to_half)
        mean_s = mean(stages_to_half)
        std_s  = length(stages_to_half) > 1 ? std(stages_to_half) : 0.0
    else
        mean_s = NaN
        std_s  = NaN
    end

    println("Results: mean_stages=$(round(mean_s, digits=4)), std_stages=$(round(std_s, digits=4)), gelled=$(length(stages_to_half))/$(length(trial_results))")

    # ── Save outputs ─────────────────────────────────────────────────
    out_dir = joinpath(args["out_dir"], tag)
    mkpath(out_dir)

    # Full trial-level results
    results_json_path = joinpath(out_dir, "results_all.json")
    open(results_json_path, "w") do f
        JSON3.write(f, Dict(tag => trial_results))
    end
    println("Saved trial results to $results_json_path")

    # Summary CSV
    summary_csv_path = joinpath(out_dir, "summary.csv")
    open(summary_csv_path, "w") do f
        write(f, "L,mring,nring,mlin,nlin,mean_stages,std_stages,n_trials,n_gelled\n")
        write(f, "$(sys_cfg["L"]),$(sys_cfg["mring"]),$(sys_cfg["nring"]),$(sys_cfg["mlin"]),$(sys_cfg["nlin"]),$(round(mean_s, digits=4)),$(round(std_s, digits=4)),$(length(trial_results)),$(length(stages_to_half))\n")
    end
    println("Saved summary to $summary_csv_path")
end

end

# Entry point
if abspath(PROGRAM_FILE) == @__FILE__
    RunSingle.main()
end
