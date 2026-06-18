println("Initializing simulation modules (first run may compile)...")
include("PolymerUtils.jl")
include("DSMC.jl")
include("Network.jl")
include("Main.jl")
println("Modules loaded.")

module RunSingle

using ArgParse
using JSON3
using Statistics
using Dates
using ..MainLigMC

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--L"
            help = "Reference box size"
            arg_type = Float64
            required = true
        "--mring"
            help = "Initial number of rings"
            arg_type = Int
            required = true
        "--nring"
            help = "Ring polymer length"
            arg_type = Int
            required = true
        "--mlin"
            help = "Linear chains injected per stage"
            arg_type = Int
            required = true
        "--nlin"
            help = "Linear chain length"
            arg_type = Int
            required = true
        "--trials"
            help = "Number of independent Monte Carlo trials"
            arg_type = Int
            default = 1000
        "--n_stages"
            help = "Number of growth stages"
            arg_type = Int
            default = 100
        "--k1"
            help = "Merge rate constant"
            arg_type = Float64
            default = MainLigMC.FITTED_K1_DEFAULT
        "--k2"
            help = "Cyclisation rate constant"
            arg_type = Float64
            default = MainLigMC.FITTED_K2_DEFAULT
        "--alpha"
            help = "Smoluchowski alpha parameter"
            arg_type = Float64
            default = 1.0
        "--nu"
            help = "Flory scaling exponent"
            arg_type = Float64
            default = 0.5
        "--val_A"
            help = "Valence model prefactor A"
            arg_type = Float64
            default = MainLigMC.FITTED_A_DEFAULT
        "--max_steps"
            help = "Maximum SSA events per stage"
            arg_type = Int
            default = 1000000
        "--out_dir"
            help = "Directory for outputs"
            arg_type = String
            default = "./results"
    end

    return parse_args(s)
end

function make_tag(args)
    return "L$(args[\"L\"])_mring$(args[\"mring\"])_nring$(args[\"nring\"])_mlin$(args[\"mlin\"])_nlin$(args[\"nlin\"])"
end

function summarize(trial_results)
    stages_to_half = Float64[]
    for r in trial_results
        if r["stages_to_half"] !== nothing
            push!(stages_to_half, Float64(r["stages_to_half"]))
        end
    end

    mean_stages = isempty(stages_to_half) ? NaN : mean(stages_to_half)
    std_stages = length(stages_to_half) > 1 ? std(stages_to_half) : 0.0

    return Dict(
        "mean_stages" => mean_stages,
        "std_stages" => std_stages,
        "n_trials" => length(trial_results),
        "n_gelled" => length(stages_to_half)
    )
end

function main()
    args = parse_commandline()

    sys_cfg = Dict{String, Any}(
        "L" => args["L"],
        "mring" => args["mring"],
        "nring" => args["nring"],
        "mlin" => args["mlin"],
        "nlin" => args["nlin"]
    )

    println("Running system: ", sys_cfg)
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

    summary = summarize(trial_results)

    tag = make_tag(args)
    out_system_dir = joinpath(args["out_dir"], tag)
    mkpath(out_system_dir)

    results_path = joinpath(out_system_dir, "results_all.json")
    open(results_path, "w") do f
        JSON3.write(f, Dict(tag => trial_results))
    end

    summary_path = joinpath(out_system_dir, "summary.csv")
    open(summary_path, "w") do f
        write(f, "timestamp,L,mring,nring,mlin,nlin,mean_stages,std_stages,n_trials,n_gelled\n")
        write(
            f,
            string(
                Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"), ",",
                args["L"], ",",
                args["mring"], ",",
                args["nring"], ",",
                args["mlin"], ",",
                args["nlin"], ",",
                summary["mean_stages"], ",",
                summary["std_stages"], ",",
                summary["n_trials"], ",",
                summary["n_gelled"],
                "\n"
            )
        )
    end

    println("Saved: ", results_path)
    println("Saved: ", summary_path)
    println("Summary: ", summary)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    RunSingle.main()
end
