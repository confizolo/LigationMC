println("Initializing PolymerUtils...")
include("PolymerUtils.jl")
println("PolymerUtils loaded.")

module ScalingMrNr

using ArgParse
using Statistics
using Printf
using ..PolymerUtils

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--nr_min"
            help = "Minimum ring length"
            arg_type = Int
            default = 64
        "--nr_max"
            help = "Maximum ring length"
            arg_type = Int
            default = 4096
        "--nr_step"
            help = "Step for ring length"
            arg_type = Int
            default = 16
        "--L"
            help = "Reference box length"
            arg_type = Float64
            default = 80.0
        "--ccsr"
            help = "Ring concentration multiplier (c/c*)"
            arg_type = Float64
            default = 5.0
        "--ccsl"
            help = "Linear concentration multiplier (for ml output only)"
            arg_type = Float64
            default = 0.05
        "--nlin"
            help = "Linear length used when computing ml"
            arg_type = Int
            default = 64
        "--out_csv"
            help = "Output CSV path"
            arg_type = String
            default = "./mr_nr_scaling.csv"
    end

    return parse_args(s)
end

# Fit y = C * x^beta by linear regression in log-space.
function power_law_exponent(xs::Vector{Float64}, ys::Vector{Float64})
    lx = log.(xs)
    ly = log.(ys)

    mx = mean(lx)
    my = mean(ly)

    denom = sum((x - mx)^2 for x in lx)
    if denom == 0.0
        error("Cannot fit exponent: all x values are identical.")
    end

    beta = sum((lx[i] - mx) * (ly[i] - my) for i in eachindex(lx)) / denom
    intercept = my - beta * mx
    return beta, exp(intercept)
end

function main()
    args = parse_commandline()

    nrs = collect(args["nr_min"]:args["nr_step"]:args["nr_max"])
    if isempty(nrs)
        error("No ring lengths selected. Check nr_min, nr_max, nr_step.")
    end

    rows = NamedTuple{(:nr, :mr, :ml, :mr_nr), Tuple{Int, Int, Int, Float64}}[]
    xs = Float64[]
    ys = Float64[]

    for nr in nrs
        mr, ml = calculate_polymer_numbers(
            args["ccsr"],
            args["ccsl"],
            args["L"],
            RG_RING_128,
            RG_LINEAR_128,
            nr,
            args["nlin"]
        )

        if mr > 0
            mr_nr = Float64(mr * nr)
            push!(rows, (nr=nr, mr=mr, ml=ml, mr_nr=mr_nr))
            push!(xs, Float64(nr))
            push!(ys, mr_nr)
        end
    end

    if isempty(rows)
        error("No valid rows (mr > 0). Increase L or ccsr, or lower nr range.")
    end

    beta, coeff = power_law_exponent(xs, ys)

    out_dir = dirname(args["out_csv"])
    if !isempty(out_dir)
        mkpath(out_dir)
    end

    open(args["out_csv"], "w") do f
        write(f, "nr,mr,ml,mr_nr\n")
        for r in rows
            write(f, "$(r.nr),$(r.mr),$(r.ml),$(r.mr_nr)\n")
        end
    end

    @printf("Saved scaling table to %s\n", args["out_csv"])
    @printf("Fitted power law: mr*nr = %.6g * nr^(%.6f)\n", coeff, beta)
    println("Expected asymptotic exponent from c* construction: -0.5")
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    ScalingMrNr.main()
end
