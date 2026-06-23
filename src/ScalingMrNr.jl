println("Loading PolymerUtils...")
include("PolymerUtils.jl")
println("PolymerUtils loaded.")

module ScalingMrNr

using ArgParse
using ..PolymerUtils

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--nr_min"
            help = "Minimum ring monomer count"
            arg_type = Int
            default = 64
        "--nr_max"
            help = "Maximum ring monomer count"
            arg_type = Int
            default = 4096
        "--nr_step"
            help = "Step size for ring monomer count"
            arg_type = Int
            default = 16
        "--L"
            help = "Box side length"
            arg_type = Float64
            default = 80.0
        "--ccsr"
            help = "Concentration scalar for rings (c/c*)"
            arg_type = Float64
            default = 5.0
        "--ccsl"
            help = "Concentration scalar for linears (c/c*)"
            arg_type = Float64
            default = 0.05
    end

    return parse_args(s)
end

"""
Simple least-squares fit of a degree-1 polynomial (y = a + b*x) via the normal equations.
Returns (intercept, slope).
"""
function polyfit1(x::Vector{Float64}, y::Vector{Float64})::Tuple{Float64, Float64}
    n = length(x)
    sx  = sum(x)
    sy  = sum(y)
    sxx = sum(x .* x)
    sxy = sum(x .* y)
    denom = n * sxx - sx^2
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    return (a, b)
end

function main()
    args = parse_commandline()

    nrs = collect(range(args["nr_min"], step=args["nr_step"], stop=args["nr_max"]))
    L   = args["L"]
    ccsr = args["ccsr"]
    ccsl = args["ccsl"]

    # ── Compute mr and mr*nr for each Nr ────────────────────────────
    nr_vec    = Int[]
    mr_vec    = Int[]
    mr_nr_vec = Int[]

    for nr in nrs
        (mr, _) = calculate_polymer_numbers(ccsr, ccsl, L, RG_RING_128, RG_LINEAR_128, nr, REFERENCE_N)
        if mr > 0
            push!(nr_vec, nr)
            push!(mr_vec, mr)
            push!(mr_nr_vec, mr * nr)
        end
    end

    if length(nr_vec) < 2
        println("ERROR: fewer than 2 valid data points — cannot fit power law.")
        return
    end

    # ── Log-log linear regression:  log(mr*nr) = log(A) + B*log(nr) ─
    log_nr    = log.(Float64.(nr_vec))
    log_mr_nr = log.(Float64.(mr_nr_vec))

    (log_A, B) = polyfit1(log_nr, log_mr_nr)
    A_coeff = exp(log_A)

    println("Power-law fit: mr*nr = $(round(A_coeff, digits=6)) * nr^$(round(B, digits=6))")
    println("Fitted exponent B = $(round(B, digits=6))")

    # ── Print CSV to stdout ──────────────────────────────────────────
    println()
    println("nr,mr,mr_nr")
    for i in eachindex(nr_vec)
        println("$(nr_vec[i]),$(mr_vec[i]),$(mr_nr_vec[i])")
    end
end

end

# Entry point
if abspath(PROGRAM_FILE) == @__FILE__
    ScalingMrNr.main()
end
