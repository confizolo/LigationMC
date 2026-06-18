using JSON3
using Statistics

include("simulation_jl/PolymerUtils.jl")
include("simulation_jl/DSMC.jl")
include("simulation_jl/Network.jl")
include("simulation_jl/Main.jl")

using .MainLigMC

function run_julia_test()
    sys_cfg = Dict{String, Any}(
        "L" => 80.0,
        "mring" => 27,
        "nring" => 1024,
        "mlin" => 6,
        "nlin" => 64
    )
    
    results = MainLigMC.run_trials_for_system(sys_cfg, trials=100)
    
    stages = Float64[]
    for r in results
        if r["stages_to_half"] !== nothing
            push!(stages, Float64(r["stages_to_half"]))
        end
    end
    
    mean_stages = isempty(stages) ? NaN : mean(stages)
    println("Julia mean stages: ", mean_stages)
end

run_julia_test()
