include("../src/PolymerUtils.jl")
include("../src/GillespieSSA.jl")
include("../src/ParticleDSMC.jl")
include("../src/Network.jl")
include("../src/MDMapping.jl")
include("../src/CompareDSMC.jl")

println("Running smoke test for DSMC methods...")
empty!(ARGS)
push!(ARGS, "--L", "80.0", "--mring", "10", "--nring", "128", "--mlin", "4", "--nlin", "64", "--trials", "5", "--n_stages", "20", "--out_dir", "./results/compare/smoke")

CompareDSMC.main()
println("Smoke test completed.")
