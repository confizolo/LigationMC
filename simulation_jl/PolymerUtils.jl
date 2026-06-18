module PolymerUtils

export RESULTS_DIR, REFERENCE_N, RG_RING_128, RG_LINEAR_128
export calculate_polymer_numbers, smoluchowski_kernel, cyclisation_rate, valence_model

const RESULTS_DIR = "/storage/cmstore02/groups/TAPLab/fconforto-projects/fconforto-olympic-gels-mc"

# Default Rg values from MD equilibration of 128-monomer chains.
const REFERENCE_N = 128
const RG_RING_128 = 9.9     # Rg of a 128-monomer ring polymer [σ]
const RG_LINEAR_128 = 14.0  # Rg of a 128-monomer linear polymer [σ]

"""
Compute `(mr, ml)`: ring and linear counts from c/c* constraints.

The geometric argument is:
1. Scale reference `Rg` as `sqrt(N/128)`.
2. Approximate each polymer as an occupied sphere of volume `(4/3)πRg^3`.
3. Compute overlap concentration counts from `V_box / V_polymer`.
4. Multiply by requested concentration factors `ccsr` and `ccsl`.
"""
function calculate_polymer_numbers(
    ccsr::Float64,
    ccsl::Float64,
    L::Float64 = 200.0,
    Rgr_base::Float64 = RG_RING_128,
    Rgl_base::Float64 = RG_LINEAR_128,
    Nr::Int = 1024,
    Nl::Int = 128
)::Tuple{Int, Int}
    rgr_scaled = Rgr_base * sqrt(Nr / 128.0)
    rgl_scaled = Rgl_base * sqrt(Nl / 128.0)

    v_box = L^3
    v_polymer_ring = (4.0 / 3.0) * pi * (rgr_scaled^3)
    v_polymer_linear = (4.0 / 3.0) * pi * (rgl_scaled^3)

    mrc = v_box / v_polymer_ring
    mlc = v_box / v_polymer_linear

    mr = floor(Int, mrc * ccsr)
    ml = floor(Int, mlc * ccsl)
    return (mr, ml)
end

function smoluchowski_kernel(i::Int, j::Int; alpha::Float64=1.0, nu::Float64=0.5)::Float64
    # Pair kernel used for merge propensity a_ij = k1 * n_i * n_j * K(i,j).
    return (Float64(i)^(-alpha) + Float64(j)^(-alpha)) * (Float64(i)^nu + Float64(j)^nu)
end

function cyclisation_rate(length::Int, k2::Float64; nu::Float64=0.5)::Float64
    # Larger chains cyclise more slowly through the power-law scaling exponent nu.
    return k2 * (Float64(length)^(-4.0 * nu))
end

function valence_model(l_cyc::Int, n_total::Int, A::Float64, box_volume::Float64=1.0)::Float64
    if l_cyc <= 0 || n_total <= 0 || box_volume <= 0.0
        return 0.0
    end
    concentration_proxy = (Float64(n_total) * Float64(l_cyc)) / box_volume
    return A * concentration_proxy
end

end
