"""Fit the Valence Model linking parameter A to MD data.

Minimises the Root Mean Squared Error (RMSE) between the theoretical
Poisson expected linking and the empirically observed MD link counts.

Outputs a JSON file with the fitted {A}.
"""

import json
import math
import os
import pandas as pd
from scipy.optimize import minimize_scalar

# Deterministic mapping for the 16 reference systems (L = 80.0)
NRING_TO_MRING = {
    256: 222,
    512: 78,
    768: 42,
    1024: 27
}

def main():
    md_csv = (
        "/storage/cmstore02/groups/TAPLab/fconforto-projects/"
        "fconforto-olympic-gels/results/histories/"
        "summary_all_systems_links_by_size_s0.csv"
    )
    md_df = pd.read_csv(md_csv)
    
    L = 80.0
    
    def objective(A: float) -> float:
        if A <= 0:
            return 1e6
        
        squared_errors = []
        for _, row in md_df.iterrows():
            nring = int(row["nring"])
            l_cyc = int(row["linear_size"])
            avg_links_created = float(row["avg_links_created"])
            
            if nring not in NRING_TO_MRING:
                continue
            
            # The fitting restricts cyclised sizes to <= 300 to avoid
            # extreme high-variance outliers at the tail end of the distribution.
            if l_cyc > 300:
                continue
                
            mring = NRING_TO_MRING[nring]
            
            # Expected threadings/links per target
            mu_target = A * nring * l_cyc / (L ** 3)
            # Total expected valid linked targets
            predicted_mean = mring * (1.0 - math.exp(-mu_target))
            
            error = (predicted_mean - avg_links_created) ** 2
            squared_errors.append(error)
            
        if not squared_errors:
            return 1e6
            
        rmse = math.sqrt(sum(squared_errors) / len(squared_errors))
        return rmse

    print("Starting optimization for parameter A...")
    res = minimize_scalar(
        objective, bounds=(0.05, 0.5), method="bounded"
    )

    print("\nOptimization Complete:")
    print("Best A:", res.x)
    print("Best RMSE:", res.fun)

    out_dict = {"A": res.x, "rmse": res.fun, "method": "PoissonLinkingRMSE"}

    os.makedirs("./parameters", exist_ok=True)
    out_path = "./parameters/fitted_valence_model.json"
    with open(out_path, "w") as f:
        json.dump(out_dict, f, indent=2)
    print(f"Saved fitted parameters to {out_path}")

if __name__ == "__main__":
    main()
