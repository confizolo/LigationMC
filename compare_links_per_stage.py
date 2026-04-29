#!/usr/bin/env python3
"""Compare observed links-per-stage with valence-model predictions.

Expect a CSV with one row per (stage,size) group containing at least:
 - columns for `nring`, `mring`, `nlin` (or `linear_size`),
 - `n_events` (number of cyclisation events in the group),
 - `observed_links` (total links observed for those events),
 - optionally `targets` (number of potential existing rings considered per event).

The script reads a fitted model JSON produced by `fit_valence_model.py` containing `A`.
It computes predicted links-per-event using mu = A * (basis) where basis can be
`nring*mring*nlin` or `nring*nlin` (choose via --basis). For each row the predicted
mean links-per-event is `targets * (1 - exp(-mu))` where `targets` is taken from the
`targets` column if present, otherwise defaults to `mring`.

Outputs a CSV with observed and predicted columns and an optional PNG plot.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def detect_cols(df: pd.DataFrame) -> dict:
    return {
        'nring': 'nring' if 'nring' in df.columns else None,
        'mring': 'mring' if 'mring' in df.columns else None,
        'nlin': 'nlin' if 'nlin' in df.columns else ('linear_size' if 'linear_size' in df.columns else None),
        'n_events': 'n_events' if 'n_events' in df.columns else ('n_cyclisations' if 'n_cyclisations' in df.columns else None),
        'observed_links': 'observed_links' if 'observed_links' in df.columns else ('links' if 'links' in df.columns else None),
        'targets': 'targets' if 'targets' in df.columns else None,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--csv', required=True, help='Input summary CSV (grouped by stage/size)')
    p.add_argument('--model', required=True, help='Fitted model JSON with A (output from fit_valence_model)')
    p.add_argument('--basis', choices=['nring*mring*nlin', 'nring*nlin'], default='nring*mring*nlin')
    p.add_argument('--box-volume', type=float, default=1.0)
    p.add_argument('--out-csv', help='Output CSV with predictions', default='compare_links_per_stage_out.csv')
    p.add_argument('--out-png', help='Optional output PNG (scatter observed vs predicted)')
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    with open(args.model, 'r') as fh:
        model = json.load(fh)
    A = float(model.get('A'))

    cols = detect_cols(df)
    missing = [k for k, v in cols.items() if v is None and k in ('nring', 'nlin', 'n_events', 'observed_links')]
    if missing:
        raise SystemExit(f'Missing required columns in CSV: {missing}. Detected columns: {df.columns.tolist()}')

    # compute derived columns
    nring = df[cols['nring']]
    nlin = df[cols['nlin']]
    mring = df[cols['mring']] if cols['mring'] is not None else pd.Series(np.ones(len(df)), index=df.index)
    targets = df[cols['targets']] if cols['targets'] is not None else mring
    n_events = df[cols['n_events']]
    observed_links = df[cols['observed_links']]

    if args.basis == 'nring*mring*nlin':
        basis = nring * mring * nlin / float(args.box_volume)
    else:
        basis = nring * nlin / float(args.box_volume)

    mu = A * basis
    p_link = 1.0 - np.exp(-mu)
    predicted_links_per_event = targets * p_link

    df_out = df.copy()
    df_out['mu'] = mu
    df_out['predicted_p_link'] = p_link
    df_out['predicted_links_per_event'] = predicted_links_per_event
    df_out['observed_links_per_event'] = observed_links / n_events.replace({0: np.nan})

    outpath = Path(args.out_csv)
    df_out.to_csv(outpath, index=False)
    print(f'Wrote predictions to {outpath}')

    if args.out_png:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 6))
        x = df_out['observed_links_per_event']
        y = df_out['predicted_links_per_event']
        plt.scatter(x, y)
        lim = max(np.nanmax(x), np.nanmax(y))
        plt.plot([0, lim], [0, lim], 'k--', alpha=0.6)
        plt.xlabel('Observed links per event')
        plt.ylabel('Predicted links per event')
        plt.title(f'Observed vs Predicted links per event (A={A:.3g})')
        plt.tight_layout()
        plt.savefig(args.out_png)
        print(f'Wrote plot to {args.out_png}')


if __name__ == '__main__':
    main()
