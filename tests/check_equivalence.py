"""Compares a run of the example study with the baseline taken before the YAML change.

    python tests/check_equivalence.py examples/hydroformylation/results

Tolerances come from COPASI's own run-to-run noise, measured over the 10 pairs of 5 fresh runs on
2026-10-02 (two of them of the old code, the others of the new or transitional code; two runs of the
identical old code differ as much as old and new runs do). Each limit is about 1.5 times the largest
difference seen, where that difference is not far below a round number:
- barrier tables: deterministic, compared exactly;
- Ea of the rate-determining steps (r08, r13): largest difference 6.3e-5 kcal/mol, limit 1e-3;
- Ea of product formation up to 1e-2 M PMe3: largest difference 0.020 kcal/mol, limit 0.03. Above
  1e-2 M it varies more: two GitHub runs of the same commit gave 0.081 and 0.008 at 3.2e-2 M
  (2026-10-02), limit 0.12;
- DRC of the three key steps up to 1e-2 M PMe3: moved 5.5e-4, limit 1e-3. Above 1e-2 M less than 1 %
  of the catalyst is in the 0L cycle and the finite difference is ill-conditioned: the same input
  gave 0.9655, 0.9666 or 0.9875 for I3_1L = I4_1L at 0.1 M depending on key order, number type and
  what ran before in the process. There the largest difference was 0.022, limit 0.03. For every
  step, including near-equilibrium ones: largest difference 0.041, limit 0.06;
- catalyst amounts at 1 h: moved 3.6e-6 relative, limit 1e-4; times to 99 %: identical, limit 1e-6.
Near-equilibrium steps (net flux a small difference of large rates) and which near-zero rows are
skipped for sign changes vary from run to run and are not compared.
"""

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "hydroformylation_baseline"
)
KEY_DRC_STEPS = ["I8_0L = I9_0L", "I3_1L = I4_1L", "I3_0L = I4_0L"]
failures = []


def same_text(series):
    return series.astype(str).str.split().str.join(" ")


def report(name, worst, limit):
    ok = worst <= limit
    print(
        f"{'ok  ' if ok else 'FAIL'} {name}: max difference {worst:.3g} (limit {limit:g})"
    )
    if not ok:
        failures.append(name)


def merged(results, kind, names):
    old = pd.read_csv(os.path.join(BASE, f"df_{kind}.csv"))
    new = pd.read_csv(
        glob.glob(
            os.path.join(results, "microkinetics_simulations", f"df_{kind}_*.csv")
        )[0]
    )
    column = next(c for c in old.columns if c.startswith("c0("))
    # The concentrations come from a log spacing whose last bit differs between machines
    # (1e-10 here, 9.999999999999999e-11 on GitHub's runners), so match them rounded
    for df in (old, new):
        df[column] = df[column].map(lambda c: float(f"{c:.12g}"))
    both = old.merge(new, on=["name", column], suffixes=("_old", "_new"))
    both = both[both["name"].isin(names)]
    expected = old[old["name"].isin(names)]
    assert len(both) == len(expected), (
        f"df_{kind}: {names} rows missing from the new run"
    )
    return both


def main(results):
    for base_file in sorted(glob.glob(os.path.join(BASE, "reaction_df_*.csv"))):
        name = os.path.basename(base_file)
        new = pd.read_csv(os.path.join(results, "G_values_of_reactions", name))
        old = pd.read_csv(base_file)
        assert (same_text(new["Rx"]) == same_text(old["Rx"])).all(), (
            f"{name}: steps differ"
        )
        difference = np.abs(
            new[["Gdir", "Ginv"]].to_numpy() - old[["Gdir", "Ginv"]].to_numpy()
        ).max()
        report(f"barriers {name}", float(difference), 1e-9)

    both = merged(results, "flux", ["r08.Flux", "r13.Flux"])
    report(
        "Ea of the rate-determining steps",
        float((both["Ea_old"] - both["Ea_new"]).abs().max()),
        1e-3,
    )
    both = merged(results, "rate", ["prod.Rate"])
    difference = (both["Ea_old"] - both["Ea_new"]).abs()
    low = both[next(c for c in both.columns if c.startswith("c0("))] <= 1e-2
    report("Ea of product formation up to 1e-2 M", float(difference[low].max()), 3e-2)
    report(
        "Ea of product formation above 1e-2 M",
        float(difference[~low].max()),
        0.12,
    )

    def drc(path):
        return pd.read_csv(path).rename(columns=lambda c: " ".join(c.split()))

    old = drc(os.path.join(BASE, "df_drc.csv"))
    new = drc(
        glob.glob(os.path.join(results, "microkinetics_simulations", "df_drc_*.csv"))[0]
    )
    steps = [c for c in old.columns if "=" in c]
    conditioned = old[next(c for c in old.columns if c.startswith("c0("))] <= 1e-2
    key = np.abs(old[KEY_DRC_STEPS].to_numpy() - new[KEY_DRC_STEPS].to_numpy())
    report(
        "DRC of the key steps up to 1e-2 M",
        float(key[conditioned.to_numpy()].max()),
        1e-3,
    )
    report(
        "DRC of the key steps above 1e-2 M",
        float(key[~conditioned.to_numpy()].max(initial=0)),
        3e-2,
    )
    report(
        "DRC of every step",
        float(np.abs(old[steps].to_numpy() - new[steps].to_numpy()).max()),
        6e-2,
    )

    with open(os.path.join(BASE, "microkinetics.json")) as f:
        old = json.load(f)
    with open(os.path.join(results, "microkinetics.json")) as f:
        new = json.load(f)
    for cycle in old["catalyst_M"]:
        a, b = np.array(old["catalyst_M"][cycle]), np.array(new["catalyst_M"][cycle])
        report(
            f"catalyst in {cycle} at 1 h",
            float(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30))),
            1e-4,
        )
    a, b = np.array(old["t99_h"]), np.array(new["t99_h"])
    report("time to 99 % conversion", float(np.max(np.abs(a - b) / a)), 1e-6)

    if failures:
        sys.exit(f"{len(failures)} comparison(s) failed: {', '.join(failures)}")
    print("Run matches the baseline.")


if __name__ == "__main__":
    main(sys.argv[1])
