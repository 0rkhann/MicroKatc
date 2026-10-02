"""Checks a full main.py run against the published results (Abdullayev et al., ACS Catal. 2025, 15, 4739).

Not a unit test: it reads the results of the example study, so run it after the example:
    python main.py && python tests/check_reproduction.py

Paper values come from the text and from reading Figures 4-6; tolerances cover that reading
precision and the ~0.08 DRC gap documented in the README.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from auxiliary_functions import AuxiliaryFunctions
from microkatc import initial_concentrations
from microkinetics_simulation import SimulationHandler
from readme_figures import (
    CONVERSION,
    MK,
    MK_SIMULATION_TIME,
    REACTANT,
    STUDY,
    T,
    latest,
)

os.chdir(STUDY.results_dir)

failures = []


def check(name, value, low, high):
    ok = low <= value <= high
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {value:.3f} (expected {low} to {high})")
    if not ok:
        failures.append(name)


def at(df, column, c0):
    row = df[np.isclose(df[f"c0({REACTANT})"], c0, rtol=1e-6, atol=0)]
    return row[column].item()


flux, rate, drc = latest("flux"), latest("rate"), latest("drc")
flux["reaction"] = flux["reaction"].str.split().str.join(" ")
drc = drc.rename(columns=lambda c: " ".join(c.split()))
c0 = np.sort(drc[f"c0({REACTANT})"].unique())
low, high = c0.min(), c0.max()
product = rate[rate["compound"] == "prod"]

# Figure 6: apparent activation energies (kcal/mol)
check(
    "Ea of product formation, low PMe3 (paper 23.45)",
    at(product, "Ea", low),
    23.3,
    23.6,
)
check(
    "Ea of product formation, 0.1 M PMe3 (paper 22.05)",
    at(product, "Ea", high),
    21.9,
    22.2,
)
check("minimum Ea of product formation (paper 21.35)", product["Ea"].min(), 21.2, 21.7)
for step, paper in (("I8_0L = I9_0L", 23.5), ("I3_1L = I4_1L", 17.5)):
    rows = flux[flux["reaction"] == step]
    check(
        f"Ea of {step}, low PMe3 (paper {paper})",
        at(rows, "Ea", low),
        paper - 0.15,
        paper + 0.15,
    )

# Figure 5: degree of rate control at 350 K
inhibitor = drc.set_index(f"c0({REACTANT})")["I3_0L = I4_0L"]
check("minimum DRC of I3_0L = I4_0L (paper -0.22)", inhibitor.min(), -0.26, -0.18)
check(
    "log10 c0 of that minimum (paper -3.5)", np.log10(inhibitor.idxmin()), -4.01, -2.99
)
check(
    "DRC of I8_0L = I9_0L, low PMe3 (paper ~1)",
    at(drc, "I8_0L = I9_0L", low),
    0.97,
    1.01,
)

# Figure 4: catalyst distribution at 1 h and time to 99 % conversion
handler = SimulationHandler(
    T,
    AuxiliaryFunctions.compute_pressure_value(T),
    REACTANT,
    MK_SIMULATION_TIME,
    STUDY.output_step_s,
)
sims = [
    handler.get_simulation_df(initial_concentrations(STUDY, MK.initial_M, c))
    for c in c0
]
cycles = {label: list(members) for label, members in STUDY.cycles.items()}
zero_l, one_l = AuxiliaryFunctions.get_concentrations_of_catalyst(sims, 1, cycles)
crossing = c0[np.argmax(np.array(one_l) > np.array(zero_l))]
check(
    "log10 c0 where 1L holds most catalyst (paper ~-3.2)",
    np.log10(crossing),
    -3.6,
    -2.9,
)
times = (
    AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
        sims,
        c0,
        REACTANT,
        [
            CONVERSION
            * STUDY.max_product_M(initial_concentrations(STUDY, MK.initial_M, c))
            for c in c0
        ],
        STUDY.product,
    )
)
check("time to 99 % conversion, low PMe3 / h (paper 5.66)", times[0][0], 5.55, 5.75)
check("time to 99 % conversion, 0.1 M PMe3 / h (paper 1.69)", times[-1][0], 1.6, 1.8)

if failures:
    sys.exit(f"{len(failures)} result(s) differ from the paper: {', '.join(failures)}")
print("All results match the paper.")
