"""Builds the README summary figures in pics/ from the results that main.py saves.

Run after main.py, from the same directory:
    python readme_figures.py
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from auxiliary_functions import AuxiliaryFunctions
from microkinetics_simulation import SIMULATIONS_OUTPUT_DIR_NAME, SimulationHandler
from plotting_functions import PlotFunctions

REACTANT = "PMe3"
T = 350.0
RDS = ["I8_0L = I9_0L", "I3_1L = I4_1L"]
DRC_STEPS = ["I3_0L = I4_0L", "I8_0L = I9_0L", "I3_1L = I4_1L"]
POISONING = ["I1_0L", "I7_0L", "I1_1L", "I7_1L"]
# Same as c0_1, total_simulation_time_Ea and time in main.py (apparent Ea conditions)
C0_EA = {"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 1e-6}
EA_SIMULATION_TIME = 10_000
EA_TIME_H = 2


def latest(kind):
    """Path of the newest df_<kind> table that main.py saved"""
    paths = glob.glob(os.path.join(SIMULATIONS_OUTPUT_DIR_NAME, f"df_{kind}_*.csv"))
    if not paths:
        raise FileNotFoundError(
            f"No df_{kind} table in {SIMULATIONS_OUTPUT_DIR_NAME}/: run main.py first"
        )
    return max(paths, key=os.path.getmtime)


def main():
    df_flux, df_rate, df_drc = (pd.read_csv(latest(k)) for k in ("flux", "rate", "drc"))
    df_flux["reaction"] = df_flux["reaction"].str.split().str.join(" ")
    df_drc = df_drc.rename(columns=lambda c: " ".join(c.split()))
    concentrations = np.sort(df_drc[f"c0({REACTANT})"].unique())
    plots = PlotFunctions()

    plots.plot_Ea_vs_reactant_c0(
        1, 2, df_flux, (15, 10), REACTANT, RDS, "", "ri", log_x=True
    )
    plt.savefig("pics/Ea_of_rate_determining_steps.png", dpi=200, bbox_inches="tight")
    plots.plot_drc_vs_c0(1, 3, df_drc, (15, 10), DRC_STEPS, REACTANT, "", T, log_x=True)
    plt.savefig("pics/DRC_of_influencing_steps.png", dpi=200, bbox_inches="tight")

    # Combined view: poisoning intermediates, apparent Ea and DRC against c0(PMe3)
    handler = SimulationHandler(
        T, AuxiliaryFunctions.compute_pressure_value(T), REACTANT, EA_SIMULATION_TIME
    )
    at_ea_time = [
        handler.get_simulation_df({**C0_EA, REACTANT: c})
        .set_index("time")
        .loc[EA_TIME_H]
        for c in concentrations
    ]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
    for species in POISONING:
        axes[0].plot(
            concentrations, [row[species] for row in at_ea_time], "o-", label=species
        )
    axes[0].set(
        yscale="log",
        ylabel="concentration (M)",
        title=f"Poisoning intermediates at t = {EA_TIME_H} h",
    )
    for reaction in RDS:
        rows = df_flux[df_flux["reaction"] == reaction].sort_values(f"c0({REACTANT})")
        axes[1].plot(rows[f"c0({REACTANT})"], rows["Ea"], "o-", label=reaction)
    product = df_rate[df_rate["compound"] == "prod"].sort_values(f"c0({REACTANT})")
    axes[1].plot(
        product[f"c0({REACTANT})"], product["Ea"], "^--", label="product formation"
    )
    axes[1].set(
        ylabel="apparent Ea (kcal/mol)", title="Apparent activation energy, 325-375 K"
    )
    for reaction in DRC_STEPS:
        axes[2].plot(
            concentrations,
            df_drc.sort_values(f"c0({REACTANT})")[reaction],
            "o-",
            label=reaction,
        )
    axes[2].set(ylabel="DRC", title=f"Degree of rate control, T = {T:.0f} K")
    for ax in axes:
        ax.set(xscale="log", xlabel=f"c0({REACTANT}) (M)")
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig("pics/final_pieces.png", dpi=200)
    print(
        "Saved pics/Ea_of_rate_determining_steps.png, pics/DRC_of_influencing_steps.png, pics/final_pieces.png"
    )


if __name__ == "__main__":
    main()
