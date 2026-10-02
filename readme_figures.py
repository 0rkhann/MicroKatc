"""Builds the README figures in pics/ from the results that main.py saves.

Each figure is sized for GitHub's README column (about 840 px wide) and shares one
style: colour identifies the catalytic cycle (0L blue, 1L pink, as in the paper's TOC
graphic) and the product is green with a dashed line and triangles. main.py's own
figures, with every step and every concentration, stay in microkinetics_simulations_images/.

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

REACTANT = "PMe3"
T = 350.0
# Same as c0_1 / c0_2 and the simulation times in main.py, so the cached simulations are reused
C0_EA = {"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 1e-6}
C0_MK = {"CO": 0.05, "H2": 0.05, "ete": 0.05, "prod": 0, "I1_0L": 5e-4}
EA_SIMULATION_TIME, EA_TIME_H = 10_000, 2
MK_SIMULATION_TIME, CATALYST_TIME_H, CONVERSION = 100_000, 1, 0.99

CYCLE_0L, CYCLE_1L, PRODUCT, INK, MUTED = (
    "#2a78d6",
    "#e87ba4",
    "#1baf7a",
    "#2b2b2b",
    "#8a8a85",
)
RDS_0L, RDS_1L, INHIBITOR = "I8_0L = I9_0L", "I3_1L = I4_1L", "I3_0L = I4_0L"
C0_LABEL = r"initial [PMe$_3$] (M)"
EA_LABEL = r"apparent $E_a$ (kcal mol$^{-1}$)"

plt.rcParams.update(
    {
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "legend.fontsize": 10.5,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#e6e6e3",
        "grid.linewidth": 0.8,
        "lines.linewidth": 2,
        "lines.markersize": 6,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    }
)


def step_label(reaction):
    """'I8_0L = I9_0L' -> 'I8-0L ⇌ I9-0L'"""
    return reaction.replace("_", "-").replace(" = ", " ⇌ ")


def molar(c):
    """1e-06 -> '10$^{-6}$ M' (half decades as 3 x 10$^{-4}$ M)"""
    mantissa, exponent = f"{c:.0e}".split("e")
    exponent = int(exponent)
    return (
        f"$10^{{{exponent}}}$ M"
        if mantissa == "1"
        else f"{mantissa} × $10^{{{exponent}}}$ M"
    )


def latest(kind):
    """Newest df_<kind> table that main.py saved"""
    paths = glob.glob(os.path.join(SIMULATIONS_OUTPUT_DIR_NAME, f"df_{kind}_*.csv"))
    if not paths:
        raise FileNotFoundError(
            f"No df_{kind} table in {SIMULATIONS_OUTPUT_DIR_NAME}/: run main.py first"
        )
    return pd.read_csv(max(paths, key=os.path.getmtime))


def save(fig, name):
    fig.savefig(os.path.join("pics", name), dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved pics/{name}")


def log_x(ax):
    ax.set_xscale("log")
    ax.set_xlabel(C0_LABEL)


def main():
    flux, rate, drc = latest("flux"), latest("rate"), latest("drc")
    flux["reaction"] = flux["reaction"].str.split().str.join(" ")
    drc = drc.rename(columns=lambda c: " ".join(c.split())).sort_values(
        f"c0({REACTANT})"
    )
    c0 = drc[f"c0({REACTANT})"].to_numpy()
    by_c0 = lambda df: df.sort_values(f"c0({REACTANT})")
    product = by_c0(rate[rate["compound"] == "prod"])
    rds = {r: by_c0(flux[flux["reaction"] == r]) for r in (RDS_0L, RDS_1L)}

    # 1. Arrhenius check in each regime: ln(r / r at 350 K) against 1000/T
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, conc, step, colour, regime in (
        (axes[0], c0.min(), RDS_0L, CYCLE_0L, "0L regime"),
        (axes[1], c0.max(), RDS_1L, CYCLE_1L, "1L regime"),
    ):
        for row, colour_, style, label in (
            (
                rds[step][np.isclose(rds[step][f"c0({REACTANT})"], conc)].iloc[0],
                colour,
                "o-",
                step_label(step),
            ),
            (
                product[np.isclose(product[f"c0({REACTANT})"], conc)].iloc[0],
                PRODUCT,
                "^--",
                "product formation",
            ),
        ):
            x = 1000 * np.array([row[f"1/T({i})"] for i in range(5)])
            y = np.array([row[f"ln_value({i})"] for i in range(5)])
            y_ref = row["slope"] * (1 / T) + row["intercept"]
            ax.plot(
                x,
                y - y_ref,
                style,
                color=colour_,
                label=f"{label}: $E_a$ = {row['Ea']:.1f} kcal mol$^{{-1}}$",
            )
        ax.set_title(f"{regime}, [PMe$_3$] = {molar(conc)}", color=INK)
        ax.set_xlabel(r"1000 / $T$ (K$^{-1}$)")
        ax.legend(loc="upper right")
    axes[0].set_ylabel(r"ln($r$ / $r_{350\,\mathrm{K}}$)")
    fig.tight_layout()
    save(fig, "arrhenius.png")

    # 2. Apparent Ea against [PMe3]
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    for step, colour in ((RDS_0L, CYCLE_0L), (RDS_1L, CYCLE_1L)):
        ax.plot(
            rds[step][f"c0({REACTANT})"],
            rds[step]["Ea"],
            "o-",
            color=colour,
            label=f"{step_label(step)} (rate-determining, {step[-2:]})",
        )
    ax.plot(
        product[f"c0({REACTANT})"],
        product["Ea"],
        "^--",
        color=PRODUCT,
        label="product formation",
    )
    ax.set_ylabel(EA_LABEL)
    ax.set_title("Apparent activation energy (325–375 K)", color=INK)
    log_x(ax)
    ax.legend(loc="upper left")
    save(fig, "apparent_ea.png")

    # 3. Degree of rate control against [PMe3]
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.axhline(0, color=MUTED, linewidth=1)
    for step, colour, style, note in (
        (RDS_0L, CYCLE_0L, "o-", "rate-determining, 0L"),
        (INHIBITOR, CYCLE_0L, "s--", "inhibiting, 0L"),
        (RDS_1L, CYCLE_1L, "o-", "rate-determining, 1L"),
    ):
        ax.plot(
            c0, drc[step], style, color=colour, label=f"{step_label(step)} ({note})"
        )
    ax.set_ylabel("degree of rate control")
    ax.set_title(f"Degree of rate control at {T:.0f} K", color=INK)
    log_x(ax)
    ax.legend(loc="center left")
    save(fig, "drc.png")

    # 4-6. Catalyst distribution, concentration profiles and conversion time (main.py's c0_2)
    mk = SimulationHandler(
        T, AuxiliaryFunctions.compute_pressure_value(T), REACTANT, MK_SIMULATION_TIME
    )
    sims = [mk.get_simulation_df({**C0_MK, REACTANT: c}) for c in c0]
    cycles = AuxiliaryFunctions.find_intermediates_of_cycle(["0L", "1L"])
    catalyst = AuxiliaryFunctions.get_concentrations_of_catalyst(
        sims, CATALYST_TIME_H, cycles
    )
    total = np.add(*catalyst)

    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    for share, colour, label in zip(
        catalyst, (CYCLE_0L, CYCLE_1L), ("0L cycle", "1L cycle")
    ):
        ax.plot(c0, 100 * np.array(share) / total, "o-", color=colour, label=label)
    ax.set_ylabel(f"share of catalyst at t = {CATALYST_TIME_H} h (%)")
    ax.set_ylim(-3, 103)
    ax.set_title("Catalyst distribution between the cycles", color=INK)
    log_x(ax)
    ax.legend(loc="center left")
    save(fig, "catalyst_distribution.png")

    shown = [np.argmin(abs(np.log10(c0) - v)) for v in (-6, -3.5, -1)]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.9), sharey=True)
    for ax, i in zip(axes, shown):
        sim = sims[i][sims[i]["time"] > 0]
        for species, colour, style in (
            ("I1_0L", CYCLE_0L, "-"),
            ("I7_0L", CYCLE_0L, "--"),
            ("I1_1L", CYCLE_1L, "-"),
            ("I7_1L", CYCLE_1L, "--"),
        ):
            ax.plot(
                sim["time"],
                sim[species],
                style,
                color=colour,
                label=species.replace("_", "-"),
            )
        ax.set(xscale="log", yscale="log", xlabel="time (h)", ylim=(1e-12, 1e-3))
        ax.set_title(f"[PMe$_3$] = {molar(c0[i])}", color=INK)
    axes[0].set_ylabel("concentration (M)")
    fig.tight_layout()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 1.0))
    save(fig, "concentration_profiles.png")

    times = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
        sims, c0, REACTANT, CONVERSION * min(C0_MK[r] for r in ("CO", "H2", "ete"))
    )
    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    ax.plot(
        [c for _, c in times],
        [t for t, _ in times],
        "^--",
        color=PRODUCT,
        label="time to 99 % conversion",
    )
    ax.set_ylabel(f"time to {CONVERSION:.0%} conversion (h)")
    ax.set_ylim(0, None)
    ax.set_title(f"Time to {CONVERSION:.0%} conversion", color=INK)
    log_x(ax)
    save(fig, "conversion_time.png")

    # 7. Poisoning intermediates at the Ea sampling time, against [PMe3] (main.py's c0_1)
    ea = SimulationHandler(
        T, AuxiliaryFunctions.compute_pressure_value(T), REACTANT, EA_SIMULATION_TIME
    )
    rows = [
        ea.get_simulation_df({**C0_EA, REACTANT: c}).set_index("time").loc[EA_TIME_H]
        for c in c0
    ]
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    for species, colour, style in (
        ("I1_0L", CYCLE_0L, "o-"),
        ("I7_0L", CYCLE_0L, "s--"),
        ("I1_1L", CYCLE_1L, "o-"),
        ("I7_1L", CYCLE_1L, "s--"),
    ):
        ax.plot(
            c0,
            [row[species] for row in rows],
            style,
            color=colour,
            label=species.replace("_", "-"),
        )
    ax.set(yscale="log", ylabel="concentration (M)")
    ax.set_title(f"Poisoning intermediates at t = {EA_TIME_H} h", color=INK)
    log_x(ax)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    save(fig, "poisoning_intermediates.png")


if __name__ == "__main__":
    main()
