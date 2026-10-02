"""Gibbs energies of a study's species and the barrier tables the COPASI simulations read"""

import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from auxiliary_functions import AuxiliaryFunctions
from calculating_G_for_microkinetics import (
    G_COMPOUNDS_OUTPUT_DIR_NAME,
    REACTION_DF_OUTPUT_DIR_NAME,
    GibbsEnergyCalculator,
)

HARTREE_TO_KCAL_MOL = 627.509
# A temperature/pressure correction moves G by millihartrees; more than this means thermochange failed
MAX_CORRECTION_HARTREE = 0.1


class ThermochemistryError(RuntimeError):
    """thermochange failed or gave an implausible Gibbs energy"""


def barrier_table_name(temperature_K):
    """File name the simulations look up for the barriers at temperature_K"""
    pressure = AuxiliaryFunctions.compute_pressure_value(temperature_K)
    return f"reaction_df_{temperature_K}K_{pressure:.5e}atm.csv"


def compounds_table_name(temperature_K):
    pressure = AuxiliaryFunctions.compute_pressure_value(temperature_K)
    return f"G_values_at_{temperature_K}K_{pressure:.5e}atm.csv"


def gibbs_from_file(out_file, temperature_K, correction, thermochange):
    """G (kcal/mol, 1 M) of one output file at temperature_K, through thermochange"""
    pressure = AuxiliaryFunctions.compute_pressure_value(temperature_K)
    script = Path(thermochange) / "formatters" / "formatted_energy_outputter.sh"
    flags = ["-g"] if correction == "Grimme" else []
    name = Path(out_file).stem
    command = [
        "bash",
        str(script),
        *flags,
        f"{name}.out",
        str(temperature_K),
        str(pressure),
    ]
    # thermochange writes temp_summary.temp to the current folder and does not always remove it.
    # It also uses its file argument unquoted, so it gets a bare name linked into that folder: the
    # study's own path may contain spaces.
    with tempfile.TemporaryDirectory() as scratch:
        os.symlink(Path(out_file).resolve(), Path(scratch) / f"{name}.out")
        result = subprocess.run(
            command, cwd=scratch, capture_output=True, text=True, check=False
        )
    lines = result.stdout.strip().splitlines()
    fields = lines[-1].split("\t") if lines else []
    try:
        file_g, corrected = float(fields[2]), float(fields[3])
    except (IndexError, ValueError) as error:
        raise ThermochemistryError(
            f"thermochange gave no Gibbs energy for {name}:\n{result.stdout}{result.stderr}"
        ) from error
    if corrected == 0.0 or abs(corrected - file_g) > MAX_CORRECTION_HARTREE:
        raise ThermochemistryError(
            f"thermochange gave a corrected G of {corrected} hartree for {name} "
            f"(the file's own G is {file_g}); its error output:\n{result.stderr}"
        )
    return corrected * HARTREE_TO_KCAL_MOL


def _typed(study, name, temperature_K):
    for known, value in study.energies_kcal_mol[name].items():
        if abs(known - temperature_K) < 1e-6:
            return value
    raise ThermochemistryError(f"no typed energy for {name} at {temperature_K} K")


def gibbs_energies(study, temperature_K):
    """{species or transition state: G in kcal/mol} at temperature_K"""
    names = study.species() + [
        ts for ts in study.transition_states() if ts not in study.species()
    ]
    if study.files is None:
        return {name: _typed(study, name, temperature_K) for name in names}
    thermochange = os.environ["thermochange"]
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 1) as pool:
        values = list(
            pool.map(
                lambda name: gibbs_from_file(
                    study.files / f"{name}.out",
                    temperature_K,
                    study.vibrational_correction,
                    thermochange,
                ),
                names,
            )
        )
    return dict(zip(names, values))


def barrier_table(study, gibbs):
    """Rx, TS, Gdir, Ginv for every step, in the format copasi_helper reads"""
    table = pd.DataFrame(
        {"Rx": study.reactions(), "TS": [step.ts or "-" for step in study.steps]}
    )
    energies = pd.DataFrame({"Gibbs Free Energies": pd.Series(gibbs)})
    gdir, ginv = GibbsEnergyCalculator(
        energies
    ).calculate_direct_inverse_reactions_gibbs_free_energies(table)
    table["Gdir"], table["Ginv"] = gdir, ginv
    return table


def negative_barrier_warnings(study, table):
    messages = []
    for i, (step, gdir, ginv) in enumerate(
        zip(study.steps, table["Gdir"], table["Ginv"])
    ):
        for value, direction in ((gdir, "forward"), (ginv, "reverse")):
            if value < 0:
                messages.append(
                    f"warning: steps[{i}] ({step.text}) has a negative {direction} barrier ({value:.1f} kcal/mol)"
                )
    return messages


def overall_reaction_energy(study, gibbs):
    """ΔG of the overall reaction from the species energies (kcal/mol)"""
    products = sum(n * gibbs[name] for name, n in study.overall.products)
    reactants = sum(n * gibbs[name] for name, n in study.overall.reactants)
    return products - reactants


def write_barrier_tables(study):
    """Writes the compound and barrier tables for every needed temperature into the current folder.

    Tables already present are kept (the results folder is per study, guarded by run_info.json).
    Returns the messages to show: warnings and the overall reaction energy at the working temperature.
    """
    messages = []
    for temperature in study.temperatures():
        reactions_path = Path(REACTION_DF_OUTPUT_DIR_NAME) / barrier_table_name(
            temperature
        )
        compounds_path = Path(G_COMPOUNDS_OUTPUT_DIR_NAME) / compounds_table_name(
            temperature
        )
        if reactions_path.is_file() and compounds_path.is_file():
            continue
        print(f"Computing Gibbs energies at {temperature} K...")
        gibbs = gibbs_energies(study, temperature)
        table = barrier_table(study, gibbs)
        compounds_path.parent.mkdir(exist_ok=True)
        reactions_path.parent.mkdir(exist_ok=True)
        pd.DataFrame({"Gibbs Free Energies": pd.Series(gibbs)}).rename_axis(
            "Compounds"
        ).to_csv(compounds_path)
        partial = reactions_path.with_suffix(".partial")
        table.to_csv(partial, index=False)
        partial.replace(
            reactions_path
        )  # written last and atomically: its presence means complete
        messages += negative_barrier_warnings(study, table)
        if temperature == study.temperature_K:
            delta = overall_reaction_energy(study, gibbs)
            messages.append(
                f"Overall reaction {study.overall.text}: ΔG = {delta:.1f} kcal/mol at {temperature} K"
            )
    return messages
