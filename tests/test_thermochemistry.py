"""Gibbs energies and barrier tables from typed energies and from (fake) thermochange.

python tests/test_thermochemistry.py
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from study import load_study
from thermochemistry import (
    ThermochemistryError,
    barrier_table_name,
    write_barrier_tables,
)

TYPED = {
    "species": {
        "energies_kcal_mol": {"S": 0.0, "P": -10.0, "C1": 0.0, "C2": -2.0, "TS1": 15.0},
        "product": "P",
        "overall_reaction": "S <=> P",
        "cycles": {"cat": ["C1", "C2"]},
    },
    "steps": ["C1 + S <=> C2", "C2 <=> C1 + P  via TS1"],
    "conditions": {
        "temperature_K": 300,
        "studied_species": "S",
        "studied_range_M": {"from": 0.01, "to": 0.1, "points": 3},
    },
    "analyses": {
        "microkinetics": {
            "initial_M": {"C1": 1e-3},
            "simulation_time_s": 1000,
            "catalyst_snapshot_h": 0.01,
            "conversion": 0.99,
            "plot_species": ["C1", "C2"],
        }
    },
}

FAKE_THERMOCHANGE = """#!/bin/bash
echo "$@" >> "$FAKE_LOG"
[ "$1" = "-g" ] && shift
# like the real script, $1 is unquoted, so a path with a space breaks it
ls $1 > /dev/null 2>&1 || { echo "No valid output file!"; exit 1; }
echo leftover > temp_summary.temp
if [ "$FAKE_MODE" = zero ]; then
  printf "%s\\t-1.000000\\t-1.000000\\t0.000000\\n" "$1"
else
  printf "%s\\t-1.000000\\t-1.000000\\t-1.001000\\n" "$1"
fi
"""


def study_from(doc, files=()):
    folder = Path(tempfile.mkdtemp(prefix="study "))
    for name in files:
        (folder / name).parent.mkdir(parents=True, exist_ok=True)
        (folder / name).touch()
    (folder / "study.yaml").write_text(yaml.safe_dump(doc))
    return load_study(folder / "study.yaml")


def in_new_folder():
    os.chdir(tempfile.mkdtemp())


def test_typed_barriers_match_hand_calculation():
    study = study_from(TYPED)
    in_new_folder()
    messages = write_barrier_tables(study)
    table = pd.read_csv(Path("G_values_of_reactions") / barrier_table_name(300.0))
    # C1 + S <=> C2, barrierless: TS at max(0 + 0, -2) + 4 = 4
    # C2 <=> C1 + P via TS1: Gdir = 15 - (-2) = 17, Ginv = 15 - (0 - 10) = 25
    assert table["Rx"].tolist() == ["C1 + S = C2", "C2 = C1 + P"]
    assert table["TS"].tolist() == ["-", "TS1"]
    assert table["Gdir"].tolist() == [4.0, 17.0] and table["Ginv"].tolist() == [
        6.0,
        25.0,
    ], table
    assert "Overall reaction S <=> P: ΔG = -10.0 kcal/mol at 300.0 K" in messages, (
        messages
    )


def test_barrier_table_names_use_float_temperatures():
    study = study_from(TYPED)
    in_new_folder()
    write_barrier_tables(study)
    assert os.listdir("G_values_of_reactions") == [
        "reaction_df_300.0K_2.46172e+01atm.csv"
    ]


def test_coefficient_multiplies_the_energy():
    doc = {**TYPED, "steps": ["2 S <=> C2 via TS1", "C2 <=> C1 + P"]}
    doc["species"] = {
        **TYPED["species"],
        "overall_reaction": "2 S <=> P",
        "energies_kcal_mol": {"S": 1.0, "P": -10.0, "C1": 0.0, "C2": -3.0, "TS1": 10.0},
    }
    study = study_from(doc)
    in_new_folder()
    write_barrier_tables(study)
    table = pd.read_csv(Path("G_values_of_reactions") / barrier_table_name(300.0))
    assert table.loc[0, "Rx"] == "S + S = C2"
    assert table.loc[0, "Gdir"] == 10.0 - 2 * 1.0 and table.loc[0, "Ginv"] == 10.0 - (
        -3.0
    )


def test_negative_barrier_is_a_warning():
    doc = {
        **TYPED,
        "species": {
            **TYPED["species"],
            "energies_kcal_mol": {**TYPED["species"]["energies_kcal_mol"], "TS1": -5.0},
        },
    }
    study = study_from(doc)
    in_new_folder()
    messages = write_barrier_tables(study)
    assert any(
        "steps[1] (C2 <=> C1 + P via TS1) has a negative forward barrier (-3.0 kcal/mol)"
        in m
        for m in messages
    ), messages


def fake_thermochange(mode):
    root = Path(tempfile.mkdtemp())
    (root / "formatters").mkdir()
    script = root / "formatters" / "formatted_energy_outputter.sh"
    script.write_text(FAKE_THERMOCHANGE)
    os.environ.update(
        thermochange=str(root), FAKE_MODE=mode, FAKE_LOG=str(root / "calls.log")
    )
    return root


def files_study(correction="RRHO"):
    doc = {
        **TYPED,
        "species": {
            k: v for k, v in TYPED["species"].items() if k != "energies_kcal_mol"
        },
    }
    doc["species"].update(files="out", vibrational_correction=correction)
    return study_from(
        doc, files=[f"out/{n}.out" for n in ("S", "P", "C1", "C2", "TS1")]
    )


def test_file_energies_are_converted_and_leave_no_temporary_files():
    root = fake_thermochange("ok")
    study = files_study()
    in_new_folder()
    write_barrier_tables(study)
    compounds = pd.read_csv(
        Path("G_values_of_compounds") / "G_values_at_300.0K_2.46172e+01atm.csv",
        index_col=0,
    )
    assert abs(compounds.loc["S", "Gibbs Free Energies"] - (-1.001 * 627.509)) < 1e-9
    assert sorted(os.listdir(".")) == [
        "G_values_of_compounds",
        "G_values_of_reactions",
    ], os.listdir(".")
    calls = (root / "calls.log").read_text().splitlines()
    expected = f" 300.0 {0.082057366080960 * 300.0}"  # the same strings the old shell pipeline passed
    assert len(calls) == 5 and all(c.endswith(expected) for c in calls), calls
    assert all(
        c.split()[0] in {f"{n}.out" for n in ("S", "P", "C1", "C2", "TS1")}
        for c in calls
    ), calls  # only the file name: the study folder's path has a space


def test_grimme_passes_the_g_flag():
    root = fake_thermochange("ok")
    study = files_study("Grimme")
    in_new_folder()
    write_barrier_tables(study)
    assert all(
        c.startswith("-g ") for c in (root / "calls.log").read_text().splitlines()
    )


def test_zero_energy_from_thermochange_stops_the_run():
    fake_thermochange("zero")
    study = files_study()
    in_new_folder()
    try:
        write_barrier_tables(study)
    except ThermochemistryError as error:
        assert "corrected G of 0.0 hartree" in str(error), str(error)
    else:
        raise AssertionError("expected ThermochemistryError")
    assert not Path("G_values_of_reactions").exists() or not os.listdir(
        "G_values_of_reactions"
    )


def test_rerun_reports_the_same_messages():
    doc = {
        **TYPED,
        "species": {
            **TYPED["species"],
            "energies_kcal_mol": {**TYPED["species"]["energies_kcal_mol"], "TS1": -5.0},
        },
    }
    study = study_from(doc)
    in_new_folder()
    first = write_barrier_tables(study)
    assert write_barrier_tables(study) == first and len(first) == 2, first


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
