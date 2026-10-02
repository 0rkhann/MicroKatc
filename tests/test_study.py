"""Loading and checking study files (spec sections 3 and 5).

python tests/test_study.py
"""

import copy
import os
import sys
import tempfile

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from study import load_study
from study_yaml import StudyError

TYPED = {
    "name": "Typed toy",
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


def write(doc, extra_files=()):
    folder = tempfile.mkdtemp()
    for name in extra_files:
        os.makedirs(os.path.join(folder, os.path.dirname(name)), exist_ok=True)
        open(os.path.join(folder, name), "w").close()
    path = os.path.join(folder, "study.yaml")
    with open(path, "w") as f:
        yaml.safe_dump(doc, f, allow_unicode=True)
    return path


def errors_of(doc, **kwargs):
    try:
        load_study(write(doc, **kwargs))
    except StudyError as error:
        return error.messages
    raise AssertionError("expected StudyError")


def changed(edit):
    doc = copy.deepcopy(TYPED)
    edit(doc)
    return doc


def assert_error(doc, text, **kwargs):
    messages = errors_of(doc, **kwargs)
    assert any(text in m for m in messages), (text, messages)


def test_typed_study_loads_with_derived_values():
    study = load_study(write(TYPED))
    assert study.temperature_K == 300.0 and isinstance(study.temperature_K, float)
    assert (
        len(study.studied_concentrations_M) == 3
        and study.studied_concentrations_M[0] == 0.01
    )
    assert abs(study.studied_concentrations_M[1] - 10**-1.5) < 1e-12
    assert study.species() == ["C1", "S", "C2", "P"]
    assert study.transition_states() == ["TS1"]
    assert study.energies_kcal_mol["TS1"] == {300.0: 15.0}
    assert study.reactions() == ["C1 + S = C2", "C2 = C1 + P"]
    assert study.cycles == {"cat": ("C1", "C2")}
    assert study.results_dir.name == "results"
    assert study.degree_of_rate_control is None and study.activation_energy is None


def test_max_product_uses_coefficients():
    doc = changed(lambda d: d["species"].update(overall_reaction="2 S + C1 <=> P"))
    doc["species"]["cycles"] = {"cat": ["C2"]}
    doc["analyses"]["microkinetics"]["initial_M"] = {"C2": 1e-3, "C1": 0.3}
    study = load_study(write(doc))
    assert study.max_product_M({"S": 1.0, "C1": 0.3}) == 0.3
    assert study.max_product_M({"S": 0.4, "C1": 0.3}) == 0.2


def test_unknown_key_suggests_the_closest():
    assert_error(
        changed(lambda d: d["conditions"].update(temprature_K=300)),
        "conditions.temprature_K: unknown key; did you mean temperature_K?",
    )


def test_files_and_energies_together():
    assert_error(
        changed(lambda d: d["species"].update(files="out")),
        "species: give files or energies_kcal_mol, not both",
    )


def test_required_keys():
    assert_error(
        changed(lambda d: d["species"].pop("product")), "species.product: required"
    )
    assert_error(
        changed(lambda d: d["species"].pop("overall_reaction")),
        "species.overall_reaction: required",
    )
    assert_error(changed(lambda d: d.pop("analyses")), "analyses: required")
    assert_error(
        changed(lambda d: d.update(analyses={})), "analyses: request at least one"
    )


def test_step_errors_are_located():
    assert_error(
        changed(lambda d: d["steps"].append("C2 + S C1")),
        'steps[2] "C2 + S C1": no <=> between the two sides',
    )
    assert_error(changed(lambda d: d["steps"].append("1.5 S <=> P")), "steps[2]")
    assert_error(
        changed(lambda d: d["steps"].append("C2 <=> P + C1 via TS1")),
        "steps[2] repeats steps[1]",
    )


def test_missing_energy_source():
    assert_error(
        changed(lambda d: d["species"]["energies_kcal_mol"].pop("TS1")),
        "no energy for TS1: add it to species.energies_kcal_mol",
    )
    doc = changed(lambda d: d["species"].pop("energies_kcal_mol"))
    doc["species"]["files"] = "out"
    os.environ["thermochange"] = tempfile.mkdtemp()
    messages = errors_of(
        doc, extra_files=["out/S.out", "out/P.out", "out/C1.out", "out/C2.out"]
    )
    assert "no energy for TS1: out/TS1.out not found" in messages, messages
    assert any(
        "does not contain formatters/formatted_energy_outputter.sh" in m
        for m in messages
    ), messages
    del os.environ["thermochange"]
    assert_error(doc, "$thermochange is not set", extra_files=["out/TS1.out"])
    # Scripts that only read finished results skip these checks
    from study import load_study as load

    study = load(write(doc), check_energy_sources=False)
    assert study.files.name == "out"


def test_name_consistency():
    assert_error(
        changed(lambda d: d["species"]["cycles"]["cat"].append("C9")),
        "species.cycles.cat: C9 does not appear in any step",
    )
    assert_error(
        changed(lambda d: d["species"]["cycles"].update(other=["C2"])),
        "species.cycles.other: C2 is also in cycle cat",
    )
    assert_error(
        changed(lambda d: d["species"].update(product="Q")),
        "species.product: Q does not appear in any step",
    )
    assert_error(
        changed(lambda d: d["species"].update(overall_reaction="S + X <=> P")),
        "species.overall_reaction: X does not appear in any step",
    )
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"].update(plot_species=["Z"])),
        "plot_species: Z does not appear in any step",
    )


def test_studied_species_and_product_not_in_initial_M():
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"]["initial_M"].update(S=0.1)),
        "the studied species is set by conditions.studied_range_M",
    )
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"]["initial_M"].update(P=0.0)),
        "the product always starts at 0",
    )
    assert_error(
        changed(
            lambda d: d["analyses"]["microkinetics"].update(
                initial_M={"S": 0.0, "P": 0.0}
            )
        ),
        "at least one catalyst intermediate",
    )


def test_typed_energies_need_every_activation_energy_temperature():
    def add_ea(d):
        d["analyses"]["activation_energy"] = {
            "initial_M": {"C1": 1e-6},
            "temperatures_K": {"from": 290, "to": 310, "points": 3},
            "simulation_time_s": 1000,
            "sampling_time_h": 0.1,
        }

    assert_error(
        changed(add_ea),
        "species.energies_kcal_mol.TS1: needs values at 290, 310 K for activation_energy",
    )


def test_single_values_enough_without_activation_energy():
    study = load_study(write(TYPED))
    assert study.temperatures() == [300.0]


def test_numbers():
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"].update(conversion=1.2)),
        "conversion: must be between 0 and 1",
    )
    assert_error(
        changed(lambda d: d["analyses"]["microkinetics"].update(catalyst_snapshot_h=1)),
        "1 h is after the end of the simulation (1000 s = 0.28 h)",
    )
    assert_error(
        changed(lambda d: d["conditions"].update(output_step_s=7)),
        "0.01 h is not a multiple of output_step_s (7 s)",
    )
    assert_error(
        changed(lambda d: d["conditions"].update(temperature_K=-5)),
        "conditions.temperature_K: must be positive",
    )
    assert_error(
        changed(
            lambda d: d["conditions"]["studied_range_M"].update({"from": 1, "to": 0.1})
        ),
        "from must be smaller than to",
    )


def test_unassigned_species_warning():
    doc = changed(lambda d: d["steps"].append("C2 + L <=> C3"))
    doc["species"]["energies_kcal_mol"].update(L=0.0, C3=-1.0)
    study = load_study(write(doc))
    assert study.warnings == (
        "warning: L, C3 belong to no cycle, so the catalyst distribution does not count them",
    ), study.warnings


def test_all_errors_reported_together():
    doc = changed(lambda d: d["conditions"].update(temprature_K=1))
    doc["species"]["product"] = "Q"
    messages = errors_of(doc)
    assert len(messages) >= 2, messages


def test_empty_analysis_is_an_error():
    assert_error(
        changed(lambda d: d["analyses"].update(microkinetics=None)),
        "analyses.microkinetics.simulation_time_s: must be a number, got None",
    )


def test_simulation_time_must_be_on_the_output_grid():
    assert_error(
        changed(lambda d: d["conditions"].update(output_step_s=3)),
        "analyses.microkinetics.simulation_time_s: 1000 s is not a multiple of output_step_s (3 s)",
    )


def test_overall_reactant_needs_an_initial_concentration():
    doc = changed(lambda d: d["species"].update(overall_reaction="S + B <=> P"))
    doc["steps"][1] = "C2 + B <=> C1 + P  via TS1"
    doc["species"]["energies_kcal_mol"]["B"] = 0.0
    assert_error(
        doc,
        "analyses.microkinetics.initial_M: give the overall reactant B, or no product can form",
    )


def test_reversed_and_empty_steps_are_errors():
    assert_error(
        changed(lambda d: d["steps"].append("C2 <=> S + C1")),
        "steps[2] repeats steps[0]",
    )
    assert_error(
        changed(lambda d: d["steps"].append("C2 <=> C2")),
        "steps[2] C2 <=> C2: both sides are the same",
    )


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
