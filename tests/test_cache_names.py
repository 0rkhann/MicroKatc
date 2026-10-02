"""Cache file names change whenever any input that affects their contents changes.

Uses a fake copasi_parser, so it runs without COPASI:
    python tests/test_cache_names.py   (or: pytest tests)
"""

import os
import sys
import tempfile
import types

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def test_cache_names_depend_on_all_inputs():
    calls = []

    class FakeModel:
        def self_destruct(self):
            pass

    def time_course_simulation(ch, total_time, outfile, **kwargs):
        calls.append(outfile)
        pd.DataFrame({"time": [0.0, total_time], "prod": [0.0, 1.0]}).to_csv(
            outfile, index=False
        )
        return None, True

    cpx = types.ModuleType("copasi_parser")
    cpx.prepare_copasi_model = lambda **kwargs: FakeModel()
    cpx.time_course_simulation = time_course_simulation
    cpx.read_simulation = pd.read_csv
    sys.modules["copasi_parser"] = cpx

    os.chdir(tempfile.mkdtemp())
    from apparent_activation_energy import ApparentEaAnalysis, DRCAnalysis
    from microkinetics_simulation import SimulationHandler

    # simulation files: a different time step is a different file
    c0 = {"PMe3": 1e-5, "CO": 1}
    SimulationHandler(350.0, 28.7, "PMe3", 3600, time_step=1).get_simulation_df(c0)
    SimulationHandler(350.0, 28.7, "PMe3", 3600, time_step=1).get_simulation_df(c0)
    assert len(calls) == 1  # same inputs: cached
    SimulationHandler(350.0, 28.7, "PMe3", 3600, time_step=2).get_simulation_df(c0)
    assert len(calls) == 2  # new time step: recomputed

    # flux / rate tables
    base = {
        "temperature_value": 350.0,
        "T_values_array": [300.0, 350.0, 400.0],
        "reactant_concentration_array": [1e-5, 1e-4, 1e-3],
        "reactant_to_study": "PMe3",
        "c0": {"PMe3": 1e-5, "CO": 1},
        "total_simulation_time": 3600,
        "time": 1,
        "reactions": ["r1", "r2"],
        "time_step": 1,
    }

    def names(**changes):
        a = ApparentEaAnalysis(**{**base, **changes})
        return a.df_flux_filename, a.df_rate_filename

    ref = names()
    assert names() == ref  # stable
    assert ref[0].startswith("df_flux_T_range_300.0K_400.0K_C(PMe3)_1e-05M_0.001M_")
    for changes in (
        {"T_values_array": [300.0, 325.0, 350.0, 375.0, 400.0]},
        {"reactant_concentration_array": [1e-5, 1e-3]},
        {"c0": {"PMe3": 1e-5, "CO": 2}},
        {"total_simulation_time": 7200},
        {"time": 2},
        {"time_step": 2},
    ):
        new = names(**changes)
        assert new[0] != ref[0] and new[1] != ref[1], changes

    # DRC table
    drc = {
        "reactant_concentration_array": base["reactant_concentration_array"],
        "T_values_array": base["T_values_array"],
        "total_simulation_time": 3600,
        "c0": base["c0"],
        "reactant_to_study": "PMe3",
        "reactions": ["r1", "r2"],
        "e_shift": 0.1,
    }
    ref = DRCAnalysis(**drc).df_drc_filename
    assert DRCAnalysis(**drc).df_drc_filename == ref
    for changes in (
        {"T_values_array": [300.0, 350.0, 375.0, 400.0]},
        {"reactant_concentration_array": [1e-5, 1e-3]},
        {"c0": {"PMe3": 1e-5, "CO": 2}},
        {"total_simulation_time": 7200},
        {"e_shift": 0.2},
    ):
        assert DRCAnalysis(**{**drc, **changes}).df_drc_filename != ref, changes


if __name__ == "__main__":
    test_cache_names_depend_on_all_inputs()
    print("ok")
