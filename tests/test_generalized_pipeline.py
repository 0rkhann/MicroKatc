"""The analyses no longer assume a product named prod, regex-named cycles or a 1:1 yield.

Uses a fake copasi_parser, so it runs without COPASI:
    python tests/test_generalized_pipeline.py
"""

import os
import sys
import tempfile
import types

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.modules.setdefault("copasi_parser", types.ModuleType("copasi_parser"))

from auxiliary_functions import AuxiliaryFunctions
from plotting_functions import PlotFunctions


def test_conversion_threshold_is_per_concentration():
    time = np.arange(0, 11) / 10
    sims = [
        pd.DataFrame({"time": time, "P": time * 1.0}),
        pd.DataFrame({"time": time, "P": time * 2.0}),
    ]
    # First simulation must pass 0.5, second must pass 0.5 too: thresholds differ per concentration
    times = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
        sims, [0.01, 0.1], "S", [0.45, 1.5], "P"
    )
    assert times == [(0.5, 0.01), (0.8, 0.1)], times


def test_grid_sizes():
    assert PlotFunctions.grid(1) == (1, 1)
    assert PlotFunctions.grid(5) == (1, 5)
    assert PlotFunctions.grid(6) == (2, 5)
    assert PlotFunctions.grid(22) == (5, 5)
    assert PlotFunctions.grid(0) == (1, 1)


def test_results_follow_the_current_folder():
    from file_operations import FileOperations, run_directory

    folder = tempfile.mkdtemp()
    os.chdir(folder)
    assert run_directory() == os.getcwd()
    open("table.csv", "w").close()
    FileOperations.move_to_output_directory("out", "table.csv")
    assert os.path.isfile(os.path.join(folder, "out", "table.csv"))


def test_drc_uses_the_product_and_time_step():
    calls = []

    class FakeModel:
        def self_destruct(self):
            pass

    cpx = sys.modules["copasi_parser"]
    cpx.prepare_copasi_model = lambda **kwargs: FakeModel()
    cpx.drc_calc = lambda base_model, **kwargs: (
        calls.append(kwargs) or (np.array([1.0]), 1.0)
    )
    os.chdir(tempfile.mkdtemp())
    from apparent_activation_energy import DRCAnalysis

    DRCAnalysis(
        np.array([1e-3]),
        np.array([300.0]),
        100,
        {"C1": 1e-3},
        "S",
        ["A = B"],
        0.1,
        1,
        product="P",
        time_step=5,
    ).calculate_degree_of_rate_control()
    assert {c["target_spc"] for c in calls} == {"P"} and {
        c["time_step"] for c in calls
    } == {5}, calls


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
