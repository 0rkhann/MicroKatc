"""DRC uses a central difference: +e_shift and -e_shift runs are averaged, cancelling first-order error.

Uses a fake copasi_parser, so it runs without COPASI:
    python tests/test_drc_central_difference.py   (or: pytest tests)
"""

import os
import sys
import tempfile
import types

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

TRUE_DRC = np.array([0.7, 0.4, -0.1])
FIRST_ORDER = np.array([2.0, -1.0, 0.5])  # one-sided error per kcal/mol of shift


def test_central_difference_cancels_first_order_error():
    shifts = []

    class FakeModel:
        def self_destruct(self):
            pass

    def drc_calc(base_model, e_shift, **kwargs):
        shifts.append(e_shift)
        return TRUE_DRC + FIRST_ORDER * e_shift, 1.0

    cpx = types.ModuleType("copasi_parser")
    cpx.prepare_copasi_model = lambda **kwargs: FakeModel()
    cpx.drc_calc = drc_calc
    sys.modules["copasi_parser"] = cpx

    os.chdir(tempfile.mkdtemp())
    from apparent_activation_energy import DRCAnalysis

    reactions = ["a = b", "b = c", "c = d"]
    drc = DRCAnalysis(
        np.array([1e-4]), np.array([350.0]), 10_000, {"X": 0}, "X", reactions, 0.1
    ).calculate_degree_of_rate_control()

    assert sorted(shifts) == [-0.1, 0.1]
    assert np.allclose(drc[reactions].to_numpy()[0], TRUE_DRC)


if __name__ == "__main__":
    test_central_difference_cancels_first_order_error()
    print("ok")
