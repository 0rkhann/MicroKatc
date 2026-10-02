"""Fluxes are attached to reactions by COPASI name (rNN = row NN); a mismatch must raise, not mislabel.

Uses a fake copasi_parser, so it runs without COPASI:
    python tests/test_reaction_matching.py   (or: pytest tests)
"""

import os
import sys
import types

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.modules["copasi_parser"] = types.ModuleType("copasi_parser")


def _flux_df(n_fluxes, order=None):
    rows = [
        {
            "name": f"r{k:02d}.Flux",
            "c0(A)": 1.0,
            "1/T": 1 / T,
            "ln_value": -k * 1000 / T,
        }
        for k in (order or range(1, n_fluxes + 1))
        for T in (300.0, 310.0, 320.0)
    ]
    return pd.DataFrame(rows)


def _calc(df, reactions):
    from apparent_activation_energy import ReactionParameterCalculator

    return ReactionParameterCalculator.calculate_reaction_parameters(
        df, [300.0, 310.0, 320.0], np.array([1.0]), "A", "ri", reactions=reactions
    )


def test_matching_counts():
    out = _calc(_flux_df(2), ["a = b", "b = c"])
    assert list(out["reaction"]) == ["a = b", "b = c"]


def test_fluxes_out_of_order_keep_their_reaction():
    out = _calc(_flux_df(2, order=[2, 1]), ["a = b", "b = c"])
    assert dict(zip(out["name"], out["reaction"])) == {
        "r01.Flux": "a = b",
        "r02.Flux": "b = c",
    }


def test_flux_count_mismatch_raises():
    try:
        _calc(_flux_df(3), ["a = b", "b = c"])
    except ValueError as e:
        assert "3" in str(e) and "2" in str(e)
    else:
        raise AssertionError("expected ValueError")


if __name__ == "__main__":
    test_matching_counts()
    test_fluxes_out_of_order_keep_their_reaction()
    test_flux_count_mismatch_raises()
    print("ok")
