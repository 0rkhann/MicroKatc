"""Apparent Ea fit: zero / sign-changing fluxes are skipped, signed fluxes keep their sign, R2 survives flat data.

Uses a fake copasi_parser, so it runs without COPASI:
    python tests/test_ea_flux_sign.py   (or: pytest tests)
"""

import os
import sys
import types

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def test_flux_sign_handling():
    sys.modules["copasi_parser"] = types.ModuleType("copasi_parser")
    from apparent_activation_energy import (
        ReactionDataHandler,
        ReactionParameterCalculator,
    )
    from plotting_functions import PlotFunctions

    T = [300.0, 310.0, 320.0]
    fluxes = {
        "r01.Flux": [1e-3, 2e-3, 4e-3],
        "r02.Flux": [-1e-3, -2e-3, -4e-3],
        "r03.Flux": [1e-3, -2e-3, 4e-3],
        "r04.Flux": [1e-3, 0.0, 4e-3],
    }
    rows = []
    for i, t in enumerate(T):
        sim_df = pd.DataFrame({"time": [1.0], **{k: [v[i]] for k, v in fluxes.items()}})
        rows += ReactionDataHandler.collect_log_values(
            sim_df, list(fluxes), 1.0, 0.1, t, "X"
        )
    df = pd.DataFrame(rows)

    out = ReactionParameterCalculator.calculate_reaction_parameters(
        df, T, [0.1], "X", "ri", reactions=["fwd", "rev", "flip", "zero"]
    )
    assert list(out["name"]) == ["r01.Flux", "r02.Flux"]  # flip and zero skipped
    fwd, rev = out.iloc[0], out.iloc[1]
    assert (fwd["sign"], rev["sign"]) == (1, -1)
    assert np.isclose(fwd["Ea"], rev["Ea"]) and np.isclose(fwd["slope"], rev["slope"])
    assert np.isfinite(fwd["R2"])

    flat = np.array([1.0, 1.0, 1.0])
    assert np.isnan(PlotFunctions.calculate_R2_values(np.array(T), flat, 0.0, 1.0))


if __name__ == "__main__":
    test_flux_sign_handling()
    print("ok")
