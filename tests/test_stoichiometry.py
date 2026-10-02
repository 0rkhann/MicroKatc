"""COPASI follows the rate law and mass balance of a step with a coefficient (real COPASI).

For 2 A -> B (reverse barrier so high that the step is irreversible), mass action gives
rate = k[A]^2 with two A consumed per event, so A(t) = A0 / (1 + 2 k A0 t) and A + 2 B = A0.
Both spellings MicroKatc can hand to copasi_helper ("A + A = B" and "2*A = B") must follow it.
    python tests/test_stoichiometry.py
"""

import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

T, A0, GDIR, GINV, END = 350.0, 0.1, 20.0, 60.0, 2000


def simulate(equation):
    import copasi_parser as cpx

    os.chdir(tempfile.mkdtemp())
    pd.DataFrame(
        {"Rx": [equation], "TS": ["-"], "Gdir": [GDIR], "Ginv": [GINV]}
    ).to_csv("rx.csv", index=False)
    model = cpx.prepare_copasi_model(
        reactions_file="rx.csv",
        temp=T,
        initial_concentrations={"A": A0, "B": 0},
        csv_delim=",",
    )
    cpx.time_course_simulation(model, END, 1, outfile="sim.txt")
    model.self_destruct()
    return cpx.read_simulation("sim.txt"), cpx.k_calc_Eyring(GDIR, T)


def test_second_order_rate_law_and_mass_balance():
    for equation in ("A + A = B", "2*A = B"):
        df, k = simulate(equation)
        t, a, b = df["time"].to_numpy(), df["A"].to_numpy(), df["B"].to_numpy()
        analytic = A0 / (1 + 2 * k * A0 * t)
        assert np.abs(a - analytic).max() / A0 < 1e-5, (
            equation,
            np.abs(a - analytic).max() / A0,
        )
        # COPASI writes 6 significant digits, so A + 2 B rounds to within 1e-7 M (measured 1.0e-7)
        assert np.abs(a + 2 * b - A0).max() < 5e-7, equation
        first_order = A0 * np.exp(
            -k * t
        )  # the law a missing coefficient would give: must differ
        assert np.abs(a - first_order).max() / A0 > 0.1, equation


if __name__ == "__main__":
    test_second_order_rate_law_and_mass_balance()
    print("ok")
