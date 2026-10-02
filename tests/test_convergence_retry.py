"""COPASI retry path: a run that stops early is re-simulated and time is converted to hours once.

Uses a fake copasi_parser, so it runs without COPASI:
    python tests/test_convergence_retry.py   (or: pytest tests)
"""

import os
import sys
import tempfile
import types

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def test_retry_after_early_stop():
    calls = []

    class FakeModel:
        def self_destruct(self):
            pass

    def time_course_simulation(ch, total_time, outfile, **kwargs):
        calls.append(outfile)
        end = total_time if len(calls) > 1 else total_time / 2  # first run stops early
        pd.DataFrame({"time": [0.0, end], "prod": [0.0, 1.0]}).to_csv(
            outfile, index=False
        )
        return None, True

    cpx = types.ModuleType("copasi_parser")
    cpx.prepare_copasi_model = lambda **kwargs: FakeModel()
    cpx.time_course_simulation = time_course_simulation
    cpx.read_simulation = pd.read_csv
    sys.modules["copasi_parser"] = cpx

    os.chdir(tempfile.mkdtemp())
    from microkinetics_simulation import SimulationHandler

    handler = SimulationHandler(350.0, 28.7, "PMe3", 7200)
    df = handler.get_simulation_df({"PMe3": 1e-5, "CO": 1})
    assert len(calls) == 2
    assert df["time"].iloc[-1] == 2.0  # 7200 s -> 2 h

    cached = handler.get_simulation_df({"PMe3": 1e-5, "CO": 1})
    assert len(calls) == 2
    assert cached["time"].iloc[-1] == 2.0


if __name__ == "__main__":
    test_retry_after_early_stop()
    print("ok")
