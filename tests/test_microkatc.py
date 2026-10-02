"""The microkatc.py command on the typed-energy example (real COPASI, runs in a few seconds).

python tests/test_microkatc.py
"""

import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import microkatc

EXAMPLE = ROOT / "examples" / "typed_energies" / "study.yaml"


def copy_example():
    folder = Path(
        tempfile.mkdtemp(prefix="my study ")
    )  # a space in the path on purpose
    shutil.copy(EXAMPLE, folder / "study.yaml")
    return folder / "study.yaml"


def run(*argv):
    """microkatc.py in its own process, as users run it (one study per process, see microkatc.py)"""
    result = subprocess.run(
        [sys.executable, str(ROOT / "microkatc.py"), *argv],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.returncode, result.stdout, result.stderr


def run_here(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = microkatc.main(list(argv))
    return code, out.getvalue(), err.getvalue()


def test_check_reports_ok_and_errors():
    study = copy_example()
    code, out, _ = run("check", str(study))
    assert code == 0 and "OK (2 steps, 4 species)" in out, out
    study.write_text(study.read_text().replace("product: P", "product: Q"))
    code, _, err = run("check", str(study))
    assert code == 1 and "species.product: Q does not appear in any step" in err, err
    assert not (study.parent / "results").exists()  # check never writes results


def test_run_from_other_folder():
    study = copy_example()
    started_in = tempfile.mkdtemp()
    result = subprocess.run(
        [sys.executable, str(ROOT / "microkatc.py"), "run", str(study)],
        cwd=started_in,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    results = study.parent / "results"
    assert sorted(p.name for p in results.iterdir()) == [
        "G_values_of_compounds",
        "G_values_of_reactions",
        "microkinetics.json",
        "microkinetics_simulations",
        "microkinetics_simulations_images",
        "model.cps",  # COPASI model of the last simulation, as before
        "run_info.json",
    ], sorted(p.name for p in results.iterdir())
    assert os.listdir(started_in) == []  # nothing written where the command was started
    summary = json.loads((results / "microkinetics.json").read_text())
    assert len(summary["t99_h"]) == 3 and all(
        abs(sum(x) - 1e-3) < 1e-9 for x in zip(*summary["catalyst_M"].values())
    )
    info = json.loads((results / "run_info.json").read_text())
    for key in (
        "study_sha256",
        "microkatc_commit",
        "python",
        "python-copasi",
        "numpy",
        "scipy",
        "started",
        "finished",
        "duration_s",
    ):
        assert key in info, key
    assert info["finished"] is not None


def test_rerun_same_study_is_allowed():
    study = copy_example()
    assert run("run", str(study))[0] == 0
    assert run("run", str(study))[0] == 0


def test_edited_study_stops_unless_fresh():
    study = copy_example()
    assert run("run", str(study))[0] == 0
    study.write_text(study.read_text().replace("TS1: 15.0", "TS1: 16.0"))
    code, _, err = run("run", str(study))
    assert (
        code == 1
        and "was made from a different version of study.yaml; delete it or run with --fresh"
        in err
    ), err
    assert run("run", str(study), "--fresh")[0] == 0


def test_failure_during_the_run_exits_2():
    study = copy_example()
    real = microkatc.run_analyses

    def broken(_study):
        raise RuntimeError("COPASI exploded")

    microkatc.run_analyses = broken
    try:
        code, _, err = run_here("run", str(study))
    finally:
        microkatc.run_analyses = real
    assert code == 2 and "run failed: COPASI exploded" in err, err


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
